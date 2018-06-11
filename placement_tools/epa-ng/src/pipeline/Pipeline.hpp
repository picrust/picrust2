#pragma once

#include <functional>
#include <vector>
#include <stdexcept>
#include <type_traits>
#include <tuple>
#include <memory>

#include "pipeline/Stage.hpp"
#include "pipeline/Token.hpp"
#include "util/Timer.hpp"
#include "net/Intercom.hpp"
#include "pipeline/schedule.hpp"
#include "util/function_traits.hpp"
#include "util/template_magic.hpp"

/**
 * Building a Stage Tuple out of a bunch of lambda functions/functors
 */
template < class I, class... lambdas>
struct stage_types_base;

template < std::size_t... I, class... lambdas >
struct stage_types_base<std::index_sequence<I...>, lambdas...>
{
  using types = typename std::tuple< Typed_Stage<I, lambdas>... >;
};

template < class... lambdas >
struct stage_types 
  : stage_types_base<std::make_index_sequence<sizeof...(lambdas) >, lambdas...>
{
};

/**
 * Basic Pipeline Class. Runs all stages in serial.
 */
template <class... lambdas>
class Pipeline
{
  using stack_type      = typename stage_types< lambdas... >::types;
  using token_set_type  = typename token_types< stack_type >::types;

public:
  using hook_type       = std::function<void()>;

  Pipeline( const stack_type& stages, 
            const hook_type& per_loop_hook,
            const hook_type& init_hook,
            const hook_type& final_hook)
    : stages_(stages)
    , per_loop_hook_(per_loop_hook)
    , init_hook_(init_hook)
    , final_hook_(final_hook)
    , icom_(std::tuple_size<stack_type>::value)
  { 
     init_pipeline_();
  }

  ~Pipeline() = default;

  template <class Function>
  auto push(const Function& f) const
  {
    constexpr size_t num_stages = sizeof...(lambdas) + 1u;
    constexpr auto new_stage_id = num_stages - 1u;

    using stage_type = Typed_Stage<new_stage_id, Function>;
    using new_stack_type = typename stage_types<lambdas..., Function>::types;

    new_stack_type stage_tuple
      = std::tuple_cat(stages_, std::make_tuple(stage_type(f)));
    
    return Pipeline<lambdas..., Function>(stage_tuple, per_loop_hook_, init_hook_, final_hook_);
  }

  void process()
  {
    token_set_type tokens;
    // "last" token that is still used on the particular MPI-Rank (or thread or...)
    Token const * last_token = nullptr;
    
    constexpr size_t dedicated_write = std::tuple_element_t<std::tuple_size<decltype(stages_)>::value - 1u
                                              , decltype(stages_)>::id();
    LOG_DBG1 << "dedicated_write: " << dedicated_write;

    if (icom_.stage_active(dedicated_write)) {
      init_hook_();
    }

    size_t chunk_num = 1;
    do { 
      elapsed_time_.start();

      // per-loop pre-hook
      per_loop_hook_();

      for_each(stages_, [&](const auto s) {

        if (s.exec()) {

          constexpr auto stage_id = s.id();

          auto& in_token = std::get<stage_id>(tokens);
          auto& out_token = std::get<stage_id+1u>(tokens);

          s.accept(in_token); // noop if shared mem, mpi_merge_receive if mpi

          if (in_token.valid()) {
            LOG_DBG1 << "in_token size: " << in_token.size();
            out_token = s.process(in_token); // do the actual work
          } else {
            LOG_DBG1 << std::to_string(icom_.rank()) << " received end token. Terminating.";
            out_token.is_last(true);
          }

          if (stage_id != 0) {
            // carry over the token status
            out_token.status(in_token.status());
          }

          LOG_DBG1 << "out_token size: " << out_token.size();

          s.put(out_token); // noop if shared mem, mpi_split_send if mpi

          last_token = &out_token;

        }
      });

      elapsed_time_.stop();

      // if(last_token->rebalance()) {
      if (rebalance_on_(chunk_num)) {
        // do the kansas city shuffle...
        // calculate new schedule
        icom_.rebalance(elapsed_time_);

        init_pipeline_();

        advance_rebalance_check_();
      }

      // clear the tokens at the end of each chunk
      for_each(tokens, [](auto& t) {
        t.clear();
      });

      ++chunk_num;
    } while (last_token->valid()); //returns valid if data token or default initialized

    if (icom_.stage_active(dedicated_write)) {
      final_hook_();
    }

    icom_.barrier();
  }

private:

  void init_pipeline_()
  {
    // reassign the local per-stage execution status
    assign_exec_status_();

    // reassign the inter-stage communication operations
    assign_comm_ops_();
  }

  bool rebalance_on_(const size_t chunk)
  {
    return chunk == next_rebalance_chunk_;
  }

  void advance_rebalance_check_()
  {
    rebalance_delta_ *= 2;
    next_rebalance_chunk_ += rebalance_delta_;
  }

  void assign_exec_status_()
  {
    for_each(stages_, [&](auto& s) {
      s.exec(icom_.stage_active(s.id()));
    });
  }

  /**
   * Assign appropriate communication operations between the stages according to the
   * current schedule.
   *
   */
  void assign_comm_ops_()
  {
    /**
     * Contiguous active stages communicate via no-op (aka. they share the same memory),
     * whereas between such clusters MPI communication routines are called
     */
    #ifdef __MPI
    for_each_pair(stages_, [&](auto& lhs, auto& rhs) {
      using lhs_t = std::remove_reference_t<decltype(lhs)>;
      using rhs_t = std::remove_reference_t<decltype(rhs)>;
      using put_arg_t     = typename lhs_t::out_type;
      using accept_arg_t  = typename rhs_t::in_type;
      using put_func_t    = typename lhs_t::put_func_type;
      using accept_func_t = typename rhs_t::accept_func_type;

      put_func_t put        = [](put_arg_t&){};
      accept_func_t accept  = [](accept_arg_t&){};

      using namespace std::placeholders;

      constexpr auto src = lhs.id();
      constexpr auto dst = rhs.id();
      
      if (lhs.exec()) {
        if (not rhs.exec()) {
          put = std::bind(
            epa_mpi_split_send<put_arg_t>, 
            _1, 
            std::ref(icom_.schedule(dst)), 
            MPI_COMM_WORLD,
            std::ref(icom_.previous_requests()),
            std::ref(elapsed_time_)
          );
        }
      } else {
        if (rhs.exec()) {
          accept = std::bind(
            epa_mpi_receive_merge<accept_arg_t>, 
            _1, 
            std::ref(icom_.schedule(src)),
            MPI_COMM_WORLD,
            // std::ref(icom_.previous_requests()),
            std::ref(elapsed_time_)
          );
        }
      }

      lhs.set_put_func(put);
      rhs.set_accept_func(accept); 
    });
    #endif //__MPI
  }

  stack_type stages_;
  hook_type per_loop_hook_;
  hook_type init_hook_;
  hook_type final_hook_;
  Intercom icom_;
  Timer<> elapsed_time_;

  size_t next_rebalance_chunk_ = 3;
  size_t rebalance_delta_ = next_rebalance_chunk_;

};

template <class stage_f>
auto make_pipeline( const stage_f& first_stage, 
                    const typename Pipeline<stage_f>::hook_type& per_loop_hook,
                    const typename Pipeline<stage_f>::hook_type& init_hook,
                    const typename Pipeline<stage_f>::hook_type& final_hook) 
{
  return Pipeline<stage_f>( std::make_tuple(Typed_Stage<0u, stage_f>(first_stage)), 
                            per_loop_hook, 
                            init_hook, 
                            final_hook);
}
