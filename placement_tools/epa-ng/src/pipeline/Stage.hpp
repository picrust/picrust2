#pragma once

#include <functional>
#include <type_traits>
#include <memory>

#include "pipeline/Token.hpp"
#include "util/function_traits.hpp"
#include "net/epa_mpi_util.hpp"


/**
 * Interface for Stages
 */
class Stage
{

public:
  Stage() = default;
  ~Stage() = default;

  virtual bool exec() const final
  {
    return exec_;
  }

  virtual void exec(const bool b) final
  {
    exec_ = b;
  }

private:
  bool exec_ = true;
  
};

/**
 * Templated interface for Stages
 */
template <size_t I, typename lambda>
class Typed_Stage : public Stage
{
public:
  static constexpr size_t id() {return I;}

  using traits = thrill::common::FunctionTraits<lambda>;

  using base_in_type  = typename traits::template arg<0>;
  using base_out_type = typename traits::result_type;
  using in_type   = std::remove_reference_t<base_in_type>;
  using out_type  = std::remove_reference_t<base_out_type>;
  // using in_type   = std::add_lvalue_reference_t<in_type_noref>;
  // using out_type  = std::add_lvalue_reference_t<out_type_noref>;

  using accept_func_type  = std::function<void(in_type&)>;
  using process_func_type = std::function<out_type(in_type&)>;
  using put_func_type     = std::function<void(out_type&)>;

  /**
   * Constructs a Stage with the three supplied functions
   */
  Typed_Stage(accept_func_type   accept,
              process_func_type  process,
              put_func_type      put)
    : accept_(accept)
    , process_(process)
    , put_(put)
  { }

  /**
   * Constructs a Stage with the supplied core process function. 
   * accept and put are defaulted to noops.
   */
  Typed_Stage(process_func_type process)
    : accept_([](in_type) { })
    , process_(process)
    , put_([](out_type) { })
  { }

  Typed_Stage()		= delete;
  ~Typed_Stage() 	= default;

  inline void accept(in_type& arg) const
  {
    accept_(arg);
  }

  inline out_type process(in_type& arg) const
  {
    return process_(arg);
  }

  inline void put(out_type& arg) const
  {
    put_(arg);
  }

  void set_accept_func(accept_func_type& f)
  {
    accept_ = f;
  }

  void set_put_func(put_func_type& f)
  {
    put_ = f;
  }
  

private:

  accept_func_type  accept_;
  process_func_type process_;
  put_func_type     put_;

};
