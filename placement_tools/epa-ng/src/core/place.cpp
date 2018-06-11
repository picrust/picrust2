#include "core/place.hpp"

#include <fstream>
#include <string>
#include <memory>
#include <functional>
#include <limits>

#ifdef __OMP
#include <omp.h>
#endif

#include "io/file_io.hpp"
#include "io/jplace_util.hpp"
#include "io/msa_reader.hpp"
#include "io/Binary_Fasta.hpp"
#include "io/jplace_writer.hpp"
#include "util/stringify.hpp"
#include "util/logging.hpp"
#include "util/Timer.hpp"
#include "tree/Tiny_Tree.hpp"
#include "net/mpihead.hpp"
#include "pipeline/schedule.hpp"
#include "pipeline/Pipeline.hpp"
#include "seq/MSA.hpp"
#include "core/pll/pll_util.hpp"
#include "core/pll/epa_pll_util.hpp"
#include "core/Work.hpp"
#include "core/Lookup_Store.hpp"
#include "core/Work.hpp"
#include "core/heuristics.hpp"
#include "sample/Sample.hpp"
#include "set_manipulators.hpp"

#ifdef __MPI
#include "net/epa_mpi_util.hpp"
#endif

using mytimer = Timer<std::chrono::milliseconds>;

template <class T>
static void place(MSA& msa,
                  Tree& reference_tree,
                  const std::vector<pll_unode_t *>& branches,
                  Sample<T>& sample,
                  const Options& options,
                  std::shared_ptr<Lookup_Store>& lookup_store,
                  mytimer* time=nullptr)
{

#ifdef __OMP
  const unsigned int num_threads  = options.num_threads
                                  ? options.num_threads
                                  : omp_get_max_threads();
  omp_set_num_threads(num_threads);
  LOG_DBG << "Using threads: " << num_threads;
  LOG_DBG << "Max threads: " << omp_get_max_threads();
#else
  const unsigned int num_threads = 1;
#endif

  const size_t num_sequences  = msa.size();
  const size_t num_branches   = branches.size();

  std::vector<std::unique_ptr<Tiny_Tree>> branch_ptrs(num_threads);
  auto prev_branch_id = std::numeric_limits<size_t>::max();
  
  if (time){
    time->start();
  }
#ifdef __OMP
  #pragma omp parallel for schedule(guided, 10000), firstprivate(prev_branch_id)
#endif
  for (size_t i = 0; i < num_sequences * num_branches; ++i) {

#ifdef __OMP
    const auto tid = omp_get_thread_num();
#else
    const auto tid = 0;
#endif
    // reference to the threadlocal branch
    auto& branch = branch_ptrs[tid];

    const auto branch_id = static_cast<size_t>(i) / num_sequences;
    const auto seq_id = i % num_sequences;

    // get a tiny tree representing the current branch,
    // IF the branch has changed. Overwriting the old variable ensures
    // the now unused previous tiny tree is deallocated
    if ((branch_id != prev_branch_id) or not branch) {
      // as make_unique produces an rvalue, this is a move assignment and thus legal
      branch = std::make_unique<Tiny_Tree>(branches[branch_id],
                                           branch_id,
                                           reference_tree,
                                           false,
                                           options,
                                           lookup_store);
    }

    sample[seq_id][branch_id] = branch->place(msa[seq_id]);

    prev_branch_id = branch_id;
  }
  if (time){
    time->stop();
  }
}

template <class T>
static void place_thorough(const Work& to_place,
                  MSA& msa,
                  Tree& reference_tree,
                  const std::vector<pll_unode_t *>& branches,
                  Sample<T>& sample,
                  const Options& options,
                  std::shared_ptr<Lookup_Store>& lookup_store,
                  const size_t seq_id_offset=0,
                  mytimer* time=nullptr)
{
  
#ifdef __OMP
  const unsigned int num_threads  = options.num_threads
                                  ? options.num_threads
                                  : omp_get_max_threads();
  omp_set_num_threads(num_threads);
  LOG_DBG << "Using threads: " << num_threads;
  LOG_DBG << "Max threads: " << omp_get_max_threads();
#else
  const unsigned int num_threads = 1;
#endif

  // split the sample structure such that the parts are thread-local
  std::vector<Sample<T>> sample_parts(num_threads);

  // build vector of elements
  std::vector<Work::Work_Pair> id;
  for(auto it = to_place.begin(); it != to_place.end(); ++it) {
    id.push_back(*it);
  }

  // Map from sequence indices to indices in the pquery vector.
  auto seq_lookup_vec = std::vector<std::unordered_map<size_t, size_t>>(num_threads);

  std::vector<std::unique_ptr<Tiny_Tree>> branch_ptrs(num_threads);
  auto prev_branch_id = std::numeric_limits<size_t>::max();
  
  // work seperately
  if (time){
    time->start();
  }
#ifdef __OMP
  #pragma omp parallel for schedule(dynamic), firstprivate(prev_branch_id)
#endif
  for (size_t i = 0; i < id.size(); ++i) {

#ifdef __OMP
    const auto tid = omp_get_thread_num();
#else
    const auto tid = 0;
#endif
    auto& local_sample = sample_parts[tid];
    auto& seq_lookup = seq_lookup_vec[tid];

    const auto branch_id = id[i].branch_id;
    const auto seq_id = id[i].sequence_id;
    const auto& seq = msa[seq_id];

    // get a tiny tree representing the current branch,
    // IF the branch has changed. Overwriting the old variable ensures
    // the now unused previous tiny tree is deallocated
    if ((branch_id != prev_branch_id) or not branch_ptrs[tid]) {
      // as make_unique produces an rvalue, this is a move assignment and thus legal
      branch_ptrs[tid] = std::make_unique<Tiny_Tree>(branches[branch_id],
                                           branch_id,
                                           reference_tree,
                                           true,
                                           options,
                                           lookup_store);
    }

    if (seq_lookup.count( seq_id ) == 0) {
      auto const new_idx = local_sample.add_pquery( seq_id_offset + seq_id, seq.header() );
      seq_lookup[ seq_id ] = new_idx;
    }
    assert( seq_lookup.count( seq_id ) > 0 );
    local_sample[ seq_lookup[ seq_id ] ].emplace_back( branch_ptrs[tid]->place(seq) );

    prev_branch_id = branch_id;
  }
  if (time){
    time->stop();
  }
  // merge samples back
  merge(sample, std::move(sample_parts));
  collapse(sample);
}

void simple_mpi(Tree& reference_tree,
                const std::string& query_file,
                const MSA_Info& msa_info,
                const std::string& outdir,
                const Options& options,
                const std::string& invocation)
{
  const auto num_branches = reference_tree.nums().branches;

  // get all edges
  std::vector<pll_unode_t *> branches(num_branches);
  auto num_traversed_branches = utree_query_branches(reference_tree.tree(), &branches[0]);
  if (num_traversed_branches != num_branches) {
    throw std::runtime_error{"Traversing the utree went wrong during pipeline startup!"};
  }

  auto lookups = 
    std::make_shared<Lookup_Store>(num_branches, reference_tree.partition()->states);

  // some MPI prep
  int local_rank = 0;
  int num_ranks = 1;

  MPI_COMM_RANK(MPI_COMM_WORLD, &local_rank);
  MPI_COMM_SIZE(MPI_COMM_WORLD, &num_ranks);

  LOG_INFO << "Number of ranks: " << num_ranks;

  auto reader = make_msa_reader(query_file, msa_info, options.premasking);

  size_t local_rank_seq_offset = 0;

  if (num_ranks > 1) {
    // how many should each rank read?
    const size_t part_size = ceil(reader->num_sequences() / static_cast<double>(num_ranks));
    LOG_INFO << "Number of sequences per rank: " << part_size;

    // read only the locally relevant part of the queries
    // ... by skipping the appropriate amount
    local_rank_seq_offset = part_size * local_rank;
    reader->skip_to_sequence( local_rank_seq_offset );
    // and limiting the reading to the given window
    reader->constrain(part_size);
  }

  size_t num_sequences = 0;
  Work all_work(std::make_pair(0, num_branches), std::make_pair(0, options.chunk_size));

  Work blo_work;

  size_t chunk_num = 1;

  using Sample = Sample<Placement>;
  std::future<void> prev_gather;
  MSA chunk;
  size_t sequences_done = 0; // not just for info output!

  // prepare output file
  LOG_INFO << "Output file: " << outdir + "epa_result.jplace";
  jplace_writer jplace( outdir, "epa_result.jplace",
                        get_numbered_newick_string(reference_tree.tree()),
                        invocation);


  Sample preplace(options.chunk_size, num_branches);

  while ( (num_sequences = reader->read_next(chunk, options.chunk_size)) ) {

    assert(chunk.size() == num_sequences);

    LOG_DBG << "num_sequences: " << num_sequences << std::endl;

    const size_t seq_id_offset = sequences_done + local_rank_seq_offset;

    if (num_sequences < options.chunk_size) {
      all_work = Work(std::make_pair(0, num_branches), std::make_pair(0, num_sequences));
      preplace = Sample(num_sequences, num_branches);
    }

    if (options.prescoring) {

      LOG_DBG << "Preplacement." << std::endl;
      place(chunk,
            reference_tree,
            branches,
            preplace,
            options,
            lookups);

      LOG_DBG << "Selecting candidates." << std::endl;

      blo_work = apply_heuristic(preplace, options);

    } else {
      blo_work = all_work;
    }

    Sample blo_sample;

    LOG_DBG << "BLO Placement." << std::endl;
    place_thorough( blo_work,
                    chunk,
                    reference_tree,
                    branches,
                    blo_sample,
                    options,
                    lookups,
                    seq_id_offset);

    // Output
    compute_and_set_lwr(blo_sample);
    filter(blo_sample, options);

    // pass the result chunk to the writer
    jplace.write( blo_sample );

    sequences_done += num_sequences;
    LOG_INFO << sequences_done  << " Sequences done!";
    ++chunk_num;
  }

  jplace.wait();

  MPI_BARRIER(MPI_COMM_WORLD);
}

