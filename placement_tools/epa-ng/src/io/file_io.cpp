#include "io/file_io.hpp"

#include <stdexcept>
#include <fstream>
#include <functional>

#include "core/pll/pllhead.hpp"
#include "core/pll/pll_util.hpp"
#include "io/msa_reader.hpp"
#include "seq/MSA.hpp"
#include "util/constants.hpp"
#include "util/logging.hpp"


typedef struct fasta_record_s {
  pll_fasta_t* file;
  char * sequence = NULL;
  char * header = NULL;
} fasta_record_t;

int pll_fasta_fseek(pll_fasta_t* fd, const long int offset, const int whence)
{
  auto status = fseek(fd->fp, offset, whence);

  /* reset stripped char frequencies */
  fd->stripped_count = 0;
  for(size_t i=0; i<256; i++) {
    fd->stripped[i] = 0;
  }

  fd->line[0] = 0;
  if (!fgets(fd->line, PLL_LINEALLOC, fd->fp)) {
    pll_errno = PLL_ERROR_FILE_SEEK;
    snprintf(pll_errmsg, 200, "Unable to rewind and cache data");
    return -1;
  }
  fd->lineno = 1;

  return status;
}

MSA build_MSA_from_file(const std::string& msa_file,
                        const MSA_Info& info,
                        const bool premasking)
{
  MSA msa;
  auto reader = make_msa_reader(msa_file, info, premasking); 
  reader->read_next(msa, std::numeric_limits<size_t>::max());

  return msa;
}

pll_utree_s * build_tree_from_file(const std::string& tree_file, Tree_Numbers& nums)
{
  pll_utree_t * tree;
  pll_rtree_t * rtree;

  // load the tree unrooted
  if (!(rtree = pll_rtree_parse_newick(tree_file.c_str()))) {
    if (!(tree = pll_utree_parse_newick(tree_file.c_str()))) {
      throw std::runtime_error{std::string("Treeparsing failed! ") + pll_errmsg};
    }
  } else {
    tree = pll_rtree_unroot(rtree);
    pll_rtree_destroy(rtree, nullptr);

    /* optional step if using default PLL clv/pmatrix index assignments */
    pll_utree_reset_template_indices(get_root(tree), tree->tip_count);
  }

  if (tree->tip_count < 3) {
    throw std::runtime_error{"Number of tip nodes too small"};
  }

  nums = Tree_Numbers(tree->tip_count);

  set_missing_branch_lengths(tree, DEFAULT_BRANCH_LENGTH);

  return tree;
}

static unsigned int simd_autodetect()
{
  if (PLL_STAT(avx2_present))
    return PLL_ATTRIB_ARCH_AVX2;
  else if (PLL_STAT(avx_present))
    return PLL_ATTRIB_ARCH_AVX;
  else if (PLL_STAT(sse3_present))
    return PLL_ATTRIB_ARCH_SSE;
  else
    return PLL_ATTRIB_ARCH_CPU;
}

pll_partition_t *  build_partition_from_file( const raxml::Model& model, 
                                              Tree_Numbers& nums, 
                                              const int num_sites,
                                              const bool repeats)
{
  assert(nums.tip_nodes); // nums must have been initialized correctly

  auto attributes = simd_autodetect();

  if (repeats) {
    attributes |= PLL_ATTRIB_SITE_REPEATS;
  } else {
    attributes |= PLL_ATTRIB_PATTERN_TIP;
  }

  auto partition = pll_partition_create(nums.tip_nodes,
           nums.inner_nodes * 3, //number of extra clv buffers: 3 for every direction on the node
           model.num_states(),
           num_sites,
           1,
           nums.branches,
           model.num_ratecats(),
           (nums.inner_nodes * 3) + nums.tip_nodes, /* number of scaler buffers */
           attributes);

  if (!partition) {
    throw std::runtime_error{std::string("Could not create partition (build_partition_from_file). pll_errmsg: ") + pll_errmsg};
  }

  std::vector<double> rate_cats(model.num_ratecats(), 0.0);

  /* compute the discretized category rates from a gamma distribution
     with alpha shape */
  pll_compute_gamma_cats( model.alpha(),
                          model.num_ratecats(),
                          &rate_cats[0],
                          PLL_GAMMA_RATES_MEAN);
  pll_set_frequencies(partition,
                      0,
                      &(model.base_freqs(0)[0]));
  pll_set_subst_params( partition,
                        0,
                        &(model.subst_rates(0)[0]));
  pll_set_category_rates( partition,
                          &rate_cats[0]);

  // if (repeats) {
  //   pll_resize_repeats_lookup(partition, ( REPEATS_LOOKUP_SIZE ) * 10);
  // }

  return partition;

}

void file_check(const std::string& file_path)
{
  std::ifstream file(file_path.c_str());
  if (!file.good()) {
    throw std::runtime_error{std::string("file_check failed: ") + file_path};
  }

  file.close();
}
