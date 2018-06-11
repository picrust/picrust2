#include "io/Binary.hpp"

#include <stdexcept>
#include <algorithm>

#include "util/constants.hpp"
#include "tree/Tree.hpp"

int safe_fclose(FILE* fptr) { return fptr ? fclose(fptr) : 0; }

Binary::Binary(Binary && other) 
  : bin_fptr_(nullptr, safe_fclose)
{
  std::swap(bin_fptr_, other.bin_fptr_);
  std::swap(map_, other.map_);
}

Binary& Binary::operator=(Binary && other)
{
  bin_fptr_ = std::move(other.bin_fptr_);
  map_ = std::move(other.map_);
  return *this;
}

Binary::Binary(const std::string& binary_file_path) 
  : bin_fptr_(nullptr, safe_fclose)
{
  // open the binary file
  pll_binary_header_t header;
  bin_fptr_ = file_ptr_type(pllmod_binary_open(binary_file_path.c_str(), &header), safe_fclose);

  if (!bin_fptr_) {
    throw std::runtime_error{"Could not open binary file for reading."};
  }

  if (header.access_type != PLLMOD_BIN_ACCESS_RANDOM) {
    throw std::runtime_error{"Binary file must be random access enabled."};
  }

  if (header.n_blocks <= 0) {
    throw std::runtime_error{std::string("Binary file header must have nonzero positive number of blocks: ")
              + std::to_string(header.n_blocks)};
  }

  // proccess the random access map
  unsigned int n_blocks;
  pll_block_map_t* block_map;
  block_map = pllmod_binary_get_map(bin_fptr_.get(), &n_blocks);

  assert(block_map);
  assert(n_blocks);

  for (size_t i = 0; i < n_blocks; i++) {
    map_.push_back(block_map[i]);
  }

  free(block_map);
}

static long int get_offset(std::vector<pll_block_map_t>& map, const int block_id)
{
  auto item = map.begin();
  while (item != map.end()) {
    if (item->block_id == block_id) {
      break;
    }
    item++;
  }

  if(item == map.end()) {
    throw std::runtime_error{std::string("Map does not contain block_id: ") + std::to_string(block_id)};
  }
  return item->block_offset;
}

void Binary::load_clv(pll_partition_t * partition,
                      const unsigned int clv_index)
{
  assert(bin_fptr_);
  assert(clv_index < partition->clv_buffers + partition->tips);
  if (partition->attributes & PLL_ATTRIB_PATTERN_TIP) {
    assert(clv_index >= partition->tips);
  }

  if (!(partition->clv[clv_index])) {
    const size_t clv_size = pll_get_clv_size(partition, clv_index) * sizeof(double);

    partition->clv[clv_index] = static_cast<double*>(pll_aligned_alloc(clv_size, partition->alignment));
    if (!partition->clv[clv_index]) {
      throw std::runtime_error{"Could not allocate CLV memory"};
    }
  }

  {
    unsigned int attributes;
    std::lock_guard<std::mutex> lock(file_mutex_);
    auto err = pllmod_binary_clv_load(bin_fptr_.get(),
                                    0,
                                    partition,
                                    clv_index,
                                    &attributes,
                                    get_offset(map_, clv_index));
    if (err != PLL_SUCCESS) {
      throw std::runtime_error{std::string("Loading CLV failed: ") 
                              + pll_errmsg 
                              + std::string(". CLV index: ") 
                              + std::to_string(clv_index)};
    }
  }
}

void Binary::load_tipchars( pll_partition_t * partition,
                            const unsigned int tipchars_index)
{
  assert(bin_fptr_);
  assert(tipchars_index < partition->tips);
  assert(partition->attributes & PLL_ATTRIB_PATTERN_TIP);

  unsigned int type = 0;
  unsigned int attributes = 0;
  size_t size = 0;

  {
    std::lock_guard<std::mutex> lock(file_mutex_);
    auto ptr = pllmod_binary_custom_load( bin_fptr_.get(),
                                          0,
                                          &size,
                                          &type,
                                          &attributes,
                                          get_offset(map_, tipchars_index));
    if (!ptr) {
      throw std::runtime_error{std::string("Loading tipchar failed: ") + pll_errmsg};
    }
    partition->tipchars[tipchars_index] = static_cast<unsigned char*>(ptr);
  }
}

void Binary::load_scaler( pll_partition_t * partition, 
                          const unsigned int scaler_index)
{
  assert(bin_fptr_);
  assert(scaler_index < partition->scale_buffers);

  auto block_offset = partition->clv_buffers + partition->tips;

  unsigned int type, attributes;
  size_t size;

  {
    std::lock_guard<std::mutex> lock(file_mutex_);
    auto ptr = pllmod_binary_custom_load( bin_fptr_.get(), 
                                          0,
                                          &size, 
                                          &type, 
                                          &attributes, 
                                          get_offset(map_, block_offset + scaler_index));
    if (!ptr) {
      throw std::runtime_error{std::string("Loading scaler failed: ") + pll_errmsg};
    }
  
    partition->scale_buffer[scaler_index] = static_cast<unsigned int*>(ptr);
  }
}

pll_partition_t* Binary::load_partition()
{
  std::lock_guard<std::mutex> lock(file_mutex_);
  // make skeleton partition that only allocates the pointers to the clv/tipchar buffers
  unsigned int part_attribs = PLLMOD_BIN_ATTRIB_PARTITION_LOAD_SKELETON;
  auto partition =  pllmod_binary_partition_load( bin_fptr_.get(), 
                                                  0, 
                                                  nullptr, 
                                                  &part_attribs, 
                                                  get_offset(map_, -1));

  if (!partition) {
    throw std::runtime_error{std::string("Error loading partition: ") + pll_errmsg};
  }

  if (partition->attributes & PLL_ATTRIB_SITE_REPEATS) {
    unsigned int repeats_attribs = 0;
    if (pllmod_binary_repeats_load( bin_fptr_.get(), 
                                    0, 
                                    partition, 
                                    &repeats_attribs, 
                                    get_offset(map_, -3))
        != PLL_SUCCESS) {
      throw std::runtime_error{std::string("Error loading repeats: ") + pll_errmsg};   
    }
  }

  return partition;
}

pll_utree_t* Binary::load_utree(const unsigned int num_tips)
{
  std::lock_guard<std::mutex> lock(file_mutex_);
  unsigned int attributes = 0;
  auto root =  pllmod_binary_utree_load(bin_fptr_.get(), 
                                        0, 
                                        &attributes, 
                                        get_offset(map_, -2));
  if (!root) {
    throw std::runtime_error{std::string("Loading tree: ") + pll_errmsg};
  }

  return pll_utree_wraptree(root, num_tips);
}

static int full_trav(pll_unode_t*)
{
  return 1;
}

static auto create_scaler_to_clv_map(Tree& tree)
{
  const auto num_scalers = tree.partition()->scale_buffers;
  std::vector<unsigned int> map(num_scalers);

  std::vector<pll_unode_t*> travbuffer(tree.nums().nodes);
  unsigned int trav_size = 0;
  pll_utree_traverse( get_root(tree.tree()),
                      PLL_TREE_TRAVERSE_POSTORDER,
                      full_trav, 
                      &travbuffer[0], 
                      &trav_size);

  for (auto& n : travbuffer) {
    if (n->scaler_index != PLL_SCALE_BUFFER_NONE) {
      map[n->scaler_index] = n->clv_index;
    }

    // node is an inner?
    if (n->next) {
      auto nn = n->next;
      if (nn->scaler_index != PLL_SCALE_BUFFER_NONE) {
        map[nn->scaler_index] = nn->clv_index;
      }
      nn = nn->next;
      if (nn->scaler_index != PLL_SCALE_BUFFER_NONE) {
        map[nn->scaler_index] = nn->clv_index;
      }
    }
  }

  return map;
}

/**
  Writes the structures and data encapsulated in Tree to the specified file in the binary format.
  Writes them in such a way that the Binary class can read them.
*/
void dump_to_binary(Tree& tree, const std::string& file)
{
  const auto num_clvs = tree.partition()->clv_buffers;
  const auto num_tips = tree.partition()->tips;
  const auto num_scalers = tree.partition()->scale_buffers;
  const auto max_clv_index = num_clvs + num_tips;
  
  const bool use_tipchars = tree.partition()->attributes & PLL_ATTRIB_PATTERN_TIP;
  const bool use_repeats = tree.partition()->attributes & PLL_ATTRIB_SITE_REPEATS;

  int block_id = use_repeats ? -3 : -2;

  const unsigned int num_blocks = abs(block_id) + num_clvs + num_tips + num_scalers;

  pll_binary_header_t header;
  auto fptr =  pllmod_binary_create(
    file.c_str(),
    &header,
    PLLMOD_BIN_ACCESS_RANDOM,
    num_blocks);

  if(!fptr) {
    throw std::runtime_error{std::string("Error opening binary file for writing: ") + pll_errmsg};
  }

  const auto attributes = PLLMOD_BIN_ATTRIB_UPDATE_MAP | PLLMOD_BIN_ATTRIB_PARTITION_DUMP_WGT;


  if (use_repeats and not
      pllmod_binary_repeats_dump(fptr, block_id++, tree.partition(), attributes)) {
    throw std::runtime_error{std::string("Error dumping the repeats: ") + pll_errmsg};
  }

  // dump the utree structure
  if(!pllmod_binary_utree_dump(fptr, block_id++, get_root(tree.tree()), num_tips, attributes)) {
    throw std::runtime_error{std::string("Error dumping the utree to binary: ") + pll_errmsg};
  }

  // dump the partition
  if(!pllmod_binary_partition_dump(fptr, block_id++, tree.partition(), attributes)) {
    throw std::runtime_error{std::string("Error dumping partition to binary: ") + pll_errmsg};
  }

  // dump the tipchars, but only if partition uses them
  size_t tip_index = 0;
  if (use_tipchars) {
    for (tip_index = 0; tip_index < num_tips; tip_index++) {
      if(!pllmod_binary_custom_dump(fptr, 
                                    block_id++, 
                                    tree.partition()->tipchars[tip_index],
                                    tree.partition()->sites * sizeof(unsigned char), 
                                    attributes)) {
        throw std::runtime_error{std::string("Error dumping tipchars to binary: ") + pll_errmsg};
      }
    }
  }

  // dump the clvs
  for (size_t clv_index = tip_index; clv_index < max_clv_index; clv_index++) {
    if(!pllmod_binary_clv_dump( fptr, 
                                block_id++, 
                                tree.partition(), 
                                clv_index, 
                                attributes)) {
      throw std::runtime_error{std::string("Error dumping clvs to binary: ") + pll_errmsg};
    }
  }

  const auto scaler_to_clv = create_scaler_to_clv_map(tree);

  const auto scaler_ptr = tree.partition()->scale_buffer;
  
  for (size_t scaler_index = 0; scaler_index < num_scalers; scaler_index++) {
    
    const auto scaler_size =  pll_get_sites_number( tree.partition(), 
                                                    scaler_to_clv[scaler_index]);

    // with the repeats the scale buffers might not be allocated. dirty fix: 
    // allocate them in this case, just so they can be written and later used
    if (scaler_ptr[scaler_index] == nullptr) {
      scaler_ptr[scaler_index] = static_cast<unsigned int *>(calloc(scaler_size, sizeof(unsigned int)));
    }

    if(!pllmod_binary_custom_dump(fptr, 
                                  block_id++, 
                                  scaler_ptr[scaler_index],
                                  scaler_size * sizeof(unsigned int), 
                                  attributes)) {
      throw std::runtime_error{std::string("Error dumping scalers to binary: ") + pll_errmsg};
    }
  }

  fclose(fptr);
}
