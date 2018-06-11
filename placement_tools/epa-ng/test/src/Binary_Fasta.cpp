#include "Epatest.hpp"

#include "io/Binary_Fasta.hpp"
#include "io/file_io.hpp"
#include "seq/MSA.hpp"

#include "genesis/utils/core/options.hpp"

static void compare_msas(const MSA& lhs, const MSA& rhs)
{
  ASSERT_EQ(lhs.size(), rhs.size());

  for (size_t i = 0; i < lhs.size(); ++i) {
    EXPECT_STREQ(lhs[i].header().c_str(), rhs[i].header().c_str());
    EXPECT_STREQ(lhs[i].sequence().c_str(), rhs[i].sequence().c_str());
  }
}

TEST(Binary_Fasta, 4bit_store_and_load)
{
  genesis::utils::Options::get().allow_file_overwriting(true);

  const std::string orig_file(env->combined_file);
  const std::string binfile_name(orig_file + ".bin");

  auto msa = build_MSA_from_file(orig_file, MSA_Info(env->combined_file));

  Binary_Fasta::save(msa, binfile_name);

  auto read_msa = Binary_Fasta::load(binfile_name, false);

  compare_msas(msa, read_msa);
}

TEST(Binary_Fasta, reader)
{
  genesis::utils::Options::get().allow_file_overwriting(true);

  const std::string orig_file(env->combined_file);
  const std::string binfile_name(orig_file + ".bin");

  MSA_Info info(env->combined_file);

  auto msa = build_MSA_from_file(orig_file, info);

  Binary_Fasta::save(msa, binfile_name);

  Binary_Fasta_Reader reader(binfile_name, info);

  const size_t skip = 3;

  reader.skip_to_sequence(skip);

  MSA read_msa;
  size_t i = skip;
  const size_t chunksize = 5;
  size_t num_sequences = 0;
  while ( (num_sequences = reader.read_next(read_msa, chunksize)) ) {

    ASSERT_EQ(num_sequences, read_msa.size()) << "bad size at i=" << i;
  
    for (size_t k = 0; k < num_sequences; ++k) {
      EXPECT_STREQ(msa[i+k].header().c_str(), read_msa[k].header().c_str());
      EXPECT_STREQ(msa[i+k].sequence().c_str(), read_msa[k].sequence().c_str());
    }

    i+=num_sequences;

  }
}
