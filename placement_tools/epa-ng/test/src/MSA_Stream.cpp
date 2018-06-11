#include "Epatest.hpp"

#include "seq/MSA_Stream.hpp"
#include "seq/MSA.hpp"
#include "seq/MSA_Info.hpp"
#include "io/file_io.hpp"

#include <string>

using namespace std;

TEST(MSA_Stream, reading)
{
  MSA_Info info(env->combined_file);
  MSA complete_msa = build_MSA_from_file(env->combined_file, info);
  const size_t chunk_size = 3;
  MSA read_msa;
  MSA_Stream streamed_msa(env->combined_file, info, false);

  for (size_t i = 0; i < complete_msa.size(); i++)
  {
    if ((i % chunk_size) == 0)
    {
      streamed_msa.read_next(read_msa, chunk_size);
    }
    EXPECT_EQ(complete_msa[i], read_msa[i % chunk_size]);
  }
  MSA_Stream dummy;
}

TEST(MSA_Stream, reading_masked)
{
  MSA_Info info(env->combined_file);
  MSA complete_msa = build_MSA_from_file(env->combined_file, info, true);
  const auto chunk_size = 3;
  MSA read_msa;
  MSA_Stream streamed_msa(env->combined_file, info, true);

  for (size_t i = 0; i < complete_msa.size(); i++)
  {
    if ((i % chunk_size) == 0)
    {
      streamed_msa.read_next(read_msa, chunk_size);
    }
    EXPECT_EQ(complete_msa[i], read_msa[i % chunk_size]);
  }
  MSA_Stream dummy;
}