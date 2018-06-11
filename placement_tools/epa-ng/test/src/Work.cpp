#include "Epatest.hpp"

#include "core/Work.hpp"
#include "sample/Sample.hpp"

using namespace std;

TEST(Work, create_from_range)
{
  size_t upper_branch = 10;
  size_t upper_sequences = 12;
  Work work(make_pair(0,10), make_pair(0,12));

  // printf("\nWork");
  // for (auto i = work.begin(); i != work.end(); ++i)
  // {
  //   printf("\nbranch %d: ", i->first);
  //   for (auto& seq_id : i->second)
  //   {
  //     printf(" %d ", seq_id);
  //   }
  // }
  // printf("\n");

  EXPECT_EQ( upper_branch * upper_sequences, work.size() );
}
