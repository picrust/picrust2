#pragma once

#include <gtest/gtest.h>

#include "core/raxml/Model.hpp"
#include "util/Options.hpp"

// The testing environment
class Epatest : public ::testing::Environment {
public:
  // You can remove any or all of the following functions if its body
  // is empty.

  Epatest() {};

  virtual ~Epatest() {
    // You can do clean-up work that doesn't throw exceptions here.
  }

  // If the constructor and destructor are not enough for setting up
  // and cleaning up each test, you can define the following methods:

  virtual void SetUp() {
    // Code here will be called immediately after the constructor (right
    // before each test).
  }

  virtual void TearDown() {
    // Code here will be called immediately after each test (right
    // before the destructor).
  }

  // Objects declared here can be used by all tests in the test case for Foo.
  std::string data_dir;
  std::string out_dir;
  std::string tree_file;
  std::string tree_file_rooted;
  std::string reference_file;
  std::string combined_file;
  std::string query_file;
  std::string binary_file;
  std::string info_file;
  raxml::Model model = raxml::Model("GTR+G");
  Options options;

};

extern Epatest* env;

#include <cmath>
// #include <cmath>

#define COMPL_REPEATS       (1 << 0)
#define COMPL_OPTIMIZE      (1 << 1)
#define COMPL_SLIDING_BLO   (1 << 2)
#define COMPL_PRESCORING    (1 << 3)
#define COMPL_MASKING       (1 << 4)

static inline Options get_options_config(const unsigned int d)
{
  Options o;
  if (d & COMPL_REPEATS) {
    o.repeats = not o.repeats;
  }
  if (d & COMPL_OPTIMIZE) {
    o.opt_branches = not o.opt_branches;
    o.opt_model = not o.opt_model;
  }
  if (d & COMPL_SLIDING_BLO) {
    o.sliding_blo = not o.sliding_blo;
  }
  if (d & COMPL_PRESCORING) {
    o.prescoring = not o.prescoring;
  }
  if (d & COMPL_MASKING) {
    o.premasking = not o.premasking;
  }
  return o;
}

template <class Func>
void all_combinations(Func f, bool verbose=false)
{
  for (size_t i = 0; i < pow(2, 4); ++i) {
    auto o = get_options_config(i);
    if (verbose) {
      printf("\nrepeats\toptim\tsliding\tprescore\tmasking\n");
      printf( "%d\t%d\t%d\t%d\t%d\n", 
              o.repeats,
              o.opt_model,
              o.sliding_blo,
              o.prescoring,
              o.premasking);
    }
    f(o);
  }
}
