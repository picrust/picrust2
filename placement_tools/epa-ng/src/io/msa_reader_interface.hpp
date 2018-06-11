#pragma once

#include "seq/MSA.hpp"

class msa_reader
{

public:
  msa_reader() = default;
  virtual ~msa_reader() = default;

  virtual void constrain(const size_t max_read) = 0;
  virtual void skip_to_sequence(const size_t n) = 0;
  virtual size_t num_sequences() = 0;
  virtual size_t read_next(MSA& result, const size_t number) = 0;

};
