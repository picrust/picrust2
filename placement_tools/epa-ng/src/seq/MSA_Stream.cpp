#include "seq/MSA_Stream.hpp"

#include <chrono>

#include "util/logging.hpp"

static void read_chunk( MSA_Stream::file_type& iter,
                        const MSA_Info& info,
                        const bool premasking,
                        const size_t number, 
                        MSA_Stream::container_type& prefetch_buffer,
                        const size_t max_read,
                        size_t& num_read)
{
  prefetch_buffer.clear();

  auto length = info.sites();
  size_t number_left = std::min(number, max_read - num_read);

  while (number_left > 0 and iter)
  {
    const auto sequence_length = iter->length();
    if ( length and (length != sequence_length) ) {
      throw std::runtime_error{"MSA file does not contain equal size sequences"};
    }

    if (!length) length = sequence_length;

    prefetch_buffer.append( iter->label(), 
                            (premasking
                              ? subset_sequence(iter->sites(), info.gap_mask())
                              : iter->sites())
                          );

    length = sequence_length;
    number_left--;
    ++iter;
  }

  num_read += prefetch_buffer.size();
}

MSA_Stream::MSA_Stream( const std::string& msa_file,
                        const MSA_Info& info,
                        const bool premasking)
  : info_(info)
  , premasking_(premasking)
{
  iter_.from_file(msa_file);

  // ensure sequences are uniformly upper case
  iter_.reader().site_casing(genesis::sequence::FastaReader::SiteCasing::kToUpper);

  if (!iter_) {
    throw std::runtime_error{std::string("Cannot open file: ") + msa_file};
  }
}

size_t MSA_Stream::read_next( MSA_Stream::container_type& result, 
                              const size_t number)
{
  if (first_) {//...this is the first chunk
    // then read the first chunk and kick off the next asynchronously
    read_chunk(iter_, info_, premasking_, number, prefetch_chunk_, max_read_, num_read_);
    first_ = false;
  }
#ifdef __PREFETCH
  // join prefetching thread to ensure new chunk exists
  if (prefetcher_.valid()) {
    prefetcher_.wait();
  }
#endif
  // perform pointer swap to data
  std::swap(result, prefetch_chunk_);
  
  // start request next chunk from prefetcher (async)
#ifdef __PREFETCH
  prefetcher_ = std::async( std::launch::async,
                            read_chunk,
                            std::ref(iter_),
                            std::ref(info_),
                            premasking_,
                            number, 
                            std::ref(prefetch_chunk_),
                            max_read_,
                            std::ref(num_read_));
#else
  read_chunk(iter_, info_, premasking_, number, prefetch_chunk_, max_read_, num_read_);
#endif
  // return size of current buffer
  return result.size();
}

MSA_Stream::~MSA_Stream() 
{
#ifdef __PREFETCH
  // avoid dangling threads
  if (prefetcher_.valid()) {
    prefetcher_.wait();
  }
#endif
}

void MSA_Stream::constrain(const size_t max_read)
{
  max_read_ = max_read;
}

void MSA_Stream::skip_to_sequence(const size_t n)
{
  // this function is too dirty, disallow usage after first read
  if (not first_) {
    throw std::runtime_error{"Skipping currently not allowed after first read!"};
  }

 #ifdef __PREFETCH
  // join prefetching thread to ensure new chunk exists
  if (prefetcher_.valid()) {
    prefetcher_.wait();
  }
#endif 

  if (n >= num_sequences()) {
    throw std::runtime_error{"Trying to skip out of bounds!"};
  }

  // if (num_read_ == initial_size_) {
  //   num_read_ = 0;
  // }

  if (n < num_read_) {
    throw std::runtime_error{"Trying to skip behind!"};
  }

  size_t offset = n - num_read_;

  // seek the fileptr
  skip_ahead(offset);

  // kick off reading a chunk
// #ifdef __PREFETCH
//   prefetcher_ = std::async( std::launch::async,
//                             read_chunk, 
//                             std::ref(iter_),
//                             std::ref(info_),
//                             premasking_,
//                             initial_size_, 
//                             std::ref(prefetch_chunk_),
//                             max_read_,
//                             std::ref(num_read_));
// #else
//   read_chunk(iter_, info_, premasking_, n, prefetch_chunk_, max_read_, num_read_);
// #endif
}

void MSA_Stream::skip_ahead(const size_t n)
{
  std::advance(iter_, n);
}

size_t MSA_Stream::num_sequences()
{
  return info_.sequences();
}
