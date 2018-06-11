#pragma once

#include <string>
#include <future>
#include <memory>
#include <sstream>
#include <cassert>

#include "sample/Sample.hpp"
#include "util/logging.hpp"
#include "io/jplace_util.hpp"

#ifdef __MPI
#include "net/epa_mpi_util.hpp"
#endif

class jplace_writer
{
public:
  jplace_writer() = default;
  jplace_writer(const std::string& out_dir,
                const std::string& file_name,
                const std::string& tree_string,
                const std::string& invocation_string)
  : tree_string_(tree_string)
  , invocation_(invocation_string)
  {
    init_mpi_();
    init_file_(out_dir, file_name);
  }

  ~jplace_writer()
  {
    // ensure last write/gather was completed
    wait();

    // finalize and close
    #ifdef __MPI

    if (local_rank_ == 0) {
      const auto trailing = finalize_jplace_string(invocation_);
      MPI_File_seek(shared_file_, 0, MPI_SEEK_END);
      MPI_File_write(shared_file_, trailing.c_str(), trailing.size(),
                      MPI_CHAR, MPI_STATUS_IGNORE);
    }
    MPI_File_close(&shared_file_);

    #else

    if (file_) {
      *file_ << finalize_jplace_string(invocation_);
      file_->close();
    }

    #endif
  }

  // jplace_writer& operator=( jplace_writer&& other )
  // {
  //   this->invocation_ = std::move(other.invocation_);
  //   this->file_ = std::move(other.file_);
  //   other.file_ = nullptr;
  //   this->prev_gather_ = std::move(other.prev_gather_);
  //   return *this;
  // }

  void write( Sample<>& chunk )
  {
    #ifdef __PREFETCH
    // ensure the last write has finished
    if (prev_gather_.valid()) {
      prev_gather_.get();
    }
    prev_gather_ = std::async(std::launch::async,
      [chunk = chunk, this]() mutable {
        this->write_(chunk);
      });
    #else
    write_(chunk);
    #endif
  }

  void wait()
  {
    #ifdef __PREFETCH
    if (prev_gather_.valid()) {
      prev_gather_.get();
    }
    #endif
  }

protected:

  void write_( Sample<>& chunk )
  {
    #ifdef __MPI // ========== MPI ==============

    if (shared_file_) {
      // serialize the sample
      std::stringstream buffer;
      if (first_){
        // account for the leading string
        if (local_rank_ == 0) {
          buffer << init_jplace_string(tree_string_);
        }
        first_ = false;
      } else {
        buffer << ",\n";
      }
      buffer << sample_to_jplace_string(chunk);

      // how much this rank intends to write this turn
      const auto buffer_str = buffer.str();
      size_t num_bytes = buffer_str.size();
      buffer.clear();

      // make the displacements known to all
      std::vector<size_t> block_sizes( all_ranks_.size() );
      MPI_Allgather(&num_bytes, 1, MPI_SIZE_T, block_sizes.data(), 1, MPI_SIZE_T, MPI_COMM_WORLD);
      size_t displacement = bytes_written_; // displacement is from cur file view
      size_t total_written = 0;
      for (size_t i = 0; i < block_sizes.size(); ++i) {
        if ( i < static_cast<size_t>(local_rank_) ) {
          displacement += block_sizes[i];
        }
        total_written += block_sizes[i];
      }

      // write the local chunk
      MPI_File_write_at_all(shared_file_,
                            displacement,
                            buffer_str.c_str(),
                            buffer_str.size(),
                            MPI_CHAR,
                            MPI_STATUS_IGNORE);

      bytes_written_ += total_written;
    }

    #else // ========== NOT MPI ==============

    if (file_) {
      if (first_){
        first_ = false;
        *file_ << init_jplace_string(tree_string_);
      } else {
        *file_ << ",\n";
      }

      *file_ << sample_to_jplace_string(chunk);
    }

    #endif
  }

  virtual void init_file_(const std::string& out_dir,
                          const std::string& file_name)
  {
    const auto file_path = out_dir + file_name;
    #ifdef __MPI
    MPI_File_open(MPI_COMM_WORLD,
              file_path.c_str(),
              MPI_MODE_WRONLY | MPI_MODE_CREATE,
              MPI_INFO_NULL,
              &shared_file_);
    #else
    file_ = std::make_unique<std::fstream>();
    file_->open(file_path,
                std::fstream::in | std::fstream::out | std::fstream::trunc);

    if (not file_->is_open()) {
      throw std::runtime_error{file_path + ": could not open!"};
    }
    #endif
  }

  void init_mpi_()
  {
    #ifdef __MPI // then have one outfile per rank
    int num_ranks = 0;
    MPI_COMM_RANK(MPI_COMM_WORLD, &local_rank_);
    MPI_COMM_SIZE(MPI_COMM_WORLD, &num_ranks);

    all_ranks_.resize(num_ranks);
    for (int i = 0; i < num_ranks; ++i) {
      all_ranks_[i] = i;
    }

    // ensure non-rank-0 blocks start with a comma
    if (local_rank_ != 0) {
      first_ = false;
    }
    #endif
  }

protected:
  std::string tree_string_;
  std::string invocation_;
  std::future<void> prev_gather_;
  bool first_ = true;
  
  #ifdef __MPI
  MPI_File shared_file_;
  size_t bytes_written_ = 0;
  int local_rank_ = 0;
  std::vector<int> all_ranks_ = {0};
  #else
  std::unique_ptr<std::fstream> file_ = nullptr;
  #endif
};