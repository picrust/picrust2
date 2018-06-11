#pragma once

#ifdef __MPI

#include "net/mpihead.hpp"
#include "util/Timer.hpp"
#include "util/logging.hpp"

#include <sstream>
#include <memory>
#include <unordered_map>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>

#include <stdint.h>
#include <limits.h>

#if SIZE_MAX == UCHAR_MAX
   #define MPI_SIZE_T MPI_UNSIGNED_CHAR
#elif SIZE_MAX == USHRT_MAX
   #define MPI_SIZE_T MPI_UNSIGNED_SHORT
#elif SIZE_MAX == UINT_MAX
   #define MPI_SIZE_T MPI_UNSIGNED
#elif SIZE_MAX == ULONG_MAX
   #define MPI_SIZE_T MPI_UNSIGNED_LONG
#elif SIZE_MAX == ULLONG_MAX
   #define MPI_SIZE_T MPI_UNSIGNED_LONG_LONG
#else
   #error "what is happening here?"
#endif

// inline auto epa_mpi_file( const int local_rank,
//                           const size_t section_size)
// {


//   return ;
// }

// types to keep track of previous async sends
typedef struct
{
  MPI_Request req;
  char*       buf = nullptr;
} request_tuple;

using previous_request_storage_t = typename std::unordered_map<int, request_tuple>;

static void err_check(int errval)
{
  if (errval != MPI_SUCCESS) {
    std::stringstream msg;
    msg << "Failed MPI call: ";
    switch (errval) {
      case MPI_ERR_COMM:
        msg << "Invalid communicator.";
        break;
      case MPI_ERR_TYPE:
        msg << "Invalid datatype argument.";
        break;
      case MPI_ERR_COUNT:
        msg << "Invalid count argument.";
        break;
      case MPI_ERR_TAG:
        msg << "Invalid tag argument.";
        break;
      case MPI_ERR_RANK:
        msg << "Invalid source or destination rank.";
        break;
      default:
        msg << "unknown";
    }
    msg << std::endl;
    throw std::runtime_error{msg.str()};
  }
}

void epa_mpi_waitall(previous_request_storage_t& reqs)
{
  for (auto& pair : reqs) {
    auto& r = pair.second;
    if (r.req) {
      MPI_Status status;
      err_check(MPI_Wait(&r.req, &status));
      delete[] r.buf;
      r.buf = nullptr;
    }
  }
}

template <typename T>
void epa_mpi_send(T& obj,
                  const int dest_rank,
                  const MPI_Comm comm)
{
  // serialize the obj
  std::stringstream ss;
  cereal::BinaryOutputArchive out_archive(ss);
  out_archive(obj);

  // send obj to specified node
  std::string data = ss.str();
  auto buffer = new char[data.size()];
  memcpy(buffer, data.c_str(), data.size() * sizeof(char));
  err_check( MPI_Send(buffer,
                      data.size(),
                      MPI_CHAR,
                      dest_rank,
                      0,
                      comm));
  delete[] buffer;
}

template <typename T>
void epa_mpi_isend( T& obj,
                    const int dest_rank,
                    const MPI_Comm comm,
                    request_tuple& prev_req,
                    Timer<>& timer)
{
  // wait for completion of previous send
  if (prev_req.req) {

    MPI_Status status;
    timer.pause();
    LOG_DBG2 << "previous request detected, calling wait...";
    err_check(MPI_Wait(&prev_req.req, &status));
    LOG_DBG2 << "Done!";
    timer.resume();
    delete[] prev_req.buf;
    // free previous request?
  }

  // serialize the obj
  std::stringstream ss;
  cereal::BinaryOutputArchive out_archive(ss);
  out_archive(obj);

  // send obj to specified node
  std::string data = ss.str();
  auto buffer = new char[data.size()];
  memcpy(buffer, data.c_str(), data.size() * sizeof(char));
  err_check( MPI_Issend(buffer,
                        data.size(),
                        MPI_CHAR,
                        dest_rank,
                        0,
                        comm,
                        &prev_req.req));
  
  prev_req.buf = buffer;
}

template <typename T>
void epa_mpi_receive( T& obj,
                      const int src_rank,
                      const MPI_Comm comm,
                      Timer<>& timer)
{
  // probe to find out the message size
  MPI_Status status;
  int size = 0;
  timer.pause();
  err_check( MPI_Probe( src_rank,
                        MPI_ANY_TAG,
                        comm,
                        &status) );
  timer.resume();

  LOG_DBG1 << "Receiving data from rank "
          << status.MPI_SOURCE;

  MPI_Get_count(&status, MPI_CHAR, &size);

  LOG_DBG1 << "of size: "
          << size << " bytes";

  // prepare buffer
  auto buffer = new char[size];

  //  get the actual payload
  err_check( MPI_Recv(buffer,
                      size,
                      MPI_CHAR,
                      status.MPI_SOURCE,
                      status.MPI_TAG,
                      comm,
                      &status) );

  LOG_DBG1 << "Done!";

  // deserialization
  LOG_DBG1 << "Deserializing...";
  std::stringstream ss;
  ss.write(buffer, size);
  cereal::BinaryInputArchive in_archive(ss);
  // build the object
  in_archive(obj);
  LOG_DBG1 << "Done!";

  delete[] buffer;
}

template <typename T>
static inline void isend_all( const std::vector<T>& parts,
                              const std::vector<int>& dest_ranks,
                              const MPI_Comm comm,
                              previous_request_storage_t& prev_reqs,
                              Timer<>& timer)
{
  for (size_t i = 0; i < parts.size(); ++i) {
    auto dest = dest_ranks[i];
    epa_mpi_isend(parts[i],
                  dest,
                  comm,
                  prev_reqs[dest],
                  timer);
  }
}

template <typename T>
void epa_mpi_split_send(T& obj,
                        const std::vector<int>& dest_ranks,
                        const MPI_Comm comm,
                        previous_request_storage_t& prev_reqs,
                        Timer<>& timer)
{
  LOG_DBG1 << "Sending...";

  std::vector<T> parts;
  split(obj, parts, dest_ranks.size());

  // TODO only relevant if T conforms with Token
  for (auto& p : parts) {
    p.status(obj.status());
  }

  isend_all(parts, dest_ranks, comm, prev_reqs, timer);

  LOG_DBG1 << "Done!";
}

template <typename T>
void epa_mpi_receive_merge( T& obj,
                            const std::vector<int>& src_ranks,
                            const MPI_Comm comm,
                            Timer<>& timer)
{
  for (const auto rank : src_ranks) {
    T remote_obj;
    epa_mpi_receive(remote_obj, rank, comm, timer);
    merge(obj, remote_obj);

    // TODO only relevant if T conforms with Token
    obj.status(remote_obj.status());
  }
}

template <typename T>
void epa_mpi_gather(T& obj,
                    const int dest_rank,
                    const std::vector<int>& src_ranks,
                    const int local_rank,
                    Timer<>& timer)
{
  if (local_rank == dest_rank) {
    for (const auto src_rank : src_ranks) {
      if (local_rank == src_rank) {
        continue;
      }
      T remote_obj;
      LOG_DBG1 << "Gather loop rank: " << src_rank;
      epa_mpi_receive(remote_obj,
                      MPI_ANY_SOURCE,
                      MPI_COMM_WORLD,
                      timer);
      merge(obj, std::move(remote_obj));
    }
    // check for equal entries and collapse them into one
    collapse(obj);
  } else {
    LOG_DBG1 << "Gather send to rank: " << dest_rank;
    epa_mpi_send(obj, dest_rank, MPI_COMM_WORLD);
  }
}

template <typename T>
void epa_mpi_bcast( T& obj,
                    const int src_rank,
                    const std::vector<int>& dest_ranks,
                    const int local_rank,
                    Timer<>& timer)
{
  if (src_rank == local_rank) {
    for (auto dest_rank : dest_ranks) {
      if (local_rank == dest_rank) {
        continue;
      }
      epa_mpi_send(obj, dest_rank, MPI_COMM_WORLD);
    }
  } else {
    epa_mpi_receive(obj, src_rank, MPI_COMM_WORLD, timer);
  }
}

#endif
