#pragma once

#include <numeric>
#include <map>
#include <cereal/types/map.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/base_class.hpp>

#include "sample/Sample.hpp"
#include "pipeline/Token.hpp"

// forward declaration
class WorkIterator;

/**
 * Container to hold sequence id's and their corresponding branches:
 * work[branch_id] = {seq_id_1. seq_id_2, ...}
 *
 * Meant as a structure that can be used by nodes to figure out what to compute.
 */
class Work : public Token
{
public:
  using key_type              = size_t;
  using value_type            = size_t;
  using container_value_type  = std::vector<value_type>;
  using container_type        = std::map<key_type, container_value_type>;
  using const_iterator        = WorkIterator;
  // typedef typename container_type::iterator         iterator;
  
  struct Work_Pair
  {
      key_type    branch_id;
      value_type  sequence_id;
  };

  /**
   * Create work object from a Sample: all entries are seen as placements to be recomputed
   */
  template<class T>
  Work(Sample<T>& sample)
  {
    for (auto& pq : sample)
    {
      const auto seq_id = pq.sequence_id();
      for (auto& placement : pq)
      {
        const auto branch_id = placement.branch_id();
        work_set_[branch_id].push_back(seq_id);
      }
    }
  }

  /**
   * Create a work object covering all sequences in [seq_range.first, seq_range.second)
   * for every branch ID in [branch_range.first, branch_range.second).
   */
  Work(std::pair<key_type, key_type>&& branch_range, std::pair<value_type, value_type>&& seq_range)
  {
    for (key_type branch_id = branch_range.first; branch_id < branch_range.second; ++branch_id) {
      for (value_type seq_id = seq_range.first; seq_id < seq_range.second; ++seq_id) {
        work_set_[branch_id].push_back(seq_id);      
      }
    }
  }

  Work(Work const& other) = default;
  Work(Work && other) = default;

  Work& operator=(Work const&) = default;
  Work& operator=(Work &&) = default;

  Work() = default;       
  ~Work() = default;

  // methods
  void clear() {work_set_.clear();}

  size_t size() const 
  {
    return std::accumulate(work_set_.begin(), work_set_.end(), 0,
      [](size_t a, std::pair<key_type, container_value_type> b){return a + b.second.size();}
      );
  }

  inline void add(key_type branch_id, value_type seq_id)
  {
    work_set_[branch_id].emplace_back(seq_id);
  }

  inline void add(Work_Pair& it);

  // template< class InputIt >
  // void insert(InputIt first, InputIt last) {work_set_.insert(first, last);}

  // Iterator Compatibility
  auto bin_begin() { return work_set_.begin(); }
  auto bin_end() { return work_set_.end(); }
  auto bin_cbegin() const { return work_set_.begin(); }
  auto bin_cend() const { return work_set_.end(); }
  const_iterator begin() const;
  const_iterator end() const;
  // const_iterator cbegin() { return work_set_.cbegin(); }
  // const_iterator cend() { return work_set_.cend(); }

  // Operator overloads
  const container_value_type& at (const key_type index) const { return work_set_.at(index); }
  container_value_type& operator[] (const key_type index) { return work_set_[index]; }
  
  // serialization
  template <class Archive>
  void serialize(Archive & ar) 
  { ar( *static_cast<Token*>( this ), work_set_ ); }


private:
  container_type work_set_;
};

class WorkIterator
{
public:
    // -----------------------------------------------------
    //     Typedefs
    // -----------------------------------------------------
    using WorkType      = Work::container_type;
    using MapItType     = typename WorkType::const_iterator;
    using VecItType     = typename WorkType::mapped_type::const_iterator;
    using self_type     = WorkIterator;
    using element_type  = Work::Work_Pair;
    using iterator_tag  = std::forward_iterator_tag;

    // -----------------------------------------------------
    //     Constructors and Rule of Five
    // -----------------------------------------------------

    WorkIterator() = delete;

    WorkIterator( WorkType const& target, bool is_end )
        : map_it( std::begin( target ) )
        , end_it( std::end( target ) )
        , vec_it()
    {
        if( is_end ) {
            map_it = end_it;
        }
        if( map_it != end_it && ! target.empty() ) {
            while( map_it != end_it && map_it->second.empty() ) {
                ++map_it;
            }
            if( map_it != end_it && ! map_it->second.empty() ) {
                vec_it = std::begin( map_it->second );
            }
        }
    }

    ~WorkIterator() = default;

    WorkIterator( WorkIterator const& ) = default;
    WorkIterator( WorkIterator&& ) = default;
    WorkIterator& operator= ( WorkIterator const& ) = default;
    WorkIterator& operator= ( WorkIterator&& ) = default;

    // -----------------------------------------------------
    //     Operators
    // -----------------------------------------------------

    element_type operator * ()
    {
        return { map_it->first, *vec_it };
    }

    // self_type const* operator -> ()
    // {
    //     return this;
    // }

    size_t current_branch_id()
    {

      return map_it->first;
    }

    size_t current_sequence_id()
    {
        return *vec_it;
    }

    self_type operator ++ ()
    {
        ++vec_it;
        if( vec_it == map_it->second.end() ) {
            do {
                ++map_it;
                if( map_it == end_it ) {
                    break;
                } else {
                    vec_it = std::begin( map_it->second );
                }
            } while( vec_it == std::end( map_it->second ) );
        }
        return *this;
    }

    self_type operator ++ (int)
    {
        self_type tmp = *this;
        ++(*this);
        return tmp;
    }

    bool operator == (const self_type &other) const
    {
        return ( other.map_it == map_it && other.end_it == end_it ) &&
               ( map_it == end_it || other.vec_it == vec_it );
    }

    bool operator != (const self_type &other) const
    {
        return !(other == *this);
    }

    // -----------------------------------------------------
    //     Data Members
    // -----------------------------------------------------

private:

    MapItType map_it;
    MapItType end_it;

    VecItType vec_it;
};

inline void Work::add(Work_Pair& it)
{
  add(it.branch_id, it.sequence_id);
}

inline Work::const_iterator Work::begin() const
{
    return WorkIterator( work_set_, false );
}

inline Work::const_iterator Work::end() const
{
    return WorkIterator( work_set_, true );
}