#pragma once

#include <string>
#include <vector>

class Sequence
{
public:
  Sequence()  = default;
  ~Sequence() = default;
  Sequence(std::string header, std::string sequence) 
    : sequence_(sequence) 
  {
    header_.push_back(header);
  }
  Sequence(const Sequence& s) = default;
  Sequence(Sequence&& s)      = default;

  // operator overloads
  Sequence& operator = (const Sequence& s)  = default;
  Sequence& operator = (Sequence&& s)       = default;
  bool operator==(const Sequence& other) {return sequence_.compare(other.sequence()) == 0;}
  bool operator==(const Sequence& other) const {return sequence_.compare(other.sequence()) == 0;}  
  // TODO doesn't merge in the full list (very tailored to the collapse func)
  void merge(const Sequence& other) {header_.push_back(other.header());}

  // member access
  const std::string& header() const {return header_.front();}
  const std::vector<std::string>& header_list() const {return header_;}
  const std::string& sequence() const {return sequence_;}

private:
  std::vector<std::string> header_;
  std::string sequence_;

};
