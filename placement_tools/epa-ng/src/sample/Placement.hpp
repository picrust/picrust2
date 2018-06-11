#pragma once

#include <cstddef>

class Slim_Placement;

class Placement {
public:
  Placement() = default;

  Placement(const size_t branch_id,
            const double likelihood, 
            const double pendant_length, 
            const double distal_length)
    : branch_id_(branch_id)
    , likelihood_(likelihood)
    , lwr_(0.0)
    , pendant_length_(pendant_length)
    , distal_length_(distal_length)
  { }

  Placement(const Slim_Placement& other);

  Placement(Placement const& other) = default;
  Placement(Placement&& other)      = default;

  Placement& operator= (Placement const& other) = default;
  Placement& operator= (Placement && other)     = default;

  ~Placement() = default;

  // getters
  double lwr() const {return lwr_;};
  double likelihood() const {return likelihood_;}
  double pendant_length() const {return pendant_length_;};
  double distal_length() const {return distal_length_;};
  size_t branch_id() const {return branch_id_;}

  // setters
  void lwr(double value) {lwr_ = value;};
  void likelihood(double value) {likelihood_ = value;};
  void pendant_length(double value) {pendant_length_ = value;};
  void distal_length(double value) {distal_length_ = value;};

  // serialization
  template<class Archive>
  void serialize(Archive& ar) { ar(branch_id_, likelihood_, lwr_, pendant_length_, distal_length_); }

private:
  size_t branch_id_;
  double likelihood_;
  double lwr_;
  double pendant_length_;
  double distal_length_;
};

class Slim_Placement
{
public:
  Slim_Placement() = default;

  Slim_Placement( const size_t branch_id,
                  const double likelihood,
                  const double, 
                  const double)
    : branch_id_(branch_id)
    , likelihood_(likelihood)
  { }
  Slim_Placement(const Placement& other);

  ~Slim_Placement() = default;

  double likelihood() const {return likelihood_;}
  size_t branch_id() const {return branch_id_;}

  // serialization
  template<class Archive>
  void serialize(Archive& ar) { ar(branch_id_, likelihood_); }

private:
  size_t branch_id_;
  double likelihood_;
  
};
