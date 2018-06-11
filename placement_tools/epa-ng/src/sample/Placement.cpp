#include "sample/Placement.hpp"

Placement::Placement(const Slim_Placement& other)
    : branch_id_(other.branch_id())
    , likelihood_(other.likelihood())
    , lwr_(0.0)
    , pendant_length_(0.0)
    , distal_length_(0.0)
  { }

Slim_Placement::Slim_Placement(const Placement& other)
    : branch_id_(other.branch_id())
    , likelihood_(other.likelihood()) 
  { }
