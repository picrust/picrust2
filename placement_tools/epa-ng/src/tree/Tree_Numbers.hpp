#pragma once

class Tree_Numbers
{
public:
  Tree_Numbers()  = default;
  ~Tree_Numbers() = default;

  Tree_Numbers(unsigned int tn)
    : tip_nodes(tn)
    , inner_nodes(tn - 2)
    , nodes(inner_nodes + tn)
    , branches(nodes -1)
  { }
  
  unsigned int tip_nodes = 0;
  unsigned int inner_nodes = 0;
  unsigned int nodes = 0;
  unsigned int branches = 0;
};
