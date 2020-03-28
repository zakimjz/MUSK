#ifndef _LATTICE_NODE_H
#define _LATTICE_NODE_H

#include "graph_iso_check.h"
#include <algorithm>
#include "helper_funs.h"
#include "pattern_factory.h"

template <typename PAT >
struct lattice_node
{

  typedef lattice_node<PAT> L_NODE;
  typedef typename PAT::VERTEX_T V_T;
  typedef typename PAT::EDGE_T E_T;
  typedef pair<int, int> EDGE;

  string get_key() {
    const typename PAT::CAN_CODE& cc = _pat->canonical_code();
    std::string min_dfs_cc = cc.to_string();
    return min_dfs_cc;
  }
 
  lattice_node(PAT* p) {
    _pat = p;
    _is_frequent = false;
    _is_maximal = false;
    _super_cnt =0; // number of super-neighbors
  }
  
  bool _is_maximal;
  bool _is_frequent;
  bool _is_processed; // it is true, when we know all neighbors and their status of this pattern
  PAT* _pat;
  vector<L_NODE*> _neighbors;
  vector<int> _neighbor_status;  // -1 for infrequent, 1 for frequent, #degree for maximal
  int  _super_cnt; // number of super-neighbors
};

#endif
