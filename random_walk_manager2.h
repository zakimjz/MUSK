#ifndef _RANDOM_WALK_MANAGER2_H_
#define _RANDOM_WALK_MANAGER2_H_

#include <unordered_map>
#include "graph_iso_check.h"
#include <algorithm>
#include "helper_funs.h"
#include "pattern_factory.h"
#include "lattice_node.h"
#include "random.h"

/**
 * Structure storing all the lattice nodes
 */
template<class PAT>
class RandomWalkManager2
{

  public:

  typedef lattice_node<PAT> LATTICE_NODE;
  typedef unordered_map<string, LATTICE_NODE*> NODE_MAP;
  typedef typename NODE_MAP::iterator NS_IT;
  typedef typename NODE_MAP::const_iterator CONST_NS_IT;
  typedef Database<PAT> DATABASE;
  typedef PatternFactory<PAT> PATTERN_FACTORY;

  RandomWalkManager2(double c, int t, DATABASE* d) {
    _C = c; 
    _last = 0;
    _pf = PATTERN_FACTORY::instance(d);
    _max_cnt = t;
    _cnt = 0;
    _fcnt = 0;
  }

  PatternFactory<PAT>* get_pat_factory() const {
    return _pf;
  }

  // random walk manager initialize with a maximal pattern node
  LATTICE_NODE* initialize() {
    PAT* p =_pf->get_one_maximal_pattern();
    // exit(1);
    const typename PAT::CAN_CODE& cc = check_isomorphism(p);
    p->set_canonical_code(cc);
    // cout << p->get_canonical_code().to_string();
    // p->print_vat();
    _pf->get_database()->verify_vat(p);
    //exit(1);

#ifdef PRINT
    cout << "In initialize: got one maximal pattern\n";
    if (p->_is_frequent) {
      cout << "Yes, it is frequent\n";
    }
    if (p->_is_maximal) {
      cout << "First maximal pattern found\n";
    }
    else {
      cout << "ERROR:this pattern is infrquent\n";
      exit(1);
    }
#endif
    LATTICE_NODE* ln = create_lattice_node(p);
    process_node(ln);
    return ln;
  }

  void walk(LATTICE_NODE* current) {
    string c_code;
    while (true) {
      process_node(current);
      cout << "\nF ("<<++_fcnt<<"):";
      // exit(1);
      PAT* p = current->_pat;
      if (p->is_max()) {
       	c_code = p->get_canonical_code().to_string();
        pair<set<string>::iterator, bool> r = _max_patterns.insert(c_code);
        if (r.second == true) {
        	cout << "\nM ("<<++_cnt<<"): ";
        	cout << c_code;
        	if (_cnt == _max_cnt) {return;}
        }
      }
      // cout << "In walk: looking for new node to visit" << endl;
      LATTICE_NODE* next = get_next(current);
      _last = p;
      current = next;
      _uniform = false;
    }
  }

  LATTICE_NODE* get_next(LATTICE_NODE* current) const {
    int total = current->_neighbors.size();
    assert(current->_neighbor_status.size() == total);
    vector<double> weight(total,0.0);
    vector<double> prob(total,0.0);
    if (_uniform == true) { // currently visiting pattern is a maximal pattern
      int idx = boost_get_a_random_number(0, total);
      // cout << "choosing " << idx << " from " <<  total << " neighbors" << endl;
      return current->_neighbors[idx];  
    }
    else {
      int max_maximal_degree = 1;
      for (int i=0; i<total; i++) {
        int val = current->_neighbor_status[i]; 
        if (val > max_maximal_degree) max_maximal_degree = val;
      }

      if (max_maximal_degree == 1) {
        for (int i=0; i<total; i++) {
          weight[i] = (i<current->_super_cnt)? 2.0 : 1.0; // we go up with double probability than going down
        }
      }
      else {
        for (int i=0; i<total; i++) {
          int val = current->_neighbor_status[i]; 
          if (val > 1)  // maximal pattern
            weight[i] = _C/val;  // The sum of of incoming weights for every maximal patterns is equal to _C
          else 
            weight[i] = (i<current->_super_cnt)? _C/(3*max_maximal_degree) : _C/(6*max_maximal_degree);
          // cout << weight[i] << " ";
        }
        // cout << endl;
      }
      std::partial_sum(weight.begin(), weight.end(), prob.begin());
      double divisor = prob[prob.size()-1];
      for (int i=0; i<prob.size(); i++) {
        prob[i] = prob[i] / divisor; 
      }
      int idx = randomWithDiscreteProbability(prob);
      // cout << "choosing " << idx << " from " << prob.size() << " neighbors" << endl;
      return current->_neighbors[idx];  
    }
  }

  LATTICE_NODE* create_lattice_node(PAT*& p) {
    const typename PAT::CAN_CODE& cc = check_isomorphism(p);
    p->set_canonical_code(cc);
    std::string min_dfs_cc = cc.to_string();
    LATTICE_NODE* node = exists(min_dfs_cc);
    if (node == 0) {  // new pattern
      node = new LATTICE_NODE(p);
      node->_is_processed = false;
      insert_lattice_node(min_dfs_cc, node);
    }
    else {
      delete p;
      p = node->_pat;
    }
    return node;
  }

  LATTICE_NODE* exists(string p) {;
    CONST_NS_IT it = _node_map.find(p);
    return (it != _node_map.end())? it->second : 0;
  }

  void insert_lattice_node(string p, LATTICE_NODE* ln) {
    _node_map.insert(make_pair(p, ln));
  }

  void process_node(LATTICE_NODE* n) {
    // cout << "In process node, uniform=" << _uniform << endl;
    if (n->_is_processed) return;
    PAT* p = n->_pat;
    assert(p->get_sup_ok());
    // cout << "Current pattern:\n";
    // cout << *p;
    vector<PAT*> neighbors;
    if (p->is_max()) {
      _pf->get_sub_patterns(p, neighbors); 
      _uniform = true;
    }
    else if (p->is_frequent()) {
      _pf->get_freq_super_patterns(p, neighbors); 
      n->_super_cnt = neighbors.size();
      _pf->get_sub_patterns(p, neighbors); 
    }
    else {
      cout << "ERROR: In this walk we are traversing only frequent patterns (1)\n";
      exit(1);
    }
#ifdef PRINT
    cout << "Its neighbors:" << endl;
#endif
    for (int i=0; i<neighbors.size(); i++) {
      PAT* one_neighbor = neighbors[i];
#ifdef PRINT
      cout << *one_neighbor;
#endif
      LATTICE_NODE* ln = create_lattice_node(one_neighbor);
      int status;
      n->_neighbors.push_back(ln);
      if (_uniform == true) {
        ln->_is_frequent = true;
        ln->_is_maximal = false;
        one_neighbor->set_status_known();
      }
      else {
        _pf->get_status_optimal(one_neighbor, ln->_is_frequent, ln->_is_maximal);
        
        const typename PAT::CAN_CODE& cc = check_isomorphism(one_neighbor);
        one_neighbor->set_canonical_code(cc);
        // cout << one_neighbor->get_canonical_code().to_string();
        // one_neighbor->print_vat();
      }
      if (ln->_is_maximal) {
#ifdef PRINT
        cout << "Neighbor " << i << " is maximal" << endl;
#endif
        int its_degree = one_neighbor->find_removable_edges();
        n->_neighbor_status.push_back(its_degree);
        _uniform = false;
      }
      else if (ln->_is_frequent) {
#ifdef PRINT
        cout << "This neighbor is frequent" << endl;
#endif
        n->_neighbor_status.push_back(1);
      }
      else {// infrequent
        cout << "ERROR: This neighbor " << i << " is infrequent, we are traversing only frequent patterns (2)\n";
        exit(1);
      }
    }  
    n->_is_processed=true;
  }

  private:
  NODE_MAP _node_map;
  PatternFactory<PAT>* _pf;
  double _C;
  PAT* _last;
  bool _uniform;
  unsigned int _max_cnt;
  unsigned int _cnt;
  unsigned int _fcnt;
  set<string> _max_patterns;
};

#endif
