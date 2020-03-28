#ifndef _PATTERN_H_
#define _PATTERN_H_

#include "adj_matrix.h"
#include "graph_can_code.h"
#include <algorithm>

// For graph pattern that does not have any edge label
template <typename V_T>
class Pattern 
{
  public:
    typedef AdjMatrix<V_T> ST;
    typedef std::pair<V_T, V_T> EDGE;
    Pattern(const vector<V_T>& labels) { 
      _matrix = new ST(labels);
    }

  private:
    vector<int> _vat;
    multiset<EDGE> _edges;
    ST* _matrix;
}; //forward declaration for friends
template <typename V_T, typename E_T>
class ExPattern;

template <typename V_T, typename E_T>
ostream& operator<< (ostream& ostr, const ExPattern<V_T, E_T>& );


// for graph pattern that has both vertex and edge label
template <typename V_T, typename E_T>
class ExPattern
{

  friend ostream& operator<< <>(ostream& osr, const ExPattern<V_T, E_T>&);

  public:
    typedef ExPattern<V_T, E_T> PAT;
    typedef V_T VERTEX_T;
    typedef E_T EDGE_T;
    typedef FullLabelAdjMatrix<V_T, E_T> ST;
    typedef std::pair<pair<V_T, V_T>, E_T> EDGE;
    typedef typename multiset<EDGE>::iterator EDGE_IT;
    typedef typename multiset<EDGE>::const_iterator EDGE_CIT;
    typedef canonical_code<V_T, E_T> CAN_CODE;
    typedef typename CAN_CODE::INIT_TYPE CC_INIT_TYPE;
    typedef typename CAN_CODE::COMPARISON_FUNC CC_COMPARISON_FUNC;
    typedef set<pair<V_T, E_T> > NEIGHBORS;
    typedef typename NEIGHBORS::const_iterator NCIT;
    typedef typename NEIGHBORS::iterator NIT;
    
    ExPattern() {
      _matrix = 0;
      _sup_ok = false;    
      _is_frequent = false;
      _is_maximal=false;
      _status_known=false;
    }

    ExPattern(const vector<V_T>& labels) {
      _matrix = new ST(labels);
      _sup_ok = false;
      _is_frequent = false;
      _is_maximal=false;
      _status_known=false;
    }

    ExPattern(ST* p) {
      ST* st = new ST(*_matrix);
      _matrix = st;
      _sup_ok = false;    
      _is_frequent = false;
      _is_maximal=false;
      _status_known=false;
   }

   ~ExPattern() {
     delete _matrix;
   }

   void add_edge(unsigned int i, unsigned int j, E_T e) {
      _matrix->add_edge(i,j,e);
      V_T v1 = _matrix->label(i);
      V_T v2 = _matrix->label(j);
      pair<V_T, V_T> v_pair = (v1 < v2) ? make_pair(v1, v2) : make_pair(v2,v1);
      EDGE edge = make_pair(v_pair, e);
      _edges.insert(edge);
      _sup_ok = false;  // once you add edges, you need to update vat
    }
    PAT* make_null_pattern(int max_id_in_vat) const {
      PAT* clone = new PAT();
      clone->_vat.resize(max_id_in_vat); 
      for (int i=0; i<=clone->_vat.size();i++)
        clone->_vat[i]=i;
      clone->_sup_ok=true;
      clone->_is_maximal=false;
      clone->_is_frequent=true;
      clone->_status_known=true;
      return clone;
    }
    PAT* clone() const {
#ifdef PRINT
      cout << "In clone:\n";
#endif
      PAT* clone = new PAT();
      clone->_matrix = this->_matrix->clone();
      clone->_vat = _vat; 
      clone->_edges = _edges;
      clone->_is_frequent= false;
      clone->_is_maximal = false;
      clone->_removable_edge_known =false;
#ifdef PRINT
      cout << "Got the new pattern, returing from clone ..." << endl;
#endif
      return clone;
    }
    
    E_T get_edge_label(unsigned int i, unsigned int j) const {
      return _matrix->get_edge_label(i, j);
    }
 
    V_T label(unsigned int i) const {
      return _matrix->label(i);
    }

    size_t size() const {
      if (_matrix)
        return _matrix->size();
      else
        return 0;
    }

    void get_vids_for_this_label(const V_T& v, vector<int>& ret_val) const {
      _matrix->get_vid_for_this_label(v, ret_val);
    }

    EDGE remove_edge(const int&a, const int& b) {  // returning the removed edge
      V_T v1 = this->label(a);
      V_T v2 = this->label(b);
      E_T e = _matrix->get_edge_label(a,b);

      _matrix->remove_edge(a,b);
      pair<V_T, V_T> v_pair = (v1 < v2) ? make_pair(v1, v2) : make_pair(v2,v1);
      EDGE edge = make_pair(v_pair, e);
      EDGE_IT eit = _edges.find(edge);
      if (eit == _edges.end()) {
        cout << "ERROR: request for non-existent edge removal\n";
        exit(1);
      }
      _edges.erase(eit);
      _sup_ok = false;
      return edge;
    }
    
    unsigned int find_removable_edges() {
      if (_removable_edge_known) return _removable_edges.size();
      AdjIterator iter(_matrix);

      for (iter.first(); !iter.is_done(); iter.next()) {
        pair<int, int> edge = iter.current();
        if (! _matrix->essential_edge(edge.first, edge.second)) {
          // cout << "edge (" << edge.first << "," << edge.second << ") is removable\n";
          _removable_edges.push_back(edge); 
        }
      }
      _removable_edge_known = true;
      return _removable_edges.size();
    }

    const vector<pair<int, int> >& get_removable_edges() const {
      return _removable_edges;
    } 

    const ST* get_adj_matrix() const {
      return _matrix;
    }

    void set_vat(const vector<int>& v) {
      _vat = v;
    }

    void set_sup_ok() {
      _sup_ok=true;
    }

    bool get_sup_ok() const {
      return _sup_ok;
    }

    bool edge_exist(int i, int j) const {
      return _matrix->at(i,j);
    }

    int add_vertex(V_T v) {
      return _matrix->add_vertex(v);
    }

    bool is_frequent() const {
      assert(_sup_ok == true);
      return _is_frequent;
    }
    
    void set_canonical_code(const CAN_CODE& c) {
      _canonical_code = c; 
    }

    const multiset<EDGE>& get_edgeset() const {
      return _edges;
    }
    bool is_super_pattern(PAT* pat) {
      const multiset<EDGE>& mset = pat->get_edgeset();
      if (includes(_edges.begin(), _edges.end(), mset.begin(), mset.end())) {
        return true;
      }
      else {
        return false;
      }
    }

    const vector<int>& get_vat() const {
      return _vat;
    }

    void print_vat() const {
      cout <<"\nVAT:Size(" << _vat.size() << ")\n";
      std::copy(_vat.begin(), _vat.end(), ostream_iterator<int>(cout," "));
      cout << endl;
    }

    void join_vat(PAT* p) {
#ifdef PRINT
      cout << "Before join: vat size:" << _vat.size() << endl;
#endif
      vector<int> out_vector;
      set_intersection(_vat.begin(), _vat.end(), p->_vat.begin(), p->_vat.end(),
                back_inserter(out_vector));
      _vat = out_vector;
#ifdef PRINT
      cout << "After join: vat size:" << _vat.size() << endl;
#endif
    }

    void set_freq() {
      _is_frequent = true;
    }

    const CAN_CODE& get_canonical_code() const { return _canonical_code;} 
 
    void add_tid_to_vat(int tid) {
      vector<int>::iterator it;
      it = lower_bound(_vat.begin(), _vat.end(), tid);
      _vat.insert(it, tid);
    }

    /* void add_tid_to_vat(vector<int>::iterator it, int tid) const { */
    /*   _vat.insert(it, tid); */
    /* } */

    int edge_counter(const EDGE& e) const {
      pair<EDGE_IT, EDGE_IT> eit_p = equal_range(_edges.begin(), _edges.end(), e);
      int cnt=0;
      for (;eit_p.first != eit_p.second;eit_p.first++) cnt++;
      return cnt;
    }

    bool status_known() const {
      return _status_known;
    }

    void set_status_known() {
      _status_known = true;
    }

    bool is_max() const {
      assert(_status_known == true);
      return _is_maximal;
    }

    void set_max() {
      _status_known = true;
      _is_maximal = true;
      _is_frequent = true;
    }
  
    int sub_pat_cnt() const {
      return _removable_edges.size();
    }
#ifdef MARGIN
    void set_last_edge(int a, int b) {
      _last_edge.first = a;
      _last_edge.second = b;
      _last_edge_exist=true;
    }
 
    pair<int, int> get_last_edge() const {
      assert (_last_edge_exist == true);
      return _last_edge;
    }
    
#endif

    int edge_cnt() const {
      return _edges.size();
    }

    
    bool _is_frequent;
    bool _is_maximal;
  private:

    vector<int> _vat;           // contains the graph ids that has at least these edges
    multiset<EDGE> _edges;
    vector<pair<int, int> > _removable_edges;  // removing these edges keeps the pattern connected
    ST* _matrix;
    bool _sup_ok;               // if the exact support of this pattern computed?
    CAN_CODE _canonical_code;
    bool _status_known;
    bool _removable_edge_known;

#ifdef MARGIN
    pair <int, int> _last_edge;
    bool _last_edge_exist;
#endif
};

template<typename V_TYPE, typename E_TYPE>
ostream& operator<< (ostream& ostr, const ExPattern<V_TYPE, E_TYPE>& p){

  typedef typename ExPattern<V_TYPE, E_TYPE>::EDGE EDGE;
  if (p.size() == 0) {
    cout << "NULL" << endl;
    return ostr;
  }
  cout << "SIZE:" <<  p.size() << endl;
  cout << "EDGE COUNT:" << p._edges.size() << endl;
  cout << "Vlabel & ADJ LIST:\n";
  cout << *(p._matrix);
/*
  if (p._sup_ok == true) { 
    cout << "Support available, value:" << p._vat.size() << endl;
  }
  else { 
    cout << "Support NOT available yet\n"; 
  }
*/
  return ostr;
}

#endif
