#ifndef _FACTORY_HELPER_H_
#define _FACTORY_HELPER_H_

#include <map>
#include "helper_funs.h"

// it remembers from which vertex, what is the other end V_T val
// with an edge E_T has tried and failed
template<typename V_T, typename E_T>
struct failed_map
{
  typedef set<E_T> EDG_L;
  typedef typename EDG_L::iterator EIT;
  typedef typename EDG_L::const_iterator CEIT;
  typedef map<V_T, EDG_L > INSIDE_MAP;
  typedef typename map<V_T, EDG_L >::iterator MIT;
  typedef typename map<V_T, EDG_L >::const_iterator CMIT;
  typedef map<int, INSIDE_MAP > OUTSIDE_MAP;
  typedef typename OUTSIDE_MAP::iterator IT;
  typedef typename OUTSIDE_MAP::const_iterator CIT;
  OUTSIDE_MAP _fm;

  void print() const {
    CIT cit = _fm.begin();
    for (; cit != _fm.end(); cit++) {
      cout << cit->first << ":\n"; 
      const INSIDE_MAP& im = cit->second; 
      CMIT cmit = im.begin();
      for (; cmit != im.end(); cmit++) {
        cout << "(" << cmit->first << "==> "; 
	CEIT ceit = cmit->second.begin();
        for (; ceit != cmit->second.end(); ceit++)
	  cout << *ceit << " ";
	cout << ")\n";
      }
    }
  }

  void insert(int vid, V_T v_l, E_T e_l) {
    IT it;
    it = _fm.find(vid);
    if (it == _fm.end()) {
      EDG_L e_set;
      e_set.insert(e_l);
      INSIDE_MAP im;
      im.insert(make_pair(v_l, e_set));
      _fm.insert(make_pair(vid,im));
      return;
    }
    else {
      INSIDE_MAP& im = it->second;
      MIT mit = im.find(v_l);
      if (mit == im.end()) {
	EDG_L e_set;
	e_set.insert(e_l);
	im.insert(make_pair(v_l, e_set));
      }
      else {
        EDG_L& edge_set = mit->second;
	pair<EIT, bool> ret = edge_set.insert(e_l);
	if (ret.second == false) {
          cout << "ERROR in failed_map:insert, this edge already present!\n";
	  exit(1);
	}
      }
    }			
  }
  bool exist(int vid, V_T v_l, E_T e_l) const {
    CIT it;
    it = _fm.find(vid);
    if (it == _fm.end())
      return false;
    else {
      const INSIDE_MAP& im = it->second;
      CMIT mit = im.find(v_l);
      if (mit == im.end())
	return false;
      else {
        const EDG_L& edge_set = mit->second;
	CEIT eit = edge_set.find(e_l);
	if (eit == edge_set.end())
	  return false;
	else 
          return true;
      }
    }
  }
};

template<typename V_T, typename E_T>
struct edge_counter 
{
  typedef pair<pair<V_T, V_T>, E_T> EDGE_T;
  typedef typename map<EDGE_T, unsigned int>::iterator IT;
  typedef typename map<EDGE_T, unsigned int>::const_iterator CIT;

  map<EDGE_T, unsigned int> _counter;

  unsigned int get_count(V_T src_l, V_T dest_l, E_T edge_l) {
    EDGE_T e;
    if (src_l < dest_l)
      e = make_pair(make_pair(src_l, dest_l), edge_l);
    else 
      e = make_pair(make_pair(dest_l, src_l), edge_l);
    CIT cit = _counter.find(e);  
    if (cit != _counter.end()) return cit->second;
    else return 0;
  }
  void insert(V_T src_l, V_T dest_l, E_T edge_l) {
    EDGE_T e;
    if (src_l < dest_l)
      e = make_pair(make_pair(src_l, dest_l), edge_l);
    else 
      e = make_pair(make_pair(dest_l, src_l), edge_l);
    IT it = _counter.find(e);
    if (it == _counter.end())
      _counter.insert(make_pair(e, 1));
    else 
      it->second++;
  }
  void print() const {
    CIT cit = _counter.begin();
    for (; cit != _counter.end(); cit++)
      
      cout << "(" << cit->first.first.first << " " << cit->first.second << " " << cit->first.first.second << "):" << cit->second << endl;
  }

  bool operator<(const edge_counter& other) const {
    const map<EDGE_T, unsigned int>& other_cnt = other._counter;
    return (_counter < other_cnt);
  }
};

#endif
