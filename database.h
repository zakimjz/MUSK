// This class stores a set of input graphs as matrix format
#ifndef _DATABASE_H_
#define _DATABASE_H_

#include <fstream>
#include <exception>
#include <iostream>
#include <map>
#include "pattern.h"
#include "element_parser.h"
#include "StringTokenizer.h"
#include "subgraph_iso.h"
#include "graph_iso_check.h"


class MinSupException: public exception
{
  virtual const char* what() const throw() {
    return "MinSupException: minimum suppport is greater than database size; No pattern is frequent\n";
  }
} min_sup_ex;

class FileOpenErrorException: public exception
{
  virtual const char* what() const throw() {
     return "FileOpenErrorException: input file could not be opened\n";
  }
} file_open_ex;

class FileFormatException: public exception
{
  virtual const char* what() const throw() {
    return "FileFormatException: input file does not has the desired format\n";
  }
} file_format_ex;

// class to sort an array to get sorted-index
template<class T> struct index_cmp {
index_cmp(const T arr) : arr(arr) {}
  bool operator()(const size_t a, const size_t b) const {
    return arr[a] < arr[b];
  }
  const T arr;
};

template <typename PAT_T>
class Database {
  public:
    typedef typename PAT_T::VERTEX_T V_T;
    typedef typename PAT_T::EDGE_T E_T;
    typedef typename PAT_T::ST MAT_T;
    typedef typename PAT_T::NEIGHBORS NEIGHBORS;
    typedef typename PAT_T::NCIT NCIT;
    typedef map<V_T, NEIGHBORS > EDGE_MAP;
    typedef typename EDGE_MAP::iterator MAPIT;
    typedef typename EDGE_MAP::const_iterator MAPCIT;
    typedef typename PAT_T::EDGE EDGE;
    typedef map<EDGE, int > EDGE_FREQ_MAP; 
    typedef map<EDGE, pair<vector<int>, int > > EDGE_INFO_MAP; // vector store which transaction it occurs,
                                                                // int stores its maximum freq. in any graph
    typedef typename EDGE_INFO_MAP::iterator INFO_IT;
    typedef typename EDGE_INFO_MAP::const_iterator INFO_CIT;
    typedef typename EDGE_FREQ_MAP::iterator FREQ_MAP_IT;
    typedef typename EDGE_FREQ_MAP::const_iterator FREQ_MAP_CIT;
    typedef typename PAT_T::EDGE_CIT EDGE_CIT;
    typedef typename PAT_T::EDGE_IT EDGE_IT;

    Database(const char* filename) {
      ifstream infile(filename, ios::in);
      if (!infile) {
        throw file_open_ex;
      }
      int pos = 0;
      int graph_no = -1;
      while (1) {
        EDGE_FREQ_MAP local_map;
        int ret_val = read_next(infile, pos, graph_no, local_map); 
        vat_and_freq_update(local_map, graph_no);
        if (ret_val == -1) break;
      }
      cout << "total graph in database:" << _graph_store.size() << endl;
    }

    const NEIGHBORS& get_neighbors(V_T src_v) const {
      MAPCIT cit = _ext_map.find(src_v);
      return cit->second;
    }

    int get_freq(EDGE e) {
      INFO_CIT cit = _edge_info.find(e);
      if (cit != _edge_info.end()) {
        return cit->second.second;
      }
      return -1;
    }

    EDGE get_a_freq_edge() const {
      INFO_CIT cit = _edge_info.begin();
      return cit->first;
    }

    const EDGE_INFO_MAP& get_all_edge_info() const {
      return _edge_info;
    }
    const vector<int>& get_edge_vat(EDGE e) const {
      INFO_CIT cit = _edge_info.find(e);
      if (cit != _edge_info.end()) {
        return cit->second.first;
      }
      return Database<PAT_T>::_no_data;
    }

    void remove_infrequent_edges() {
      // identifying infrquent edges
      INFO_IT cit=_edge_info.begin();
      for (; cit != _edge_info.end();) {
        if (cit->second.first.size() < _minsup) {
          EDGE e = cit->first;
          V_T vl1 = e.first.first;
          V_T vl2 = e.first.second;
          E_T el = e.second;
          _edge_info.erase(cit++);
          MAPIT it = _ext_map.find(vl1);
          if (it == _ext_map.end()) {
            cout << "ERROR: this v-label should have been found\n";
            exit(1);
          }
          NEIGHBORS& ngbr = it->second; 
          int removed_cnt = ngbr.erase(make_pair(vl2, el)); 
          if (removed_cnt == 0) {
            cout << "ERROR: this edge (" << vl1 << "," << vl2 << "," << el << ") should have been found\n";
            exit(1);
          }
          if (vl1 == vl2) continue;
          it = _ext_map.find(vl2);
          if (it == _ext_map.end()) {
            cout << "ERROR: this v-label should have been found\n";
            exit(1);
          }
          NEIGHBORS& ngbr2 = it->second; 
          removed_cnt = ngbr2.erase(make_pair(vl1, el)); 
          if (removed_cnt == 0) {
            cout << "ERROR: this edge (" << vl1 << "," << vl2 << "," << el << ") should have been found\n";
            exit(1);
          }
        }
        else {
          ++cit;
        }
      }
    }

    int read_next(ifstream& infile, int& pos, int& graph_no, EDGE_FREQ_MAP& local_map) {
      std::string one_line;

      vector<V_T> vlabels;
      ExPattern<V_T, E_T>* pat = 0;
      
      infile.seekg(pos);
      int id = -1;
      while (1) {
        pos = infile.tellg();
        std::getline(infile, one_line);
        // cout << one_line << endl;
        if (one_line.length()<1) {
          _graph_store.push_back(pat);
          return -1;
        }

        if (one_line.at(0) == 't') {  // second time it enters here and exit from here
          if (id > -1) { // got a new graph
            _graph_store.push_back(pat);
            return 1;
          }
          else {  // 1st time it enters here, and set the graph_id
            id = ++graph_no;
          }
        }
        else if (one_line.at(0) == 'v') {
          StringTokenizer strtok = StringTokenizer(one_line, " ");
                    
          int numToks = strtok.countTokens();
          if (numToks != 3) {
            throw file_format_ex;
          }
          strtok.nextToken(); // skipping 1st token
          int vid = strtok.nextIntToken();
          V_T v_lbl = _vl_prsr.parse_element(strtok.nextToken().c_str());
          vlabels.push_back(v_lbl);
        }
        else if (one_line.at(0) == 'e') {
          if (pat == 0) {
            pat = new ExPattern<V_T, E_T>(vlabels);
          }
          StringTokenizer strtok = StringTokenizer(one_line, " ");
                    
          int numToks = strtok.countTokens();
          if (numToks != 4) {
            throw file_format_ex;
          }
          strtok.nextToken(); // skipping 1st token
          int vid1 = strtok.nextIntToken();
          int vid2 = strtok.nextIntToken();
          E_T e_lbl = _el_prsr.parse_element(strtok.nextToken());
          pat->add_edge(vid1,vid2,e_lbl);
          V_T vl1 = pat->label(vid1);
          V_T vl2 = pat->label(vid2);
          ext_map_insert(vl1, vl2, e_lbl);
          EDGE e = (vl1<vl2)? make_pair(make_pair(vl1, vl2), e_lbl) : make_pair(make_pair(vl2, vl1), e_lbl);
          pair<FREQ_MAP_IT, bool> r = local_map.insert(make_pair(e, 1)); 
          if (r.second == false) {
            r.first->second++;
          }
        }
      }
      return -1;
    }      

    void print_database() const {
      for (int i=0; i<_graph_store.size(); i++) {
        cout << *(_graph_store[i]);
      }
    }

    void vat_and_freq_update(EDGE_FREQ_MAP local, int graph_no){
      FREQ_MAP_CIT cit;
      for (cit = local.begin(); cit != local.end(); cit++) {
        INFO_IT ret_it = _edge_info.find(cit->first);
        if (ret_it == _edge_info.end()) {
          vector<int> vat; vat.push_back(graph_no);
          _edge_info.insert(ret_it, make_pair(cit->first, make_pair(vat, cit->second)));
        }
        else {  // info about this edge exist
           ret_it->second.first.push_back(graph_no);
           if (cit->second > ret_it->second.second) {
             ret_it->second.second = cit->second;
           }
        }
      }
    }

    void ext_map_insert(V_T vl1, V_T vl2, E_T el) {
      MAPIT it;
      if ((it = _ext_map.find(vl1)) == _ext_map.end()) {
        NEIGHBORS nbrs;
        nbrs.insert(make_pair(vl2, el)); 
        _ext_map.insert(make_pair(vl1, nbrs));
      }
      else {
        it->second.insert(make_pair(vl2, el));
      }
      if (vl1 == vl2) return;
      // now doing the same for the other end (NOTE: for directed you need to do only one)
      if ((it = _ext_map.find(vl2)) == _ext_map.end()) {
        NEIGHBORS nbrs;
        nbrs.insert(make_pair(vl1, el)); 
        _ext_map.insert(make_pair(vl2, nbrs));
      }
      else {
        it->second.insert(make_pair(vl1, el));
      }
    }

    // In this case, the current VAT is a VAT of a super-pattern
    bool get_exact_sup_from_super_pat_vat(PAT_T* pat) {
#ifdef PRINT_BACKEND
      cout << "In get_exact_sup_from_super_pat_vat" << endl;
      cout << "getting support for this pattern:" << endl;
      cout << *pat << endl;
#endif
      const multiset<EDGE>& mset = pat->get_edgeset();
      EDGE_IT cit = mset.begin();
      EDGE prev = *cit;
      vector<int> sup_list = get_edge_vat(*cit);
      cit++;
      vector<int> out_list;
      for(; cit != mset.end();cit++) {
        // cout << "Inside this for loop:" << endl;
        if (prev != *cit) {
          prev = *cit;
          const vector<int>& its_vat = get_edge_vat(*cit);
          vat_join(its_vat, sup_list, out_list);
          sup_list = out_list; 
          out_list.clear();
        }
      }
#ifdef PRINT_BACKEND
      cout << "Multiset VAT:\n";
      std::copy(sup_list.begin(), sup_list.end(), ostream_iterator<int>(cout, " "));
      cout << endl;
      cout << "Size of multiSet vat:" << sup_list.size() << endl;
#endif
      const vector<int>& cur_vat = pat->get_vat();
#ifdef PRINT_BACKEND
      cout << "VAT before:" << endl;
      std::copy(cur_vat.begin(), cur_vat.end(), ostream_iterator<int>(cout, " "));
      cout << endl;
#endif
      vector<int>::const_iterator it, it2,it3;
      it2 = cur_vat.begin();
      for (it=sup_list.begin(); it <sup_list.end();it++) {
        it3=lower_bound(it2, cur_vat.end(), *it);
        if (it3 == cur_vat.end() || *it3 != *it) {  // not exists
          PAT_T* database_pat = _graph_store[*it];
          if (database_pat->is_super_pattern(pat) == true) {
            out_list.push_back(*it); 
          }
        }
        it2=it3;
      }
      for (int i=0; i<out_list.size(); i++) {
        PAT_T* database_pat = _graph_store[out_list[i]];
        Matrix m(pat->size(), database_pat->size());
        // cout << "Database graph No:" << out_list[i] << endl;
        //cout << *database_pat << endl;
        matcher(*(pat->get_adj_matrix()), *(database_pat->get_adj_matrix()), m);
        bool ret_val = UllMan_backtracking(*(pat->get_adj_matrix()), *(database_pat->get_adj_matrix()), 
                       m, false);
        if (ret_val == true){
          // cout << "adding " << out_list[i] << " to vatlist" << endl;
          pat->add_tid_to_vat(out_list[i]);
        }
      }
      const vector<int>& c_vat = pat->get_vat();
#ifdef PRINT_BACKEND
      cout << "VAT after:" << endl;
      std::copy(c_vat.begin(), c_vat.end(), ostream_iterator<int>(cout, " "));
      cout << endl;
#endif
      pat->set_sup_ok();
#ifdef PRINT_BACKEND
      cout << "Returing from get_exact_sup_from_super_pat_vat" << endl;
#endif
#ifdef DEBUG
      const vector<int> & t = pat->get_vat();
      // cout << "verifying while going down ...";
      vector<int>::const_iterator cit2 = t.begin();
      for (; cit2<t.end();cit2++) {
        PAT_T* database_pat = _graph_store[*cit2];
        Matrix m(pat->size(), database_pat->size());
        matcher(*(pat->get_adj_matrix()), *(database_pat->get_adj_matrix()), m);
        bool ret_val = UllMan_backtracking(*(pat->get_adj_matrix()), *(database_pat->get_adj_matrix()), 
                       m, false);
        if (ret_val == false) {
          cout << *cit2;
          cout << "ERROR" << endl;
          exit(1);
        }
      }
      cout << "done 1" << endl;

#endif

      if (c_vat.size() >= _minsup) {
        pat->set_freq();
        return true;
      }
      else {
        pat->set_status_known();
        return false;
      }
    }

    void verify_vat(PAT_T* pat) {
     const vector<int> & t = pat->get_vat();
      vector<int>::const_iterator cit2 = t.begin();
      for (; cit2<t.end();cit2++) {
        //  cout << "verifying:" << *cit2 << endl;
        PAT_T* database_pat = _graph_store[*cit2];
        Matrix m(pat->size(), database_pat->size());
        matcher(*(pat->get_adj_matrix()), *(database_pat->get_adj_matrix()), m);
        bool ret_val = UllMan_backtracking(*(pat->get_adj_matrix()), *(database_pat->get_adj_matrix()),
                       m, false);
        if (ret_val == false) {
          cout << *cit2;
          cout << "ERROR" << endl;
          exit(1);
        }
        else {
          // cout << "good: with the above matching matrix:" << endl;
          //cout << m << endl;
        }
      }
      // cout << "done 2" << endl;
    }

    static void vat_join(const vector<int>& v1, const vector<int>& v2, vector<int>& out_list) {
      int i=0,j=0;
      while (i<v1.size() && j<v2.size()) {
        if (v1[i] < v2[j]) {
          i++; 
        }
        else if(v2[j] < v1[i]) { 
          j++;
        }
        else {
          out_list.push_back(v1[i]);
          i++;j++;
        }
      }
    }

    // this version computes the exact vat for both frquent and
    // infrequent pattern; use this version for MARGIN lattice-space
    bool get_exact_sup(PAT_T* pat) { 
      vector<int> sup_list;
      const vector<int>& its_vat = pat->get_vat();
      vector<int>::const_iterator it;
      for (it=its_vat.begin(); it <its_vat.end();it++) {
        PAT_T* database_pat = _graph_store[*it];
        if (database_pat->is_super_pattern(pat) == false)  continue;
        sup_list.push_back(*it); 
      }
      
      int max_sup_possible = sup_list.size();
#ifdef PRINT_BACKEND
      cout << "After edge-multiset join, support is:" << max_sup_possible << endl;
#endif

      vector<int> temp;
      temp.reserve(sup_list.size());
      
      // cout << "Pattern graph:" << endl;
      // cout << *pat << endl;
      for (int i=0; i<max_sup_possible; i++) {
        PAT_T* database_pat = _graph_store[sup_list[i]];
        Matrix m(pat->size(), database_pat->size());
#ifdef PRINT_BACKEND
        cout << "Database graph No:" << sup_list[i] << ":";
#endif
        matcher(*(pat->get_adj_matrix()), *(database_pat->get_adj_matrix()), m);
        bool ret_val = UllMan_backtracking(*(pat->get_adj_matrix()), *(database_pat->get_adj_matrix()), 
                       m, false);
        if (ret_val == true)  {
          temp.push_back(sup_list[i]);  
#ifdef PRINT_BACKEND
          cout << "1\n";
#endif
        }
        else {
#ifdef PRINT_BACKEND
          cout << "0\n";
#endif
        }
      }
      pat->set_vat(temp); 
      pat->set_sup_ok();
#ifdef PRINT_BACKEND
      cout << "Support=" << temp.size() << "\n";
#endif
      if (temp.size() >= _minsup) {
        pat->set_freq();
        return true;
      }
      else {
        pat->set_status_known();
        return false;
      }
    }

    // this version exits immediately when it knows that the given pattern would not be frequent
    // so, for infrequent pattern the vat is dirty, however for frequent pattern the vat is clean
    bool get_exact_sup_optimal(PAT_T* pat) { 
      vector<int> sup_list;
      const vector<int>& its_vat = pat->get_vat();
      vector<int>::const_iterator it;
      for (it=its_vat.begin(); it <its_vat.end();it++) {
        PAT_T* database_pat = _graph_store[*it];
        if (database_pat->is_super_pattern(pat) == false)  continue;
        sup_list.push_back(*it); 
      }
      
      int max_sup_possible = sup_list.size();
      if (max_sup_possible < _minsup) return false;
 
#ifdef PRINT
      cout << "After edge-multiset join, support is:" << max_sup_possible << endl;
#endif

      int sup_require = _minsup;
      vector<int> temp;
      temp.reserve(sup_list.size());
      
      // cout << "Pattern graph:" << endl;
      // cout << *pat << endl;
      for (int i=0; i<max_sup_possible; i++) {
        PAT_T* database_pat = _graph_store[sup_list[i]];
        Matrix m(pat->size(), database_pat->size());
#ifdef PRINT
        cout << "Database graph No:" << sup_list[i] << endl;
#endif
        matcher(*(pat->get_adj_matrix()), *(database_pat->get_adj_matrix()), m);
        bool ret_val = UllMan_backtracking(*(pat->get_adj_matrix()), *(database_pat->get_adj_matrix()), 
                       m, false);
        if (ret_val == false)  {
          // cout << "This graph is not a support\n";
          int t = max_sup_possible-1-i+temp.size();
          if (t<_minsup) {
#ifdef PRINT
            cout << "Maximum support possible is:" << t << " which is less than minsup\n";
#endif
            return false;
          }
        }
        else {
          temp.push_back(sup_list[i]);  
        }
      }
      pat->set_vat(temp); 
      pat->set_sup_ok();
      pat->set_freq();
#ifdef PRINT_BACKEND
      cout << "Support=" << temp.size() << "\n";
#endif

#ifdef DEBUG
      const vector<int> & t = pat->get_vat();
      cout << "verifying while going up....\n";
      vector<int>::const_iterator cit = t.begin();
      for (; cit<t.end();cit++) {
        PAT_T* database_pat = _graph_store[*cit];
        Matrix m(pat->size(), database_pat->size());
        matcher(*(pat->get_adj_matrix()), *(database_pat->get_adj_matrix()), m);
        bool ret_val = UllMan_backtracking(*(pat->get_adj_matrix()), *(database_pat->get_adj_matrix()), 
                       m, false);
        if (ret_val == false) {
          cout << *cit;
          cout << "ERROR" << endl;
          exit(1);
        }
      }
      cout << "done 3" << endl;

#endif
      return true;
    }

    static vector<int> set_static_data() {
      vector<int> temp;
      return temp;
    }

    int get_minsup() const {return _minsup;}
    void set_minsup(int minsup) {
      if (minsup > _graph_store.size()) {
        throw min_sup_ex;
      }
      _minsup =minsup;
    }

    // This routine compares two patterns in DB and returns
    // true if the second is a sub-pattern of the first
    // if any of these two is null pattern, it returns false
    bool compare_patterns_in_db(int i, int j) const {
      //cout << "In compare_patterns_in_db" << endl;
      PAT_T* p = _graph_store[i];
      PAT_T* sub_p = _graph_store[j];
      if (p == 0 || sub_p == 0) return false;
      //cout << *sub_p;
      Matrix m(sub_p->size(), p->size());
      matcher(*(sub_p->get_adj_matrix()), *(p->get_adj_matrix()), m);
      //cout << endl << m;
      bool ret_val = UllMan_backtracking(*(sub_p->get_adj_matrix()), *(p->get_adj_matrix()),m, false);
      return ret_val; 
    }

    // This routine can be used to get the maximal patterns from the databases
    // it compares the database patterns to get all the maximal patterns, by
    // checking each against the others.
    void find_max() const {
      // for every non-maximal pattern, the vector below will keep
      // an witness, which is a super-pattern of this non-maximal pattern
      int cnt = _graph_store.size();
      vector<int> witness(cnt, -1); // initializing the witness to -1, at the end the patterns
                                    // that would have witness value=-1 are maximal patterns
      vector<int> size_vec(cnt,0);
      vector<int> sorted_index(cnt);
      for (int i=0; i<cnt; i++) {
        size_vec[i] = _graph_store[i]->edge_cnt();
        sorted_index[i] = i;
      }
      // now finding a sorted index (ascending order) based on the size
      sort(sorted_index.begin(), sorted_index.end(), index_cmp<vector<int>&>(size_vec));
      // for (int i=0; i<cnt; i++) {
        // cout << sorted_index[i] << " " << size_vec[sorted_index[i]] << endl;
      // }
      for (int i=0; i < cnt-1; i++) {
        int x = sorted_index[i];
        // cout << "processing " << x << "(" << size_vec[x] << ")\n";
        PAT_T* sub_p = _graph_store[x];
        for (int j=i+1; j < cnt; j++) {
          int y = sorted_index[j];
          if (_graph_store[x]->edge_cnt() >= _graph_store[y]->edge_cnt()) continue;
          // cout << "comparing with " << y << "(" << size_vec[y] << ")  ";
          PAT_T* p = _graph_store[y];
          Matrix m(sub_p->size(), p->size());
          matcher(*(sub_p->get_adj_matrix()), *(p->get_adj_matrix()), m);
          bool ret_val = UllMan_backtracking(*(sub_p->get_adj_matrix()), *(p->get_adj_matrix()),m, false);
          if (ret_val == true) {
            cout << "pattern(" << x << ")" << *sub_p;
            cout << " is NOT MAXIMAL, because";
            cout << " pattern(" << y << ")" << *p << " is maximal" << endl;
            witness[x] = y;
            cout << endl << m;
            break;
          }
          else {
            // cout << "return value is false\n";
          }
        }
        if (witness[x] > -1)
          cout << witness[x] << "(" << size_vec[witness[x]] << ")" << endl;
        else {
          const typename PAT_T::CAN_CODE& cc = check_isomorphism(sub_p);
          sub_p->set_canonical_code(cc);
          std::string min_dfs_cc = cc.to_string();
          cout << "MAXIMAL:" << min_dfs_cc << endl;
        }
      }
      PAT_T* sub_p = _graph_store[sorted_index[cnt-1]];
      const typename PAT_T::CAN_CODE& cc = check_isomorphism(sub_p);
      sub_p->set_canonical_code(cc);
      std::string min_dfs_cc = cc.to_string();
      cout << "MAXIMAL:" << min_dfs_cc << endl;
      // cout << "Maximal:" << endl << *(_graph_store[sorted_index[cnt-1]]) << endl;
    }

    int size() const { return _graph_store.size();}

  private:
    EDGE_MAP _ext_map;  // it remembers for each label, what are the possible other
                        // label at the other end of an edge (also store the edge label)
    EDGE_INFO_MAP _edge_info;
    vector<PAT_T* > _graph_store;  // store all the graph patterns
    element_parser<V_T> _vl_prsr;
    element_parser<E_T> _el_prsr;
    static vector<int> _no_data;       // dummy vector to return reference to null data
    int _minsup;
};

#endif
