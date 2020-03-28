#ifndef _PATTERN_FACTORY_
#define _PATTERN_FACTORY_

#include "pattern.h"
#include "database.h"
#include "graph_iso_check.h"


template <typename PAT>
class PatternFactory {

  public:

    typedef typename PAT::VERTEX_T V_T;
    typedef typename PAT::EDGE_T E_T;
    typedef typename PAT::EDGE EDGE;
    typedef Database<PAT> DATABASE;
    typedef typename DATABASE::EDGE_INFO_MAP EDGE_INFO_MAP;
    typedef typename DATABASE::INFO_CIT INFO_CIT;

    static PatternFactory<PAT>* instance(DATABASE* d) {
      if (_instance == 0) {
        _instance = new PatternFactory<PAT>(d);
      }
      return _instance;
    }  

    PAT* make_single_edge_pattern(V_T src, V_T dest, E_T e) {
      vector<V_T> vlabels(2);
      if (src < dest) { 
        vlabels[0] = src;
        vlabels[1] = dest;
      }
      else {
        vlabels[0] = dest;
        vlabels[1] = src;
      }
      PAT* p = new PAT(vlabels);
      p->add_edge(0,1,e);
      EDGE edge = make_pair(make_pair(vlabels[0], vlabels[1]), e);
      const vector<int>& vat = _d->get_edge_vat(edge); 
      p->set_vat(vat);
      p->set_sup_ok();
      int minsup = _d->get_minsup();
      if (vat.size() >= minsup)
        p->set_freq();
      return p;
    }

    // precondition: given pattern is entirely consistent with correct support count
    // and support list
    PAT* pattern_with_edge_removed(PAT* p, const int& a, const int& b) {
#ifdef PRINT
      cout << "Edge(" << a << "," << b << ") are to be removed from:" << endl;
      cout << *p << endl;
      const vector<int>& cur_vat = p->get_vat();
      cout << "its vat size:" << cur_vat.size() << endl;
      std::copy(cur_vat.begin(), cur_vat.end(), ostream_iterator<int>(cout, " "));
      cout << endl;
#endif
      PAT* clone = p->clone();
      EDGE edge = clone->remove_edge(a,b); 

      // now support list needs to be updated
      _d->get_exact_sup_from_super_pat_vat(clone);
#ifdef PRINT
      cout << "done\n";
#endif
      return clone;
    }

    //Alg: take any frequent edge and make a single edge pattern
    PAT* get_one_maximal_pattern() {
      int minsup = _d->get_minsup();
      EDGE edge = _d->get_a_freq_edge();
#ifdef PRINT
      cout << "Making single edge pattern:" << edge.first.first << " " << edge.first.second << " "
           << edge.second << endl;
#endif
      PAT* cand_pat = make_single_edge_pattern(edge.first.first, edge.first.second, edge.second);
      return extend_until_max(cand_pat, minsup);
      cand_pat->set_freq();
    }

    PAT* get_one_frequent_pattern() {
      int minsup = _d->get_minsup();
      EDGE edge = _d->get_a_freq_edge();
      PAT* cand_pat = make_single_edge_pattern(edge.first.first, edge.first.second, edge.second);
      cand_pat->set_status_known();
      return cand_pat;
    }
    int super_pat_count(PAT* pat) {
    
#ifdef PRINT
      cout << "Inside super_pat_count"<< endl;
#endif
      EDGE this_edge;
      int count = 0; 
      for (int vid=0; vid<pat->size(); vid++) {  // extentions from all vertices
      
        V_T src_v=pat->label(vid);  
        const typename DATABASE::NEIGHBORS& nbrs=_d->get_neighbors(src_v);
     
        typename DATABASE::NCIT nit = nbrs.begin();
        while(nit != nbrs.end()) {
          V_T dest_v = nit->first;
          E_T e_label = nit->second;

          if (src_v < dest_v)
            this_edge = make_pair(make_pair(src_v, dest_v), e_label);
          else 
            this_edge = make_pair(make_pair(dest_v, src_v), e_label);
          int frequency = pat->edge_counter(this_edge);
          int max_freq = _d->get_freq(this_edge);
          if (frequency > max_freq) {
            continue;
          }
    
          // now trying all the possible back-edges 
          vector<int> dest_vids;
          pat->get_vids_for_this_label(dest_v, dest_vids);
          vector<int>::iterator vit = dest_vids.begin(); 
          while (vit < dest_vids.end()) {
            if (*vit <= vid || pat->edge_exist(vid, *vit))
              dest_vids.erase(vit);
            else 
    	      vit++;
          }
          count = count + dest_vids.size() + 1;
          nit++;
        }    
      }
#ifdef PRINT
      cout << "Exiting from super_pat_count" << endl;
#endif
    }

    // Check a pattern's staus, is it freq, or maximal compute its VAT also, 
    // for infrequent VAT is INCOMPLETE
    void get_status_optimal(PAT*& pat, bool& is_freq, bool& is_max) {
      int minsup = _d->get_minsup();
      assert(pat->get_sup_ok() == true);
      if (pat->status_known() == true) {
        is_freq=pat->_is_frequent;
        is_max=pat->_is_maximal;
        return;
      }

      PAT* edge=0;
      PAT* cand_pat=0;
    
      EDGE this_edge;
      for (int vid=0; vid<pat->size(); vid++) {  // extentions from all vertices
      
        V_T src_v=pat->label(vid);  
#ifdef PRINT
        cout << "Extending from vertex-id:" << vid << " with label:" << src_v << endl;
#endif    
        const typename DATABASE::NEIGHBORS& nbrs=_d->get_neighbors(src_v);
     
        typename DATABASE::NCIT nit = nbrs.begin();
        while(nit != nbrs.end()) {
          V_T dest_v = nit->first;
          E_T e_label = nit->second;

          if (src_v < dest_v)
            this_edge = make_pair(make_pair(src_v, dest_v), e_label);
          else 
            this_edge = make_pair(make_pair(dest_v, src_v), e_label);
          int frequency = pat->edge_counter(this_edge);
          int max_freq = _d->get_freq(this_edge);
          if (frequency > max_freq) {
#ifdef PRINT
            cout << "Invalid edge: this edge exceed its frequency limit\n";
#endif    
            continue;
          }
#ifdef PRINT
          cout << "found an edge(" << src_v << ","<< dest_v <<"," << e_label<< ") to extend from vid:" << vid << endl;
#endif    
          edge = make_single_edge_pattern(src_v, dest_v, e_label);
    
          // trying all possible extensions from this vid
          // trying all fwd-extension first 
          while (true) {
            cand_pat = pat->clone();
            int lvid = cand_pat->add_vertex(dest_v);
            cand_pat->add_edge(vid, lvid, e_label);
#ifdef PRINT
            cout << "Candidate pattern:\n";
            cout << *cand_pat << endl;
#endif    
            assert(lvid == pat->size());
            assert(cand_pat->size() == lvid+1);
            cand_pat->join_vat(edge);
            if (cand_pat->get_vat().size() < minsup) {
              delete cand_pat;
              break; 
            }
            bool freq = _d->get_exact_sup_optimal(cand_pat);
            if (freq == true) {  // then the pattern is frequent
#ifdef PRINT
              cout << "Frequent, not-max returning...\n";
#endif    
              delete cand_pat;  
              delete edge;
              edge = 0;
              is_max = false;
              is_freq = true;
              pat->set_status_known();
              return;
            }
            else {
              delete cand_pat;
              break;  // no more forward extension from this vertex with this edge label is possible
            }
          } 
          // now trying all the possible back-edges 
          vector<int> dest_vids;
          pat->get_vids_for_this_label(dest_v, dest_vids);
          vector<int>::iterator vit = dest_vids.begin(); 
          while (vit < dest_vids.end()) {
            if (*vit <= vid || pat->edge_exist(vid, *vit))
              dest_vids.erase(vit);
            else 
    	      vit++;
          }
#ifdef PRINT    
          cout << "choices for dest_vid:";
          for (unsigned int i = 0; i < dest_vids.size(); i++) {
            cout << dest_vids[i] << " ";
          }
          cout << "\n";
#endif
          for (vector<int>::iterator it= dest_vids.begin(); it < dest_vids.end(); it++) {
            cand_pat = pat->clone();
            cand_pat->add_edge(vid, *it, e_label);
            cand_pat->join_vat(edge);
            if (cand_pat->get_vat().size() < minsup) {
              delete cand_pat;
              break; 
            }
            bool freq = _d->get_exact_sup_optimal(cand_pat);
            if (freq == true) { // is the pattern frequent?
#ifdef PRINT    
              cout << cand_pat << endl;
#endif
              delete cand_pat;
              is_freq = true;
              is_max = false;
              pat->set_status_known();
              return;
            }
            else {
              delete cand_pat;
            }
          }
          delete edge;
          edge = 0;
          nit++;
        }  
      }
      if (edge) {delete edge; edge=0;}
      pat->set_max();
      is_freq = true;
      is_max= true;
    }
 
    // Check a pattern's staus, is it freq, infreq or maximum
    // we assume that the pattern's VAT is exact
    void get_status(PAT*& pat, bool& is_freq, bool& is_max) {
#ifdef PRINT
      cout << "In get status:"; 
#endif
      int minsup = _d->get_minsup();
      assert(pat->get_sup_ok() == true);
      if (pat->status_known() == true) {
        is_freq=pat->_is_frequent;
        is_max=pat->_is_maximal;
        return;
      }

      PAT* edge=0;
      PAT* cand_pat=0;
    
      EDGE this_edge;
      for (int vid=0; vid<pat->size(); vid++) {  // extentions from all vertices
      
        V_T src_v=pat->label(vid);  
#ifdef PRINT
        cout << "Extending from vertex-id:" << vid << " with label:" << src_v << endl;
#endif    
        const typename DATABASE::NEIGHBORS& nbrs=_d->get_neighbors(src_v);
     
        typename DATABASE::NCIT nit = nbrs.begin();
        while(nit != nbrs.end()) {
          V_T dest_v = nit->first;
          E_T e_label = nit->second;

          if (src_v < dest_v)
            this_edge = make_pair(make_pair(src_v, dest_v), e_label);
          else 
            this_edge = make_pair(make_pair(dest_v, src_v), e_label);
          int frequency = pat->edge_counter(this_edge);
          int max_freq = _d->get_freq(this_edge);
          if (frequency > max_freq) {
#ifdef PRINT
            cout << "Invalid edge: this edge exceed its frequency limit\n";
#endif    
            continue;
          }
#ifdef PRINT
          cout << "found an edge(" << src_v << ","<< dest_v <<"," << e_label<< ") to extend from vid:" << vid << endl;
#endif    
          edge = make_single_edge_pattern(src_v, dest_v, e_label);
    
          // trying all possible extensions from this vid
          // trying all fwd-extension first 
          while (true) {
            cand_pat = pat->clone();
            int lvid = cand_pat->add_vertex(dest_v);
            cand_pat->add_edge(vid, lvid, e_label);
#ifdef PRINT
            cout << "Candidate pattern:\n";
            cout << *cand_pat << endl;
#endif    
            assert(lvid == pat->size());
            assert(cand_pat->size() == lvid+1);
            cand_pat->join_vat(edge);
            if (cand_pat->get_vat().size() < minsup) {
              delete cand_pat;
              break; 
            }
            update_vat(cand_pat);
            if (cand_pat->is_frequent() == true) {  // then the pattern is frequent
#ifdef PRINT
              cout << "Frequent, not-max returning...\n";
#endif    
              delete cand_pat;  
              delete edge;
              edge = 0;
              is_max = false;
              pat->set_status_known();
              cout << "frequent" << endl;
              return;
            }
            else {
              delete cand_pat;
              break;  // no more forward extension from this vertex with this edge label is possible
            }
          } 
          // now trying all the possible back-edges 
          vector<int> dest_vids;
          pat->get_vids_for_this_label(dest_v, dest_vids);
          vector<int>::iterator vit = dest_vids.begin(); 
          while (vit < dest_vids.end()) {
            if (*vit <= vid || pat->edge_exist(vid, *vit))
              dest_vids.erase(vit);
            else 
    	      vit++;
          }
#ifdef PRINT    
          cout << "choices for dest_vid:";
          for (unsigned int i = 0; i < dest_vids.size(); i++) {
            cout << dest_vids[i] << " ";
          }
          cout << "\n";
#endif
          for (vector<int>::iterator it= dest_vids.begin(); it < dest_vids.end(); it++) {
            cand_pat = pat->clone();
            cand_pat->add_edge(vid, *it, e_label);
            cand_pat->join_vat(edge);
            if (cand_pat->get_vat().size() < minsup) {
              delete cand_pat;
              break; 
            }
            update_vat(cand_pat);
            if (cand_pat->is_frequent()) { // is the pattern frequent?
#ifdef PRINT    
              cout << cand_pat << endl;
#endif
              delete cand_pat;
              is_max = false;
              pat->set_status_known();
              return;
            }
            else {
              delete cand_pat;
            }
          }
          delete edge;
          edge = 0;
          nit++;
        }  
      }
      if (edge) {delete edge; edge=0;}
      pat->set_max();
      cout << "MAXIMAL" << endl;
      is_max= true;
    }
   
    void get_sub_patterns(PAT* pat, vector<PAT*>& sub_patterns) {
      if (pat->size() == 0) return;
      if (pat->size() <= 2) {
        PAT* null_pat = pat->make_null_pattern(_d->size());
        sub_patterns.push_back(null_pat);
        return;
      }
      pat->find_removable_edges();
      const vector<pair<int, int> >& re = pat->get_removable_edges();
      vector<pair<int, int> >::const_iterator cit = re.begin();
      for (;cit != re.end(); cit++) {
        PAT* sub_pat = pattern_with_edge_removed(pat, cit->first, cit->second);
        sub_patterns.push_back(sub_pat);
      }
    }

    void get_freq_super_patterns(const PAT* pat, vector<PAT*>& super_patterns) {
#ifdef PRINT
      cout<<"In call to get_all_frequent_super_pattern\n";
#endif
    
      PAT* edge=0;
      PAT* cand_pat=0;
    
      EDGE this_edge;
      typedef map<EDGE, int> EDGE_FREQ;
      typedef typename EDGE_FREQ::const_iterator F_CIT;
      typedef map<std::string, int > ALL_PAT;
      typedef typename ALL_PAT::iterator APIT;
      int minsup = _d->get_minsup();
      
      if (pat->size() == 0) {
        const EDGE_INFO_MAP& eim = _d->get_all_edge_info();
        INFO_CIT cit;
        for (cit = eim.begin(); cit != eim.end(); cit++) {
          PAT* p = make_single_edge_pattern(cit->first.first.first, cit->first.first.second, cit->first.second);
          super_patterns.push_back(p); 
        }
        return;
      }

      for (int vid=0; vid<pat->size(); vid++) {  // extentions from all vertices
      
        V_T src_v=pat->label(vid);  
#ifdef PRINT
        cout << "Extending from vertex-id:" << vid << " with label:" << src_v << endl;
#endif
    
        const typename DATABASE::NEIGHBORS& nbrs=_d->get_neighbors(src_v);
     
        typename DATABASE::NCIT nit = nbrs.begin();
        for(;nit != nbrs.end();nit++) {
          V_T dest_v = nit->first;
          E_T e_label = nit->second;

          if (src_v < dest_v)
            this_edge = make_pair(make_pair(src_v, dest_v), e_label);
          else 
            this_edge = make_pair(make_pair(dest_v, src_v), e_label);
          int frequency = pat->edge_counter(this_edge);
          int max_freq = _d->get_freq(this_edge);
          if (frequency > max_freq) {
            continue;
          }
 
          edge = make_single_edge_pattern(src_v, dest_v, e_label);
    
          // trying fwd-extension from this vertex
          cand_pat = pat->clone();
          int lvid = cand_pat->add_vertex(dest_v);
          cand_pat->add_edge(vid, lvid, e_label);
          cand_pat->join_vat(edge);
          if (cand_pat->get_vat().size() < minsup) {
            delete cand_pat;
          }
          else {
            bool freq = _d->get_exact_sup_optimal(cand_pat);
            if (freq == true) { // is the pattern frequent?
#ifdef PRINT    
              cout << cand_pat << endl;
#endif
              super_patterns.push_back(cand_pat);
            }
            else {
              delete cand_pat;
            }
          }
          // now trying all the possible back-edges 
          vector<int> dest_vids;
          pat->get_vids_for_this_label(dest_v, dest_vids);
          vector<int>::iterator vit = dest_vids.begin(); 
          while (vit < dest_vids.end()) {
            if (*vit <= vid || pat->edge_exist(vid, *vit))
              vit = dest_vids.erase(vit);
            else 
    	      ++vit;
          }
       
          for (vector<int>::iterator it= dest_vids.begin(); it < dest_vids.end(); it++) {
            cand_pat = pat->clone();
            cand_pat->add_edge(vid, *it, e_label);
            cand_pat->join_vat(edge);
            if (cand_pat->get_vat().size() < minsup) {
              delete cand_pat;
            }
            else {
              bool freq = _d->get_exact_sup_optimal(cand_pat);
              if (freq == true) { // is the pattern frequent?
#ifdef PRINT    
                cout << cand_pat << endl;
#endif
                super_patterns.push_back(cand_pat);
              }
              else {
                delete cand_pat;
              }
            }
          } 
          delete edge;
          edge = 0;
        } 
      }
    }

    // get all superpattern (frequent + non-frequent)
    void get_super_patterns(const PAT* pat, vector<PAT*>& super_patterns) {
#ifdef PRINT
      cout<<"In call to get_all_frequent_super_pattern\n";
#endif
    
      PAT* edge=0;
      PAT* cand_pat=0;
    
      EDGE this_edge;
      typedef map<EDGE, int> EDGE_FREQ;
      typedef typename EDGE_FREQ::const_iterator F_CIT;
      typedef map<std::string, int > ALL_PAT;
      typedef typename ALL_PAT::iterator APIT;
    
      
      for (int vid=0; vid<pat->size(); vid++) {  // extentions from all vertices
      
        V_T src_v=pat->label(vid);  
#ifdef PRINT
        cout << "Extending from vertex-id:" << vid << " with label:" << src_v << endl;
#endif
    
        const typename DATABASE::NEIGHBORS& nbrs=_d->get_neighbors(src_v);
     
        typename DATABASE::NCIT nit = nbrs.begin();
        for(;nit != nbrs.end();nit++) {
          V_T dest_v = nit->first;
          E_T e_label = nit->second;

          if (src_v < dest_v)
            this_edge = make_pair(make_pair(src_v, dest_v), e_label);
          else 
            this_edge = make_pair(make_pair(dest_v, src_v), e_label);
          int frequency = pat->edge_counter(this_edge);
          int max_freq = _d->get_freq(this_edge);
          if (frequency > max_freq) {
            continue;
          }
 
          edge = make_single_edge_pattern(src_v, dest_v, e_label);
    
          // trying fwd-extension from this vertex
          cand_pat = pat->clone();
          int lvid = cand_pat->add_vertex(dest_v);
          cand_pat->add_edge(vid, lvid, e_label);
#ifdef MARGIN
          cand_pat->set_last_edge(vid,lvid);
#endif
          cand_pat->join_vat(edge);
          update_vat(cand_pat);
          super_patterns.push_back(cand_pat);

          // now trying all the possible back-edges 
          vector<int> dest_vids;
          pat->get_vids_for_this_label(dest_v, dest_vids);
          vector<int>::iterator vit = dest_vids.begin(); 
          while (vit < dest_vids.end()) {
            if (*vit <= vid || pat->edge_exist(vid, *vit))
              vit = dest_vids.erase(vit);
            else 
    	      ++vit;
          }
       
          for (vector<int>::iterator it= dest_vids.begin(); it < dest_vids.end(); it++) {
            cand_pat = pat->clone();
            cand_pat->add_edge(vid, *it, e_label);
#ifdef MARGIN
          cand_pat->set_last_edge(vid,*it);
#endif
            cand_pat->join_vat(edge);
            update_vat(cand_pat);
            super_patterns.push_back(cand_pat);
          } 
          delete edge;
          edge = 0;
        } 
      }
    }

#ifdef MARGIN
    // get one superpattern (frequent or non-frequent)
    PAT* get_one_super_pattern(const PAT* pat) {
#ifdef PRINT
      cout<<"In call to get_one_super_pattern\n";
#endif
    
      PAT* edge=0;
      PAT* cand_pat=0;
    
      EDGE this_edge;
      typedef map<EDGE, int> EDGE_FREQ;
      typedef typename EDGE_FREQ::const_iterator F_CIT;
      typedef map<std::string, int > ALL_PAT;
      typedef typename ALL_PAT::iterator APIT;
    
      
      for (int vid=0; vid<pat->size(); vid++) {  // extentions from all vertices
      
        V_T src_v=pat->label(vid);  
#ifdef PRINT
        cout << "Extending from vertex-id:" << vid << " with label:" << src_v << endl;
#endif
    
        const typename DATABASE::NEIGHBORS& nbrs=_d->get_neighbors(src_v);
     
        typename DATABASE::NCIT nit = nbrs.begin();
        for(;nit != nbrs.end();nit++) {
          V_T dest_v = nit->first;
          E_T e_label = nit->second;

          if (src_v < dest_v)
            this_edge = make_pair(make_pair(src_v, dest_v), e_label);
          else 
            this_edge = make_pair(make_pair(dest_v, src_v), e_label);
          int frequency = pat->edge_counter(this_edge);
          int max_freq = _d->get_freq(this_edge);
          if (frequency > max_freq) {
            continue;
          }
 
          edge = make_single_edge_pattern(src_v, dest_v, e_label);
    
          // trying fwd-extension from this vertex
          cand_pat = pat->clone();
          int lvid = cand_pat->add_vertex(dest_v);
          cand_pat->add_edge(vid, lvid, e_label);
          cand_pat->set_last_edge(vid,lvid);
          cand_pat->join_vat(edge);
          update_vat(cand_pat);
          delete edge;
          return cand_pat;
        } 
      }
    }

    PAT* get_upper_diamond_pat(PAT* pat1, PAT* pat2) const {
      pair<int, int> last_edge = pat1->get_last_edge();
      if (pat1->size() == pat2->size()) {
        if (pat2->edge_exist(last_edge.first, last_edge.second)) { return 0;} 
      }
      PAT* cand_pat = pat2->clone();
      if (pat1->size() > cand_pat->size()) { // forward edge
        V_T vlabel = pat1->label((last_edge.first > last_edge.second)? last_edge.first : last_edge.second);
        cand_pat->add_vertex(vlabel);
        E_T elabel = pat1->get_edge_label(last_edge.first, last_edge.second);
        cand_pat->add_edge(last_edge.first, last_edge.second, elabel);
      }
      else if (pat1->size() == cand_pat->size()){
        E_T elabel = pat1->get_edge_label(last_edge.first, last_edge.second);
        cand_pat->add_edge(last_edge.first, last_edge.second, elabel);
      }
      update_vat(cand_pat);
      return cand_pat;
    }
#endif

    // extend an SINGLE edge pattern until it is maximal
    PAT* extend_until_max(PAT*& pat, const int& minsup) {
    		
#ifdef PRINT
      cout<<"In call to extend_until_max\n";
#endif 
      PAT* edge=0;
      PAT* cand_pat=0;
    
      typedef map<EDGE, int> EDGE_FREQ;
      typedef typename EDGE_FREQ::const_iterator F_CIT;

      EDGE this_edge;
      
      for (int vid=0; vid<pat->size(); vid++) {  // extentions from all vertices
      
        V_T src_v=pat->label(vid);  
#ifdef PRINT
        cout << "Extending from vertex-id:" << vid << " with label:" << src_v << endl;
#endif 
    
        const typename DATABASE::NEIGHBORS& nbrs=_d->get_neighbors(src_v);
     
        typename DATABASE::NCIT nit = nbrs.begin();
        while(nit != nbrs.end()) {
          V_T dest_v = nit->first;
          E_T e_label = nit->second;
          if (src_v < dest_v)
            this_edge = make_pair(make_pair(src_v, dest_v), e_label);
          else 
            this_edge = make_pair(make_pair(dest_v, src_v), e_label);
          int frequency = pat->edge_counter(this_edge);
          int max_freq = _d->get_freq(this_edge);
          if (frequency > max_freq) {
            cout << "Invalid edge: this edge exceed its frequency limit\n";
            continue;
          }
#ifdef PRINT
          cout << "found an edge(" << src_v << ","<< dest_v <<"," << e_label<< ") to extend from vid:" << vid << endl;
#endif
          edge = make_single_edge_pattern(src_v, dest_v, e_label);
    
          // trying all possible extensions from this vid
          // trying all fwd-extension first 
          while (true) {
            cand_pat = pat->clone();
            int lvid = cand_pat->add_vertex(dest_v);
            cand_pat->add_edge(vid, lvid, e_label);
#ifdef PRINT
            cout << "Candidate pattern:\n";
            cout << *cand_pat << endl;
#endif
            assert(lvid == pat->size());
            assert(cand_pat->size() == lvid+1);
            cand_pat->join_vat(edge);
            if (cand_pat->get_vat().size() < minsup) {
              delete cand_pat;
              break;  // no more forward extension from this vertex with this edge label is possible
            }
            if (update_vat(cand_pat) == true) {  // then the pattern is frequent
              delete pat;  
              pat = cand_pat;
            }
            else {
              delete cand_pat;
              break;  // no more forward extension from this vertex with this edge label is possible
            }
          } 
          // now trying all the possible back-edges 
          vector<int> dest_vids;
          pat->get_vids_for_this_label(dest_v, dest_vids);
          vector<int>::iterator vit = dest_vids.begin(); 
          while (vit < dest_vids.end()) {
            if (*vit <= vid || pat->edge_exist(vid, *vit))
              dest_vids.erase(vit);
            else 
    	      vit++;
          }
    
#ifdef PRINT
          cout << "choices for dest_vid:";
          for (unsigned int i = 0; i < dest_vids.size(); i++) {
            cout << dest_vids[i] << " ";
          }
          cout << "\n";
#endif
        
          for (vector<int>::iterator it= dest_vids.begin(); it < dest_vids.end(); it++) {
            cand_pat = pat->clone();
            cand_pat->add_edge(vid, *it, e_label); 

            // join vat returns a vat which is a super-set of actual vat, 
            // so this vat size is an upper bound of the exact vat size
            cand_pat->join_vat(edge);
            if (cand_pat->get_vat().size() < minsup) {
              delete cand_pat;
              break;  // no more forward extension from this vertex with this edge label is possible
            }
            update_vat(cand_pat);  // update vat computes the exact vat size
            if (cand_pat->is_frequent() == true) { // is the pattern frequent?
              delete pat;
    	      pat = cand_pat;
            }
            else {
              delete cand_pat;
            }
          }
          delete edge;
          edge = 0;
          nit++;
        }  
      }
#ifdef PRINT
      cout << "exiting with the following max sub-graph:\n";
      cout << *pat << endl;
#endif
      pat->set_max();
      return pat;
    }

    static PatternFactory<PAT>* set_static_data() {
      return 0;
    }

    DATABASE* get_database() const {
      return _d;
    }
  protected:
    PatternFactory(DATABASE* d) {
      _d = d;
    }
   
  private:

  // find exact support list after joining 
  bool update_vat(PAT* pat) const {
    bool is_freq = _d->get_exact_sup(pat);
    if (is_freq ==false) return false;
    return true;
  }

  DATABASE* _d;
  static PatternFactory<PAT>* _instance;
};

#endif
