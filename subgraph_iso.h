#ifndef _SUBGRAPH_ISO_H__
#define _SUBGRAPH_ISO_H__

#include "adj_matrix.h"

using namespace std;

// we want to test whether A is  a subgraph of B, where
// M is the corresponding permutation matrix
bool subgraph_iso_test(const SqrSymMatrix& A, const SqrSymMatrix& B, const Matrix& M) {
#ifdef SUBGRAPH_ISO_DEBUG
  cout << "In final Subgraph iso test:\n";
  cout << "A size:" << A.size() << " and B size:" << B.size() << endl;
#endif
  Matrix E = M * B;
  Matrix D = Transpose (E);
  Matrix C = M * D;
#ifdef SUBGRAPH_ISO_DEBUG
  cout << "C:" << endl;
  cout << C << endl;
  cout << "A:" << endl;
  cout << A << endl;
#endif
  for (int i = 0; i < A.size(); i++) {
    for (int j = 0; j < A.size(); j++) 
      if ((A.at(i,j) == 1)  && (C.at(i,j) == 0)) {
        return 0;
      }
  }
  return 1;
}
// this code will work for undirected graph only
template <typename V_T, typename E_T>
bool UllMan_edge_Refinement(const FullLabelAdjMatrix<V_T, E_T>& A, const FullLabelAdjMatrix<V_T, E_T>& B, Matrix& M, const unsigned int& d, int k) {


  int i=0;

#ifdef SUBGRAPH_ISO_DEBUG
  cout << "In edge refinement" << endl;
#endif
  vector<unsigned int> lst;  
  A.neighbors(d, lst);
  for (int idx1=0;idx1<lst.size(); idx1++) { // for each neighbors of d in A
    int i = lst[idx1];              // i is one neighbor of d in A
    E_T e1 = A.get_edge_label(d,i);
#ifdef SUBGRAPH_ISO_DEBUG
    cout << i << " is one neighbor of " << d << endl;
    cout << "with edge label:" << e1 << endl;
#endif
    vector<unsigned int> match_i;         
    M.neighbors(i, match_i);
    for (int idx2=0; idx2<match_i.size(); idx2++) { // for each matching of i in B
      int j = match_i[idx2]; 
      // cout << j << " is one possible match of " << i << endl;
      if (B.at(k,j) == false) {
        M.set(i,j,0);
        if (M.rowset_empty(i)) {
#ifdef SUBGRAPH_ISO_DEBUG
          cout << "Refined ...\n";
#endif
          return true;
        }
      }
      else {
        E_T e2 = B.get_edge_label(k,j);
        if (e1 != e2) {
#ifdef SUBGRAPH_ISO_DEBUG
          cout << "setting (" << i << "," << j << ") to 0" << endl;
#endif
          M.set(i,j,0);
        }
        if (M.rowset_empty(i)) {
#ifdef SUBGRAPH_ISO_DEBUG
          cout << "Refined ...\n";
#endif
          return true;
        }
      }
    }
  }
#ifdef SUBGRAPH_ISO_DEBUG
  cout << "out from edge refinement" << endl;
#endif
  return false;
}

// this code will work for undirected graph only
bool UllMan_Refinement(const SqrSymMatrix& A, const SqrSymMatrix& B, Matrix& M, int d) {


  int elim = 0;
  int i=0;
  int k=0;
  int h=1;
  
#ifdef SUBGRAPH_ISO_DEBUG
  cout << "In refinement" << endl;
#endif
  for(int i=d; i<A.size(); i++) {
    vector<unsigned int> lst;                  // holds adjacent vertices of i in A
    A.neighbors(i,lst);
    vector<unsigned int> match_i;              // holds the vertices in B that matches i
    M.neighbors(i, match_i);
    for (int idx1=0;idx1<lst.size(); idx1++) { // for each neighbors of i in A
      unsigned int x = lst[idx1];              // x is one neighbor of i in A
      for (int idx2=0; idx2<match_i.size(); idx2++) { // for each matching of i in B
        unsigned int j = match_i[idx2]; 
        boost::dynamic_bitset<> res = M[x] & B[j];
        if (res.none()) {
#ifdef SUBGRAPH_ISO_DEBUG
          cout << "setting (" << i << "," << j << ") to 0" << endl;
#endif
          M.set(i,j,0);
          if (M.rowset_empty(i)) {
#ifdef SUBGRAPH_ISO_DEBUG
            cout << "Refined ...\n";
#endif
            return true;
          }
          elim++;
        }
      }
    }
    if (elim==0) break;
    else {elim=0;} 
  }
#ifdef SUBGRAPH_ISO_DEBUG
  cout << "out from refinement" << endl;
#endif
  return false;
}

void print_vector(const boost::dynamic_bitset<>& v) {
  cout << "[";
  for(int i = 0; i < v.size(); i++) {
    if (i < v.size()-1) cout << v[i] << ",";
    else cout << v[i] << "]\n";
  }
}

template <typename T>
void print_vector(vector<T>& v) {
  cout << "[";
  for(int i = 0; i < v.size(); i++) {
    if (i < v.size()-1) cout << v[i] << ",";
    else cout << v[i] << "]\n";
  }
}
// Main Loop: Ullman backtracking ===============
bool UllMan_backtracking(const SqrSymMatrix& A, const SqrSymMatrix& B, const Matrix& M_minus1, bool all_iso) {

  unsigned int alpha = M_minus1.row();   // size of smaller matrix
  unsigned int beta = M_minus1.col();    // size of larger matrix
 
  for (int i=0; i<alpha; i++) {
    if (M_minus1.rowset_empty(i) == true) {
      // cout << "Returning at the beginning\n";
      return false;
    }
  }

  vector<int> H(alpha, -1);         // H[i] = k means k'th column has been used in depth i
  boost::dynamic_bitset<> F(beta);  // set '1' for used column, 0 for unused (boost bitset default is 0)
  vector<Matrix> M_d(alpha);        // Store matching matrix at different depth

  // initialization (step 1)
  int d = 0;                      
  Matrix M = M_minus1;                 
  if (UllMan_Refinement(A, B, M, d) == true) return false; 
  int k,j;               // column scan at depth d start from k'th column

#ifdef SUBGRAPH_ISO_DEBUG
  cout << "alpha=" << alpha << " and beta=" << beta << endl;
#endif

  M_d[0] = M;
  bool extend_mode = true;
  bool same_row_next = false;
#ifdef SUBGRAPH_ISO_DEBUG
  cout << "start with:M=\n" << M << endl;
#endif
  while (true) {
#ifdef SUBGRAPH_ISO_DEBUG
    cout << "At first line after outer loop" << endl;
#endif
    if (extend_mode == true) {
      k = (d ==0) ? H[0] : -1;
      extend_mode = false;
    }
    if (same_row_next == false) {
#ifdef SUBGRAPH_ISO_DEBUG
      cout << "starting search from " << k << " th column" << endl;    
      cout << "with matrix:" <<endl;
      cout << M;
      cout << endl;
#endif
  
      for (j=k+1; j < beta; j++) {
        if (M.at(d,j)==1 && F[j] == 0) break;   // break when finds a column to choose at this depth
      }
    }
    if (j < beta || same_row_next == true) {    // j found, whose value = k, so k'th column was selected at depth d
#ifdef SUBGRAPH_ISO_DEBUG
      cout << "saving the following matrix at " << d << "th pos\n";
      cout << M << endl;
#endif
      same_row_next = false;
     
      k = j;
      M_d[d] = M; 
      M.reset(d);
      M.set(d,k, 1);   // setting M_(dj)
#ifdef SUBGRAPH_ISO_DEBUG
      cout << "column " << k << " has chosen at depth " << d << ", the matching matrix:\n"; 
      cout << M << endl; 
#endif
      bool fail = UllMan_Refinement(A, B, M, d); 
      if (fail) {
#ifdef SUBGRAPH_ISO_DEBUG
        cout << "Refined" << endl;
#endif
        for (j=k+1; j < beta; j++) {
          if (M_d[d].at(d,j)==1 && F[j] == 0) break;   // break when finds a column to choose at this depth
        }
        if (j < beta) {
          same_row_next = true;
          M = M_d[d];
          continue;
        }
      }
      if (d < alpha-1) {
        H[d] = k;
        F[k] = 1;
        d = d + 1;
        extend_mode = true;
        continue;
      }
      else {
#ifdef SUBGRAPH_ISO_DEBUG
        cout << "Now checking if following is an iso:" << endl;
        cout << M;
#endif
        if (subgraph_iso_test(A, B, M) == true) {
#ifdef SUBGRAPH_ISO_DEBUG
          cout << "YES, an iso\n";
#endif
          // cout << M << endl;
          if (!all_iso) return true;
        }
        else {
#ifdef SUBGRAPH_ISO_DEBUG
          cout << "NO, not an iso\n";
#endif
        }
        if (k < beta-1) {
          M = M_d[d];
          continue;
        }
      }
    }
    else {
#ifdef SUBGRAPH_ISO_DEBUG
      cout << "No column found at depth " << d << endl;
#endif
    }
    if (d == 0) break;
#ifdef SUBGRAPH_ISO_DEBUG
    cout << "setting F[" << k << "]=0\n";
#endif
    if (k>=0)
      F[k] = 0; 
    d = d - 1;
    M = M_d[d];
#ifdef SUBGRAPH_ISO_DEBUG
    cout << "now k=" << H[d] << endl;
#endif
    k = H[d];
    if (k>=0)
      F[k] = 0; 
#ifdef SUBGRAPH_ISO_DEBUG
    cout << "Saved the setting ..." << endl;
    cout << "F:"; print_vector(F);
    cout << "H:"; print_vector(H);
#endif
  }
  return false;
} 

// Main Loop: Ullman backtracking ===============
template <typename V_T, typename E_T>
bool UllMan_backtracking(const FullLabelAdjMatrix<V_T, E_T>& A, const FullLabelAdjMatrix<V_T, E_T>& B, const Matrix& M_minus1, bool all_iso) {

  unsigned int alpha = M_minus1.row();   // size of smaller matrix
  unsigned int beta = M_minus1.col();    // size of larger matrix
 
#ifdef SUBGRAPH_ISO_DEBUG
  cout << "\nsubgraph:" << endl;
  cout << A << endl;
  cout << "\ngraph:" << endl;
  cout << B << endl;
#endif
  for (int i=0; i<alpha; i++) {
    if (M_minus1.rowset_empty(i) == true) {
#ifdef SUBGRAPH_ISO_DEBUG
      cout << "Returning form the beginning with false" << endl;
#endif
      return false;
    }
  }

  vector<int> H(alpha, -1);         // H[i] = k means k'th column has been used in depth i
  boost::dynamic_bitset<> F(beta);  // set '1' for used column, 0 for unused (boost bitset default is 0)
  vector<Matrix> M_d(alpha);        // Store matching matrix at different depth

  // initialization (step 1)
  int d = 0;                      
  Matrix M = M_minus1;                 
  if (UllMan_Refinement(A, B, M, d) == true) return false; 
  int k,j;               // column scan at depth d start from k'th column

#ifdef SUBGRAPH_ISO_DEBUG
  cout << "alpha=" << alpha << " and beta=" << beta << endl;
#endif

  M_d[0] = M;
  bool extend_mode = true;
  bool same_row_next = false;
#ifdef SUBGRAPH_ISO_DEBUG
  cout << "start with:M=\n" << M << endl;
#endif
  while (true) {
#ifdef SUBGRAPH_ISO_DEBUG
    cout << "At first line after outer loop" << endl;
#endif
    if (extend_mode == true) {
#ifdef SUBGRAPH_ISO_DEBUG
      cout << "In extend mode, with d=" << d << " and k=" << k << endl;
#endif
      k = (d ==0) ? H[0] : -1;
      extend_mode = false;
    }
    if (same_row_next == false) {
#ifdef SUBGRAPH_ISO_DEBUG
      cout << "starting search from " << k << " th column" << endl;    
      cout << "with matrix:" <<endl;
      cout << M;
      cout << endl;
#endif
  
      for (j=k+1; j < beta; j++) {
        if (M.at(d,j)==1 && F[j] == 0) break;   // break when finds a column to choose at this depth
      }
    }
    if (j < beta || same_row_next == true) {    // j found, whose value = k, so k'th column was selected at depth d
#ifdef SUBGRAPH_ISO_DEBUG
      cout << "saving the following matrix at " << d << "th pos\n";
      cout << M << endl;
#endif
      same_row_next = false;
     
      k = j;
      M_d[d] = M; 
      M.reset(d);
      M.set(d,k, 1);   // setting M_(dj)
#ifdef SUBGRAPH_ISO_DEBUG
      cout << "column " << k << " has chosen at depth " << d << ", the matching matrix:\n"; 
      cout << M << endl; 
#endif
      bool fail = UllMan_Refinement(A, B, M, d); 
      bool fail2 = false;
      if (!fail) {
        fail2 = UllMan_edge_Refinement(A, B, M, d, k); 
      }
      if (fail || fail2) {
#ifdef SUBGRAPH_ISO_DEBUG
        cout << "Refined" << endl;
#endif
        for (j=k+1; j < beta; j++) {
          if (M_d[d].at(d,j)==1 && F[j] == 0) break;   // break when finds a column to choose at this depth
        }
        if (j < beta) {
#ifdef SUBGRAPH_ISO_DEBUG
          cout << "Refined, but another j=" << j << " is available in the same row:" << endl;
          cout << "overwriting M with the saved copy of M_d[d]:" << endl;
          cout << M_d[d];
          cout << endl;
#endif
          M = M_d[d];
          same_row_next = true;
          continue;
        }
        else {
#ifdef SUBGRAPH_ISO_DEBUG
          cout << "No j found on the same row, should backtrack now" << endl;
#endif
          extend_mode=false;
          k = beta;
          continue;
        } 
      }
      if (d < alpha-1) {
        H[d] = k;
        F[k] = 1;
        d = d + 1;
        extend_mode = true;
        continue;
      }
      else {
#ifdef SUBGRAPH_ISO_DEBUG
        cout << "Now checking if following is an iso:" << endl;
        cout << M;
        cout << "A:" << endl;
        cout << A << endl;
#endif
        if (subgraph_iso_test(A, B, M) == true) {
#ifdef SUBGRAPH_ISO_DEBUG
          cout << "YES, an iso" << endl;
#endif
          if (!all_iso) return true;
        }
        else {
#ifdef SUBGRAPH_ISO_DEBUG
          cout << "NO, not an iso\n";
#endif
        }
        if (k < beta-1) {
          M = M_d[d];
          continue;
        }
      }
    }
    else {
#ifdef SUBGRAPH_ISO_DEBUG
      cout << "No column found at depth " << d << endl;
#endif
    }
    if (d == 0) break;
    if (k>=0 && k<beta) {
#ifdef SUBGRAPH_ISO_DEBUG
      cout << "setting F[" << k << "]=0\n";
#endif
      F[k] = 0; 
    }
    d = d - 1;
    M = M_d[d];
#ifdef SUBGRAPH_ISO_DEBUG
    cout << "now k=" << H[d] << endl;
#endif
    k = H[d];
    if (k>=0)
      F[k] = 0; 
#ifdef SUBGRAPH_ISO_DEBUG
    cout << "Saved the setting ..." << endl;
    cout << "F:"; print_vector(F);
    cout << "H:"; print_vector(H);
#endif
  }
  return false;
}  

#endif
