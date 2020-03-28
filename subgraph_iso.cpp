#include <iostream>
#include <sstream>
#include <vector>
#include <cassert>
#include "adj_matrix.h"

using namespace std;

typedef FullLabelAdjMatrix<char, char> SymmetricFullLabelMatrix;
typedef SymmetricFullLabelMatrix::EDGE_T E_T;
// we want to test whether A is  a subgraph of B, where
// M is the corresponding permutation matrix
bool subgraph_iso_test(const SqrSymMatrix& A, const SqrSymMatrix& B, const Matrix& M) {
//bool subgraph_iso_test(AdjMatrix<T>& A, AdjMatrix<T>& B, Matrix& M) {
  Matrix E = M * B;
  Matrix D = Transpose (E);
  Matrix C = M * D;
  for (int i = 0; i < A.row(); i++) 
    for (int j = 0; j < A.col(); j++) 
      if ((A.at(i,j) == 1)  && (C.at(i,j) == 0))
        return 0;
  return 1;
}

bool UllMan_Refinement(const SymmetricFullLabelMatrix& A, const SymmetricFullLabelMatrix& B, Matrix& M) {


  int elim = 0;
  int i=0;
  int k=0;
  int h=1;

  for(int i=0; i<A.size(); i++) {
    vector<unsigned int> lst;                  // holds adjacent vertices of i in A
    A.neighbors(i,lst);
    vector<unsigned int> match_i;              // holds the vertices in B that matches i
    M.neighbors(i, match_i);


    for (int idx1=0;idx1<lst.size(); idx1++) { // for each neighbors of i in A
      unsigned int x = lst[idx1];
      E_T e_label1 = A.get_edge_label(i,x);
      for (int idx2=0; idx2<match_i.size(); idx2++) { // for each matching of i in B
        int j = match_i[idx2]; 
        boost::dynamic_bitset<> res = M[x] & B[j];
        if (res.none()) {
          M.set(i,j,0);
          if (M.rowset_empty(i)) {
#ifdef PRINT
            cout << "Refined ...\n";
#endif
            return false;
          }
          elim++;
        }
        else {
          vector<unsigned int> neighbor_j;              // holds the vertices in B that ar neighbors or j
          B.neighbors(j, neighbor_j);
          int idx3=0;
          for (; idx3<neighbor_j.size(); idx3++) { // for each neighbors or j in B
            if (B.get_edge_label(j, neighbor_j[idx3]) == e_label1) break;
          }
          if (idx3 == neighbor_j.size()) {
            M.set(i,j,0);
            if (M.rowset_empty(i)) {
#ifdef PRINT
              cout << "Refined ...\n";
#endif
              return false;
            }
          }
        }
      }
    }
    if (elim==0) break;
    else {elim=0;} 
  }
  return true;
}


// this code will work for undirected graph only
bool UllMan_Refinement(const SqrSymMatrix& A, const SqrSymMatrix& B, Matrix& M) {


  int elim = 0;
  int i=0;
  int k=0;
  int h=1;

  for(int i=0; i<A.size(); i++) {
    vector<unsigned int> lst;                  // holds adjacent vertices of i in A
    A.neighbors(i,lst);
    vector<unsigned int> match_i;              // holds the vertices in B that matches i
    M.neighbors(i, match_i);


    for (int idx1=0;idx1<lst.size(); idx1++) { // for each neighbors of i in A
      unsigned int x = lst[idx1];
      for (int idx2=0; idx2<match_i.size(); idx2++) { // for each matching of i in B
        unsigned int j = match_i[idx2]; 
        boost::dynamic_bitset<> res = M[x] & B[j];
        if (res.none()) {
          M.set(i,j,0);
          if (M.rowset_empty(i)) {
#ifdef PRINT
            cout << "Refined ...\n";
#endif
            return false;
          }
          elim++;
        }
      }
    }
    if (elim==0) break;
    else {elim=0;} 
  }
  return true;
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
bool UllMan_backtracking(const SqrSymMatrix& A, const SqrSymMatrix& B, const Matrix& M_0, bool all_iso) {

  bool iso_exist = false;
  unsigned int alpha = M_0.row();   // size of smaller matrix
  unsigned int beta = M_0.col();    // size of larger matrix
  vector<int> H(alpha, -1);         // H[i] = k means k'th column has been used in depth i
  boost::dynamic_bitset<> F(beta);  // set '1' for used column, 0 for unused.
  vector<Matrix> M_d(alpha);        // Store matching matrix at different depth

  // initialization
  int d = 0;                      
  Matrix M = M_0;                 
  int k = 0;               // column scan at depth d start from k'th column

#ifdef PRINT
  cout << "alpha=" << alpha << " and beta=" << beta << endl;
#endif
  do {

    int j = k;
    for (; j < beta; j++) 
      if (M.at(d,j)==1 && F[j] == 0) break;   // break when finds a column to choose at this depth
    if (j < beta) {    // j'th column was selected at depth d

      H[d] = j;
      if (k>0 && F[k-1]==1) F[k-1]=0;
      F[j] = 1; 
#ifdef PRINT
      cout << "column " << j << " has chosen at depth " << d << "\n"; 
      cout << "F:" <<  F << endl;
      cout << "H:"; print_vector(H);
#endif
      // for all columns except j of row d, making all entry to be 0 
      M.reset(d);
      M.set(d,j, 1);
      bool not_refined = UllMan_Refinement(A, B, M); 
      if (d < alpha-1 && not_refined) {   // not a leaf
        M_d[d] = M;         // saving the copy of M at depth d
        d++;
#ifdef PRINT
        cout << M;
#endif
        k = H[d]+1;
      }
      else {   // we reached at a leaf node
#ifdef PRINT
        cout << "Now checking if following is an iso:\n" << M;
#endif
        if (not_refined) { 
          if (subgraph_iso_test(A, B, M) == true) {
            iso_exist = true;
            cout << "YES, an iso\n";
            cout << M;
            if (!all_iso) return iso_exist;
          } 
          else {
            cout << "NOT an iso\n";
          }
        }
        else {
#ifdef PRINT
          cout << "This M is refined:" << endl;
          cout << M;
#endif
        }
        M = M_d[d-1];
        F[j] = 0;
        k = j+1;
#ifdef PRINT
        cout << "Search with backtrack\n";
        cout << M;
        cout << "Line 141: column scanning starts from:" << k << " for depth " << d << endl;
#endif
      }
    }
    else {
#ifdef PRINT
      cout << "Line 145: No column found at depth " << d << endl;
#endif
      cout << "H[" << d << "] = " << H[d] << endl;
      if (H[d]>=0) {
        cout << "setting H[" << d<< "] to be -1 " << "and F[" << H[d] << "] to be 0\n"; 
        F[H[d]] = 0;
        H[d] = -1;
      }
#ifdef PRINT
      cout << "F:" << F << endl;
      cout << "H:"; print_vector(H);
#endif
      d--;
      if (d == -1) return iso_exist;
      if (d>0)
        M = M_d[d-1];
      else {
        M = M_0;
        F[H[d]] = 0;
      }
      k = H[d] + 1;
     
#ifdef PRINT
      cout << "M (line 167):\n";
      cout << M;
      cout << "Line 169: column scanning starts from:" << k << " for depth " << d << endl;
#endif
    }

  } while (d > -1);
  
}

bool UllMan_backtracking(const SymmetricFullLabelMatrix& A, const SymmetricFullLabelMatrix& B, const Matrix& M_0, bool all_iso) {

  bool iso_exist = false;
  unsigned int alpha = M_0.row();   // size of smaller matrix
  unsigned int beta = M_0.col();    // size of larger matrix
  vector<int> H(alpha, -1);         // H[i] = k means k'th column has been used in depth i
  boost::dynamic_bitset<> F(beta);  // set '1' for used column, 0 for unused.
  vector<Matrix> M_d(alpha);        // Store matching matrix at different depth

  // initialization
  int d = 0;                      
  Matrix M = M_0;                 
  int k = 0;               // column scan at depth d start from k'th column

#ifdef PRINT
  cout << "alpha=" << alpha << " and beta=" << beta << endl;
#endif
  do {

    int j = k;
    for (; j < beta; j++) 
      if (M.at(d,j)==1 && F[j] == 0) break;   // break when finds a column to choose at this depth
    if (j < beta) {    // j'th column was selected at depth d

      H[d] = j;
      if (k>0 && F[k-1]==1) F[k-1]=0;
      F[j] = 1; 
#ifdef PRINT
      cout << "column " << j << " has chosen at depth " << d << "\n"; 
      cout << "F:" <<  F << endl;
      cout << "H:"; print_vector(H);
#endif
      // for all columns except j of row d, making all entry to be 0 
      M.reset(d);
      M.set(d,j, 1);
      bool not_refined = UllMan_Refinement(A, B, M); 
      if (d < alpha-1 && not_refined) {   // not a leaf
        M_d[d] = M;         // saving the copy of M at depth d
        d++;
#ifdef PRINT
        cout << M;
#endif
        k = H[d]+1;
      }
      else {   // we reached at a leaf node
#ifdef PRINT
        cout << "Now checking if following is an iso:\n" << M;
#endif
        if (not_refined) { 
          if (subgraph_iso_test(A, B, M) == true) {
            iso_exist = true;
            cout << "YES, an iso\n";
            cout << M;
            if (!all_iso) return iso_exist;
          } 
          else {
            cout << "NOT an iso\n";
          }
        }
        else {
#ifdef PRINT
          cout << "This M is refined:" << endl;
          cout << M;
#endif
        }
        M = M_d[d-1];
        F[j] = 0;
        k = j+1;
#ifdef PRINT
        cout << "Search with backtrack\n";
        cout << M;
        cout << "Line 141: column scanning starts from:" << k << " for depth " << d << endl;
#endif
      }
    }
    else {
#ifdef PRINT
      cout << "Line 145: No column found at depth " << d << endl;
#endif
      cout << "H[" << d << "] = " << H[d] << endl;
      if (H[d]>=0) {
        cout << "setting H[" << d<< "] to be -1 " << "and F[" << H[d] << "] to be 0\n"; 
        F[H[d]] = 0;
        H[d] = -1;
      }
#ifdef PRINT
      cout << "F:" << F << endl;
      cout << "H:"; print_vector(H);
#endif
      d--;
      if (d == -1) return iso_exist;
      if (d>0)
        M = M_d[d-1];
      else {
        M = M_0;
        F[H[d]] = 0;
      }
      k = H[d] + 1;
     
#ifdef PRINT
      cout << "M (line 167):\n";
      cout << M;
      cout << "Line 169: column scanning starts from:" << k << " for depth " << d << endl;
#endif
    }

  } while (d > -1);
  
}

int main() {
FullLabelAdjMatrix<char, char> A(4);
A.set_vlabel(0,'A');
A.set_vlabel(1,'A');
A.set_vlabel(2,'A');
A.set_vlabel(3,'A');
A.add_edge(0,1,'a');
A.add_edge(0,2,'b');
A.add_edge(1,2,'c');
A.add_edge(2,3,'a');
cout << "A:\n" << A;

FullLabelAdjMatrix<char, char> B(5);
B.set_vlabel(0, 'A');
B.set_vlabel(1, 'A');
B.set_vlabel(2, 'A');
B.set_vlabel(3, 'A');
B.set_vlabel(4, 'A');
B.add_edge(0,1,'a');
B.add_edge(0,2,'a');
B.add_edge(0,3,'c');
B.add_edge(1,2,'c');
B.add_edge(1,3,'a');
B.add_edge(2,3,'b');
B.add_edge(2,4,'a');
cout << "B:\n" << B;

Matrix M(4,5);
matcher(A, B, M);
cout << "M:\n" << M;
bool iso_exist = UllMan_backtracking(A, B, M, true);
}
