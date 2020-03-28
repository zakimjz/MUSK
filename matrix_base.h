#ifndef _MATRIX_BASE_H_
#define _MATRIX_BASE_H_

#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <cassert>
#include <boost/dynamic_bitset.hpp>

using namespace std;

class Matrix
{

  friend Matrix Transpose(const Matrix& );
  friend Matrix operator*(const Matrix&, const Matrix&);  //matrix mul
  friend ostream& operator<< (ostream& osr, const Matrix&);


  public:
    typedef boost::dynamic_bitset<> BITVECTOR;    
    // Default constructor   
    Matrix();

    // constructor that allocates memory, and initialize
    // all cell of the matrix as 0
    Matrix(size_t r, size_t c);

    const boost::dynamic_bitset<>& operator[](size_t i) const;

    // allocate memory for a matrix of size r X c 
    // initialize all entry with 0
    void allocate(size_t r, size_t c);

    // Return an element of the matrix [r][c]
    bool at(const unsigned int& r, const unsigned int& c) const;

    void reset(const unsigned int& r);
   
    // set the value of the element [r][c]
    void set(const unsigned int& r, const unsigned int& c, bool val);

    inline size_t row() const {return _row;}
    inline size_t col() const {return _col;}

    size_t rowset_cnt(unsigned int i);

    bool rowset_empty(unsigned int i) const;

    void neighbors(const unsigned int& i, vector<unsigned int>& ret_val) const;

  protected:
    vector<BITVECTOR> _data;
  private:
    size_t _row, _col;
};

class SqrSymMatrix :  public Matrix {
  public:
    SqrSymMatrix();

    SqrSymMatrix(size_t n);
 
    int add_vertex();

    void change_adj_matrix(int size, const vector<pair<int, int> >& adj_list);

    void remove_vertex(size_t i);

    bool at(const unsigned int& r, const unsigned int& c) const;

    void neighbors(const unsigned int& i, vector<unsigned int>& ret_val) const;

    bool essential_edge(int src, int dst) const;

    inline size_t size() const {return _size;};

    friend class AdjIterator;
    friend class NodeAdjIterator;
  protected:
   size_t _size;

  private:
 
    bool dfs_visit(int src, BITVECTOR& visited, int dst) const;
};

// Returning each edge only once
class AdjIterator {

  public:
    AdjIterator(const SqrSymMatrix* m);

    void first();

    void next();

    pair<size_t, size_t> current() const;

    bool is_done() const;

  private:
    const SqrSymMatrix* _m;    
    int _i; // row-pos
    int _j; // col-pos
    bool _is_done;
    size_t _max_i;
};

// iterate over all the edge adjacent to a given node
class NodeAdjIterator {

  public:
    NodeAdjIterator(const SqrSymMatrix* m, size_t i);

    void first();

    void next();

    size_t current() const;

    bool is_done() const;

  private:
    const SqrSymMatrix* _m;    
    int _i; // row-pos
    int _j; // col-pos
    bool _is_done;
    size_t _max_i;
};

#endif
