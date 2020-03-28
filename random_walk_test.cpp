#include <fstream>
#include "adj_matrix.h"
#include "random.h"
#include "dynamic_set.h"
#include "iterator"
#include "algorithm"

class RandomWalk {

  public:
    typedef AdjMatrix<int> MAT_TYPE;
    void walk(int s_vid, int iter) {
      int cur_iter = 0;
      while (cur_iter < iter) {
        if (_graph->label(s_vid) == 1) {
          cout << s_vid << " ";
          map<int, int>::iterator it = _stat.find(s_vid);
          it->second++;
          cur_iter++;
        }
        int degree = _graph->degree(s_vid);
        int random_num = get_a_random_number(0, degree);
        NodeAdjIterator it(_graph, s_vid);
        it.first();
        for (int i=0; i<random_num; i++) it.next();
        s_vid = it.current();
      }
      map<int, int>::iterator it = _stat.begin();
      for(;it != _stat.end(); it++) {
        cout << it->first << " " << it->second << " " << _graph->degree(it->first) << endl;
      }
    }
    void walk2(int s_vid, int iter, float C) {
      int cur_iter = 0;
      while (cur_iter < iter) {
        if (_graph->label(s_vid) == 1) {
          cout << s_vid << " ";
          map<int, int>::iterator it = _stat.find(s_vid);
          it->second++;
          cur_iter++;
        }
        s_vid = choose_with_bias_prob(C, s_vid);
      }
      map<int, int>::iterator it = _stat.begin();
      for(;it != _stat.end(); it++) {
        cout << it->first << " " << it->second << " " << _graph->degree(it->first) << endl;
      }
    }

    int choose_with_bias_prob(double C, int cur) {
      int degree = _graph->degree(cur);
      vector<double> prob(degree);
      vector<int> nodeidx(degree);
      NodeAdjIterator it(_graph, cur);
      int i=0;
      for (it.first(); !it.is_done(); it.next(),i++) {
        int j= it.current();
        nodeidx[i]=j;
        double weight;
        if (_graph->label(j) == 1) {
          weight = C/_graph->degree(j);
        }
        else {
          weight = C/10.0;
        }
        if (i==0) prob[0] = weight;
        else prob[i] = prob[i-1] + weight;
      }
      //std::copy(prob.begin(), prob.end(), ostream_iterator<double>(cout, " "));
      // cout << endl;
      for (int i=0; i< prob.size(); i++) {
        prob[i] = prob[i]/prob[prob.size()-1];
      }
/*
      std::copy(prob.begin(), prob.end(), ostream_iterator<double>(cout, " "));
      cout << endl;
      exit(1);
*/
      int ret_val = randomWithDiscreteProbability(prob);
      return nodeidx[ret_val];
    }
    void load_graph(int size, const char* filename) {
      _size = size;      
      ifstream infile(filename, ios::in);
      if (!infile) {
        cout << filename << " could not open" << endl;
        exit(1);
      }
      else {
        cout << "Input file opened successfully\n";
      }
      _graph = new MAT_TYPE(size);
      for (int i=0; i<size; i++) {
        int random_num = get_a_random_number(0, 10000);
        if (random_num <50) { 
          _graph->set_vlabel(i, 1);
          _stat.insert(make_pair(i,0));
        }
        else
          _graph->set_vlabel(i,2);
      }
      // now loading the edges
      int a, b;
      while (!infile.eof()) {
        infile >> a >> b; 
        _graph->add_edge(a,b);
      }
      infile.close();
    }

    // increase the degree of few vertex that has label 1
    void process_graph() {
      int num = 5;
      for (int i=0; i<_size; i++) {
        if (_graph->label(i) == 2) continue;
        NodeAdjIterator it(_graph, i);
        for (it.first(); !it.is_done(); it.next()) {
          int j = it.current();
          if (_graph->label(j) == 1) {
            _graph->set(i,j,0);
            _graph->set(j,i,0);
          }
        }
        /*
        int random_num = get_a_random_number(0, 100);
        if (random_num>66) continue;
        for (int j=0; j<num; j++) {
          int other = get_a_random_number(0, _size);
          if (_graph->at(i,j) == 0) {
            _graph->add_edge(i,j);
          }
        }
        */
      }
    }

    // Test the graph for connectedness, if not connected, add all different connected components
    void make_connected() {
      disjoint_set<obj> disj_set(_size); 
      for (int i=0; i<_size; i++) 
        disj_set.make_set(i);
      AdjIterator it(_graph);
      for (it.first(); !it.is_done(); it.next()) {
        pair<int, int> e = it.current();
        obj* first_obj = disj_set.id_to_ptr(e.first);
        obj* second_obj = disj_set.id_to_ptr(e.second);
        if (disj_set.find_set(first_obj) != disj_set.find_set(second_obj))
          disj_set.set_union(first_obj, second_obj);
      }
      cout << "Total set:" << disj_set.get_set_count() << endl;
      for (int i=0; i<_size; i++) {
        NodeAdjIterator  nit(_graph, i);
        cout << _graph->label(i) << "===>";
        for (nit.first(); !nit.is_done(); nit.next()) {
          int j=nit.current(); 
          cout << _graph->label(j) << " ";
        }
        cout << endl;
      }
    }
  private:

    MAT_TYPE * _graph;
    int _size;
    map<int, int> _stat;
};

int main(int argc, char *argv[]) {
  char* filename = argv[1];
  RandomWalk r;
  r.load_graph(10000, filename); 
  r.process_graph();
  //r.make_connected();
  r.walk(0, 20000);
  //r.walk2(0, 20000, 10.0);
}
