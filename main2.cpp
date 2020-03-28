#include <iostream>
#include "database.h"
#include "random_walk_manager2.h"

using namespace std;

char* datafile;
int minsup;
int pat_count;

typedef ExPattern<int, int> PAT;
template<> vector<int> Database<PAT>::_no_data = Database<PAT>::set_static_data();    
template<> PatternFactory<PAT>* PatternFactory<PAT>::_instance = PatternFactory<PAT>::set_static_data();
//template<> PatternFactory<PAT>::set_static_data();

void print_usage(char *prog) {
  cerr<<"Usage: "<<prog<<" -d data-file -c count -s minsup"<<endl;
  exit(0);
}

void parse_args(int argc, char* argv[]) {
  if(argc<7) {
    print_usage(argv[0]);
  }

  for (int i=1; i < argc; ++i){
    if (strcmp(argv[i], "-d") == 0){
      datafile=argv[++i];
    }
    else if (strcmp(argv[i], "-s") == 0){
      minsup=atoi(argv[++i]);
      cout << "Minimum Support:" << minsup << " " << endl;
    }
    else if(strcmp(argv[i],"-c") == 0){
      pat_count=atoi(argv[++i]);
    }
    else{
      print_usage(argv[0]);
    }
  }
}//end parse_args()


int main(int argc, char *argv[]) {

  parse_args(argc, argv);
  Database<PAT>* database;

  /* creating database and loading data */
  try {
    database = new Database<PAT>(datafile);
    database->set_minsup(minsup);
  }
  catch (exception& e) {
    cout << e.what() << endl;
  }
  database->remove_infrequent_edges();

  /* creating random_walk_manager and starting walk */
  double C=50.0; 
  cout << "Start Sampling.......\n"; 
  cout << "Line starting with M are maximal, and starting with F are frequent" << endl;
  cout << "Number in the parenthesis are count, for M the count is unique, for F, they are not necessarily unique" << endl;;
  RandomWalkManager2<PAT> rwm(C,pat_count,database);
  lattice_node<PAT>* start = rwm.initialize();
  rwm.walk(start);
  cout << "\nWalk terminated:" << endl;
  delete database;
}
