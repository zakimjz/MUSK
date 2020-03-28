// This program can find maximal patterns from a collection
// of frequent patterns. The frequent patterns are supplied
// in a file in dmtl format.

#include <iostream>
#include "database.h"

using namespace std;

char* datafile;
int minsup;
int uniq_pat_count;

typedef ExPattern<int, int> PAT;
template<> vector<int> Database<PAT>::_no_data = Database<PAT>::set_static_data();    

void print_usage(char *prog) {
  cerr<<"Usage: "<<prog<<" -d data-file"<<endl;
  exit(0);
}

int main(int argc, char *argv[]) {

  Database<PAT>* database;

  /* creating database and loading data */
  try {
    database = new Database<PAT>(argv[1]);
  }
  catch (exception& e) {
    cout << e.what() << endl;
  }

  /* finding maximal patterns from frequent patterns */
  database->find_max();
}
