/* Written by: Mohammad Hasan
 * This program is for debugging purpose. The objective is to validate
 * UllMan's code, by supplying two set of graphs in a single file. A
 * second parameter would be supplied to denote the size of the first
 * set. So, in the file, from the beginning upto the size belongs to the
 * first set and the remaining are the second set.
 *
 * The first set is a database graph and the second set is a set of frequent
 * pattern in that database graph, the program will output the VAT of the
 * frequent graphs.
 */

#include <iostream>
#include "database.h"

using namespace std;

char* datafile;
int minsup;
int uniq_pat_count;

typedef ExPattern<int, int> PAT;
template<> vector<int> Database<PAT>::_no_data = Database<PAT>::set_static_data();    

void print_usage(char *prog) {
  cerr<<"Usage: "<<prog<<" -d data-file -c count -s minsup"<<endl;
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

  int size1 = atoi(argv[2]);
  
  int size2 = (argc < 4)? 0 : atoi(argv[3]);
  int sup;
  for (int j = size1+size2; j < database->size(); j++) {
    sup=0;
    for (int i = 0; i < size1; i++) {
      // cout << "comparing " << i << "and " << j-size1 << endl;
      bool ret_val = database->compare_patterns_in_db(i,j);      
      if (ret_val) {
        cout << i << " ";
        sup++;
      }
    }
    cout << endl;
    cout << "SUP:" << sup << endl;
  }
}
