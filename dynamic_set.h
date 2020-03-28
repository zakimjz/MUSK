/* id should return an ordered
 * integer number from 0 to n-1
 */


struct obj {
  unsigned int val;
  obj* p;
  unsigned int rank;

  obj(unsigned int id) {
    val = id;
    p = this;
    rank = 0;
  }
};

template <typename T>
class disjoint_set {

  public:
    
    // constructor //    
    disjoint_set(size_t total) {
      ptr_array.resize(total, 0);
      set_count=0;
    }
     
    // destructor //
    ~disjoint_set() {
      for (size_t i=0; i< ptr_array.size(); i++) {
        if (ptr_array[i]) {
          delete ptr_array[i];
        }
        ptr_array[i] = 0;
      }
    }

    T* make_set(unsigned int id) {
      T* new_obj = new obj(id);    
      ptr_array[id] = new_obj;
      set_count++;
      return new_obj;
    }

    void set_union(T* x, T* y) {
      link(find_set(x), find_set(y));
      set_count--;
    }

    T* find_set(T* x) const {
      if (x != x->p) {
        x->p = find_set(x->p);
      }
      return x->p;
    }

    T* id_to_ptr(unsigned int id) const {
      return ptr_array[id];
    }

    bool is_rep(unsigned int id) const {
      if (T* obj_ptr = ptr_array[id]) {
        if (find_set(obj_ptr) == obj_ptr) {return true;}
      }
      return false;
    }

    unsigned int get_set_count() const {
      return set_count;
    }

    void print_set_membership() const {
      int total = ptr_array.size();
      for (int i=0; i < total; i++) {
        T* p = id_to_ptr(i);
        if (p) {
          T* rep = find_set(p);
          cout << rep->val << endl; // bad code, should call << operator of T class
        }
        else {
          cout << total << "\n";
        }
      }
    }
    
  private:
    vector<T*> ptr_array;
    unsigned int set_count;

    void link(T* x, T* y) {
      if (x->rank > y->rank) {
        y->p = x;
      }
      else {
        x->p = y;
        if (x->rank == y->rank) {
          y->rank = y->rank + 1;
        }
      }
    }
};
