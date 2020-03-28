CC=g++ -O3 -g
OBJ=./StringTokenizer/StringTokenizer.o 
INCLUDE-PATH=-I. -I./StringTokenizer

main2:main2.cpp database.h pattern.h adj_matrix.h pattern_factory.h matrix_base.o subgraph_iso.h graph_iso_check.h random_walk_manager2.h random.o 
	cd StringTokenizer; $(MAKE);
	$(CC) $(INCLUDE-PATH) -o main2 main2.cpp matrix_base.o random.o $(OBJ)

cmp_pat_in_db: cmp_pat_in_db.cpp database.h pattern.h adj_matrix.h matrix_base.o subgraph_iso.h
	$(CC) $(INCLUDE-PATH) -o cmp_pat_in_db cmp_pat_in_db.cpp matrix_base.o $(OBJ)

find_max:find_max.cpp database.h pattern.h adj_matrix.h matrix_base.o subgraph_iso.h
	$(CC) $(INCLUDE-PATH) -o find_max find_max.cpp matrix_base.o $(OBJ)


matrix_base.o: matrix_base.h matrix_base.cpp
	$(CC) -c matrix_base.cpp

random.o: random.h random.cpp
	$(CC) -c random.cpp
clean:
	cd StringTokenizer; $(MAKE) clean;
	rm *.o main2
