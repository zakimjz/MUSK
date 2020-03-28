# MUSK Algorithm for Uniform Sampling of Frequent Graph Patterns

MUSK samples
the frequent subgraph space and returns the maximal patterns which are
guaranteed to be sampled uniformly. 

**Relevant Publications**

[2009-musk] Mohammad Al Hasan and Mohammed J. Zaki. Musk: uniform sampling of k maximal patterns. In 9th SIAM International Conference on Data Mining. April 2009.

## HOW TO

For input graph file format, see the file GRAPH_large.dat in the dataset directory.

To Compile (you mush have boost installed in your machine):
make

To Run:

    ./main2 -d dataset/GRAPH_large.dat -s min-sup-value -c -how-many-maximal

Below is a run snapshot:

main2 -d dataset/GRAPH_large.dat -s 20 -c 20

        Minimum Support:20 
        total graph in database:1083
        Start Sampling.......
        Line starting with M are maximal, and starting with F are frequent
        Number in the parenthesis are count, for M the count is unique, for F, they are not necessarily unique

        F (1):
        M (1): 0 1 5 1 6:1 2 6 1 5:2 3 5 1 6:3 4 6 1 6:4 5 6 1 6:3 6 6 1 6:1 7 6 1 6:0 8 5 1 6
        F (2):
        F (3):
        F (4):
        F (5):
        F (6):
        M (2): 0 1 5 1 6:1 2 6 1 5:2 3 5 1 6:3 4 6 1 6:4 5 6 1 6:5 6 6 1 6:3 7 6 1 6
        F (7):
        F (8):
        F (9):
        M (3): 0 1 5 1 6:1 2 6 1 5:2 3 5 1 6:3 4 6 1 6:4 5 6 1 5:4 6 6 1 6:3 7 6 1 6:1 8 6 1 6
        F (10):
        F (11):
        F (12):
        F (13):
        F (14):
        F (15):
        F (16):
        F (17):
        F (18):
        F (19):
        F (20):
        F (21):
        F (22):
        F (23):
        F (24):
        F (25):
        F (26):

In the above, we have show part of the output of an example run. The run discovered 3 maximal patterns, 
which are shown in the line starting with M.


**Note:**

For each maximal patterns sampled, the program prints its canonical code, with
a line that proceed with the letter "M". Number in parenthesis after "M" is
the number of unique pattern found so far. 

However, the program also visits a large number of frequent patterns in between
the maximal patterns, and for those patterns, it only output "F (..)" in
each line. The number if parenthesis is number of frequent pattern nodes that
the program has visited. The uniqueness condition is not checked for frequent.

This solution guaranty uniform maximal pattern sampling, but it is generally
inefficient as the walk samples over the entire frequent pattern space. So,
sometimes it takes a long time until it finds the next maximal patterns. The
problem is really severe for very low support, for which number of maximal
patterns is much much smaller than the number of frequent patterns.

To sample maximal patterns more often, you can try to prioritize the sampling
of maximal pattern by controlling the weight ratio between maximal and frequent
pattern (Line 113 and Line 122 random_walk_manager2.h file). 


