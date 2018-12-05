lsh:
	g++ --std=c++11 -Wall -g cluster.cpp functions.cpp hash.cpp cube.cpp f_utilities.cpp -o cluster -lm

# clean house
clean:
	rm -f cluster.o functions.o hash.o cube.o f_utilities.o cluster

# do a bit of accounting
count:
	wc cluster.cpp functions.cpp hash.cpp cube.cpp f_utilities.cpp
