CC = /usr/local/intel/compiler90/bin/icc
CXX = /usr/local/intel/compiler90/bin/icpc
FC = /usr/local/intel/compiler90/bin/ifort 

mars2: mars-c.o mars-cpp.o main.o
	$(FC) mars-c.o mars-cpp.o main.o -o mars2 -cxxlib

clean:
	rm -f *.o mars2
mars-c.o: mars.c
	$(CC) -c mars.c -o mars-c.o

mars-cpp.o: mars.cpp
	$(CXX) -c mars.cpp -o mars-cpp.o

main.o: main.f90
	$(FC) -c main.f90 -o main.o
