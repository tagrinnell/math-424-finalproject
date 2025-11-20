CC = g++
FFLAGS = -Wall -fopenmp
LFLAGS = -lgomp -lm
OBJECTS = main.o
ROOT_DIR=$(pwd)

manual:
	g++ -fopenmp main.cpp -o main.o -c $<
	g++ -fopenmp -lgomp -lm main.o -o main.exe -I.

# mpicxx mpi_testing.cpp -o mpi_test.o && mpiexec -n 3 ./mpi_test.o
mpi:
	mpicxx mpi_testing.cpp -o mpi_test.o && mpiexec ./mpi_test.o -n 3
# 	mpicc main.cpp -o main.o

clean:
	rm -f $(OBJECTS) *.exe *.o
