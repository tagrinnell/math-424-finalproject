CC = g++
FFLAGS = -Wall -fopenmp
LFLAGS = -lgomp -lm
OBJECTS = main.o
ROOT_DIR=$(pwd)

manual:
	g++ -fopenmp main.cpp -o main.o -c $<
	g++ -fopenmp -lgomp -lm main.o -o main.exe -I.

# mpicxx mpi_testing.cpp -o mpi_test.o && mpiexec -n 3 ./mpi_test.o
mpi_main:
	mpicxx main.cpp -o main.o && mpiexec -np 3 ./main.o

mpitest:
	mpicxx mpi_testing.cpp -o mpi_test.o && mpiexec -n 4 ./mpi_test.o

clean:
	rm -f $(OBJECTS) *.exe *.o
