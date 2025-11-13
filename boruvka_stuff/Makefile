CC = g++
FFLAGS = -Wall -fopenmp
LFLAGS = -lgomp -lm
OBJECTS =  main.o
# OBJECTS =  boruvka.o main.o
# ROOT_DIR=$(pwd)

main.exe: $(OBJECTS)
	$(CC) $(LFLAGS) $(OBJECTS) $(FFLAGS) -o main.exe

main.o:
# 	g++ -Wall -fopenmp main.cpp -o main.o -c
	g++ -fopenmp main.cpp -o main.o -c

boruvka.o:
# 	g++ -Wall -fopenmp boruvka.cpp -o boruvka.o -c
	g++ -fopenmp boruvka.cpp -o boruvka.o -c


# %.o: %.c
# 	$(CC) $(FFLAGS) -c $<

clean:
	rm -f $(OBJECTS) *.exe *.o
