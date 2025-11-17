CC = g++
FFLAGS = -Wall -fopenmp
LFLAGS = -lgomp -lm
OBJECTS = main.o
ROOT_DIR=$(pwd)

%.o: %.c
# 	echo $ROOT_DIR
	$(CC) -c -fopenmp $<

manual:
	g++ -fopenmp main.cpp -o main.o -c $<
	g++ -fopenmp -lgomp -lm main.o -o main.exe

main.exe: $(OBJECTS)
	$(CC) $(LFLAGS) $(OBJECTS) $(FFLAGS) -o main.exe

clean:
	rm -f $(OBJECTS) *.exe *.o
