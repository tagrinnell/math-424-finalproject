CC = g++
FFLAGS = -Wall -std=c99
LFLAGS = -lgomp -lm
OBJECTS = include/graph.o
ROOT_DIR=$(pwd)

*.o: *.c
	echo $ROOT_DIR
	$(CC) $(FFLAGS) -c $<

main.exe: $(OBJECTS)
	$(CC) $(LFLAGS) $(OBJECTS) -o main.exe

clean:
	rm -f $(OBJECTS) *.exe
