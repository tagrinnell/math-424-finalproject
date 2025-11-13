CC = g++
FFLAGS = -Wall -FOPEN_MP
LFLAGS = -lgomp -lm
OBJECTS = main.o
ROOT_DIR=$(pwd)

%.o: %.c
	echo $ROOT_DIR
	$(CC) $(FFLAGS) -c $<

main.exe: $(OBJECTS)
	$(CC) $(LFLAGS) $(OBJECTS) -o main.exe

clean:
	rm -f $(OBJECTS) *.exe
