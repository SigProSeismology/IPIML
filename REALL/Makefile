#MAC openMP
CC = gcc -Os -fopenmp -W 

#Linux openMP
#CC = gcc -Os  -fopenmp -W 

#Serial
#CC = gcc -Os -W 

LIBS = -lm 

all: REAL clean

BIN = ./

REAL:REAL.o
	$(CC) -o $(BIN)/$@ $^ $(LIBS)

clean:
	rm -r *.o
