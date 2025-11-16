CC = gcc
CFLAGS = -I/mingw64/include
LDFLAGS = -L/mingw64/lib -lopenblas

all: solve

solve: main.o dgesv.o timer.o
	$(CC) -o solve main.o dgesv.o timer.o $(LDFLAGS)

main.o: main.c
	$(CC) $(CFLAGS) -c main.c

dgesv.o: dgesv.c dgesv.h
	$(CC) $(CFLAGS) -c dgesv.c

timer.o: timer.c timer.h
	$(CC) $(CFLAGS) -c timer.c

clean:
	rm -f *.o solve
