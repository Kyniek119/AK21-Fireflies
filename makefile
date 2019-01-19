CC=gcc
CFLAGS=-I. -fopenmp
DEPS = Firefly_func.h Funkcje.h

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) 

Firefly: Firefly.o Firefly_func.o Funkcje.o
	$(CC) -o Firefly.out Firefly.o Funkcje.o Firefly_func.o -fopenmp -lm
