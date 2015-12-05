CC=gcc -m64 -g
CFLAGS=-O3 -Wall -std=c99 -march=native
#CFLAGS=-O0 -Wall -std=c99
INCLUDE=-I /usr/local/include
LIB=-L /usr/local/lib -lgsl -lgslcblas -lm
#CC=icc

OBJS=hash.o main.o
EXEC = hash_test

default: hash

hash.o: hash.c
	$(CC) $(CFLAGS) -c $< -o $@

main.o: main.c
	$(CC) $(CFLAGS) -c $< -o $@

hash: $(OBJS)
	$(CC) $(CFLAGS) $(LIB) $(OBJS) -o $(EXEC)
