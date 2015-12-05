CC=gcc -m64
#CFLAGS=-O3 -Wall -std=c99 -march=native
CFLAGS=-O0 -g -Wall -std=c99
INCLUDE=-I /usr/local/include
LIB=-L /usr/local/lib -lgsl -lgslcblas -lm
#CC=icc

OBJS=hash.o pair_counts.o main.o
EXEC = hash_test

default: hash

main.o: main.c
	$(CC) $(CFLAGS) -c $< -o $@

pair_counts.o: pair_counts.c
	$(CC) $(CFLAGS) -c $< -o $@

hash.o: hash.c
	$(CC) $(CFLAGS) -c $< -o $@

hash: $(OBJS)
	$(CC) $(CFLAGS) $(LIB) $(OBJS) -o $(EXEC)
