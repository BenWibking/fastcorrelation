#CC=gcc -m64
CC=icc
CFLAGS=-O3 -Wall -march=native -vec-report=3 -std=c99
#CFLAGS=-g -Wall -std=c99 -O0
INCLUDE=-I $(HOME)/include
LIB=-L $(HOME)/lib -lgsl -lgslcblas

OBJS=hash.o pair_counts.o main.o
EXEC = hash_test

default: hash

main.o: main.c
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

pair_counts.o: pair_counts.c
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

hash.o: hash.c
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

hash: $(OBJS)
	$(CC) $(CFLAGS) $(LIB) $(OBJS) -o $(EXEC)
