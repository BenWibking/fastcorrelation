#CC=gcc -m64
CC=icc
CFLAGS=-O3 -Wall -march=native -std=c99
#CFLAGS=-g -Wall -std=c99 -O0
INCLUDE=-I $(HOME)/include
LIB=-L $(HOME)/lib -lgsl -lgslcblas

OBJS=hash.o auto_counts.o main.o
EXEC = hash_test

default: hash

main.o: main.c
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

auto_counts.o: auto_counts.c
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

hash.o: hash.c
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

hash: $(OBJS)
	$(CC) $(CFLAGS) $(LIB) $(OBJS) -o $(EXEC)
