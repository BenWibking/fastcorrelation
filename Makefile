#CC=gcc -m64
CC=icc
CFLAGS=-O3 -Wall -march=native -std=c99
#CFLAGS=-g -Wall -std=c99 -O0
INCLUDE=-I $(HOME)/include
LIB=-L $(HOME)/lib -lgsl -lgslcblas

OBJS_AUTO=hash.o auto_counts.o main.o
OBJS_CROSS=hash.o cross_counts.o main_cross.o
EXEC_AUTO = auto
EXEC_CROSS = cross

default: auto cross

clean:
	rm *.o; rm $(EXEC_AUTO) $(EXEC_CROSS)

main.o: main.c
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

main_cross.o: main_cross.c
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

auto_counts.o: auto_counts.c
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

cross_counts.o: cross_counts.c
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

hash.o: hash.c
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

auto: $(OBJS_AUTO)
	$(CC) $(CFLAGS) $(LIB) $(OBJS_AUTO) -o $(EXEC_AUTO)

cross: $(OBJS_CROSS)
	$(CC) $(CFLAGS) $(LIB) $(OBJS_CROSS) -o $(EXEC_CROSS)
