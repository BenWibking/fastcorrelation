CC=gcc -m64
CFLAGS=-g -Wall -std=c99 -O0
#CC=icc
#CFLAGS=-O3 -Wall -march=native -std=c99
INCLUDE=-I $(HOME)/include $(HDF5_C_INCLUDE) $(MPI_CFLAGS)
LIB=-L $(HOME)/lib -lgsl -lgslcblas $(HDF5_C_LIBS)

OBJS_AUTO=hash.o auto_counts.o read_hdf5.o main.o
OBJS_CROSS=hash.o cross_counts.o read_hdf5.o main_cross.o
EXEC_AUTO = auto
EXEC_CROSS = cross

default: auto cross

clean:
	rm *.o; rm $(EXEC_AUTO) $(EXEC_CROSS)

main.o: main.c
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

main_cross.o: main_cross.c
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

read_hdf5.o: read_hdf5.c
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

auto_counts.o: auto_counts.c
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

cross_counts.o: cross_counts.c
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

hash.o: hash.c
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

auto: $(OBJS_AUTO)
	$(CC) $(CFLAGS) $(OBJS_AUTO) $(LIB) -o $(EXEC_AUTO)

cross: $(OBJS_CROSS)
	$(CC) $(CFLAGS) $(OBJS_CROSS) $(LIB) -o $(EXEC_CROSS)
