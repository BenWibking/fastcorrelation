#ifndef __INCLUDE_HASH_H__
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_rng.h>

#define M_PI 3.14159265358979323846

#define testing 0

#define REALLOC_FAC 1.1

#define FLOAT double
#define INDEX(i,j,k) (k + (j + (i*ngrid))*ngrid)
#define MAX(x,y) ((x>y) ? x : y)
#define SQ(x) (x*x)
#define CUBE(x) (x*x*x)
#define PERIODIC(dx) ((dx>0.5*Lbox) ? (dx-Lbox) : dx)

typedef struct {
  double Lbox;
  int ngrid;
  size_t * counts;
  size_t * allocated;
  FLOAT ** x;
  FLOAT ** y;
  FLOAT ** z;
} GHash;


void my_free(void* block);
void* my_malloc(size_t size);
void* my_realloc(void* old_block, size_t new_size, size_t old_size);

GHash* allocate_hash(int ngrid, double Lbox, size_t npoints);
void free_hash(GHash * g);
void geometric_hash(GHash * grid, FLOAT *x, FLOAT *y, FLOAT *z, size_t npoints);

void count_pairs_naive(double * restrict x, double * restrict y, double * restrict z, size_t npoints, long int * restrict pcounts, double * restrict bin_edges_sq, int nbins, double Lbox);
void count_pairs(GHash * restrict g, long int * restrict pcounts, double * restrict bin_edges_sq, int nbins);

#define __INCLUDE_HASH_H__
#endif
