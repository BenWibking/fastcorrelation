#ifndef __INCLUDE_HASH_H__
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>

#define REALLOC_FAC 1.1

#define FLOAT double
#define INDEX(i,j,k) (k + (j + (i*ngrid))*ngrid)
#define MAX(x,y) (x>y ? x : y)

typedef struct {
  double Lbox;
  int ngrid;
  size_t * restrict counts;
  size_t * restrict allocated;
  FLOAT ** restrict x;
  FLOAT ** restrict y;
  FLOAT ** restrict z;
} GHash;


GHash* allocate_hash(int ngrid, double Lbox, size_t npoints);
void free_hash(GHash * g);
void geometric_hash(GHash * grid, FLOAT *x, FLOAT *y, FLOAT *z, size_t npoints);

#define __INCLUDE_HASH_H__
#endif
