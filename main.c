#include "hash.h"

int main(int argc, char *argv[])
{
  /* input: number of points to test, number of cells per axis */
  const long seed = 42;
  int input_npoints, input_ngrid;

  /* check inputs */
  if(argc != 3) {
    printf("./hash ngrid npoints\n");
    exit(-1);
  }

  input_ngrid = atoi(argv[1]);
  input_npoints = atoi(argv[2]);

  if(input_ngrid <= 0) {
    printf("ngrid must be positive!\n");
    exit(-1);
  }
  if(input_npoints <= 0) {
    printf("npoints must be positive\n");
    exit(-1);
  }

  size_t npoints, ngrid;
  npoints = (size_t)input_npoints;
  ngrid = (size_t)input_ngrid;

  /* generate random points (x,y,z) in unit cube */
  // separate arrays (or Fortran-style arrays) are necessary both for SIMD and cache efficiency
  FLOAT *x = (FLOAT*) malloc(npoints*sizeof(FLOAT));
  FLOAT *y = (FLOAT*) malloc(npoints*sizeof(FLOAT));
  FLOAT *z = (FLOAT*) malloc(npoints*sizeof(FLOAT));
  
  const gsl_rng_type * T;
  gsl_rng * r;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
  gsl_rng_set(r, seed); /* Seeding random distribution */
  
  size_t i;
  for(i=0;i<npoints;i++)
    {
      x[i] = gsl_rng_uniform(r); 
      y[i] = gsl_rng_uniform(r); 
      z[i] = gsl_rng_uniform(r); 
    }

  gsl_rng_free(r);

  /* hash into grid cells */
  GHash *grid = allocate_hash(ngrid, 1.0, npoints);
  if ((int)grid == 0) {
    printf("allocating grid failed!\n");
    exit(-1);
  }

  geometric_hash(grid, x, y, z, npoints);
  
  /* output array of counts in cells */
  size_t j,k;
  for(i=0;i<ngrid;i++) {
    for(j=0;j<ngrid;j++) {
      for(k=0;k<ngrid;k++) {
	printf("[%ld][%ld][%ld] = %ld\n",i,j,k,grid->counts[INDEX(i,j,k)]);
      }
    }
  }

  free_hash(grid);
  return 0;
}
