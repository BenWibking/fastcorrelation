#include "hash.h"

int main(int argc, char *argv[])
{
  /* input: number of points to test, number of cells per axis */
  const long seed = 42;
  const int nbins = 10;
  const double Lbox = 1.0;
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

  size_t npoints;
  int ngrid;
  npoints = (size_t)input_npoints;
  ngrid = (int)input_ngrid;

  /* generate random points (x,y,z) in unit cube */
  // separate arrays (or Fortran-style arrays) are necessary both for SIMD and cache efficiency
  FLOAT *x = (FLOAT*) _mm_malloc(npoints*sizeof(FLOAT),64);
  FLOAT *y = (FLOAT*) _mm_malloc(npoints*sizeof(FLOAT),64);
  FLOAT *z = (FLOAT*) _mm_malloc(npoints*sizeof(FLOAT),64);
  
  const gsl_rng_type * T;
  gsl_rng * r;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
  gsl_rng_set(r, seed); /* Seeding random distribution */
  
  size_t n;
  for(n=0;n<npoints;n++)
    {
      x[n] = gsl_rng_uniform(r); 
      y[n] = gsl_rng_uniform(r); 
      z[n] = gsl_rng_uniform(r); 
    }

  gsl_rng_free(r);

  /* hash into grid cells */
  GHash *grid = allocate_hash(ngrid, Lbox, npoints);
  if ((int)grid == 0) {
    printf("allocating grid failed!\n");
    exit(-1);
  }

  geometric_hash(grid, x, y, z, npoints);

  //  free(x); free(y); free(z);
  
  /* output array of counts in cells */
  int i,j,k;
  for(i=0;i<ngrid;i++) {
    for(j=0;j<ngrid;j++) {
      for(k=0;k<ngrid;k++) {
	//	printf("[%ld][%ld][%ld] = %ld\n",i,j,k,grid->counts[INDEX(i,j,k)]);
      }
    }
  }

  /* compute pair counts */
  double *bin_edges_sq = _mm_malloc((nbins+1)*sizeof(double),64);
  double *bin_edges = _mm_malloc((nbins+1)*sizeof(double),64);
  long int *pcounts = _mm_malloc(nbins*sizeof(long int),64);
  long int *pcounts_naive = _mm_malloc(nbins*sizeof(long int),64);
  for(i=0;i<nbins;i++) {
    pcounts[i] = (long int) 0;
    pcounts_naive[i] = (long int) 0;
  }
  double maxr = grid->Lbox/(double)(grid->ngrid);
  double minr = 0.01;
  double dlogr = (log10(maxr)-log10(minr))/(double)nbins;
  for(i=0;i<=nbins;i++) {
    double bin_edge = pow(10.0, ((double)i)*dlogr + log10(minr));
    bin_edges[i] = bin_edge;
    bin_edges_sq[i] = SQ(bin_edge);
  }

  count_pairs(grid, pcounts, bin_edges_sq, nbins);
  if(testing)
    count_pairs_naive(x,y,z, npoints, pcounts_naive, bin_edges_sq, nbins, Lbox);

  for(i=0;i<nbins;i++) {
    double ndens = npoints/CUBE(Lbox);
    double exp_counts = (2./3.)*M_PI*(CUBE(bin_edges[i+1])-CUBE(bin_edges[i]))*ndens*npoints;
    printf("pair counts between (%lf, %lf] = %ld\n",bin_edges[i],bin_edges[i+1],pcounts[i]);
    if(testing)
      printf("(naive) pair counts between (%lf, %lf] = %ld\n",bin_edges[i],bin_edges[i+1],pcounts_naive[i]);
    printf("expected pair counts between (%lf, %lf] = %lf\n\n",bin_edges[i],bin_edges[i+1],exp_counts);
  }

  _mm_free(pcounts);
  _mm_free(bin_edges);
  _mm_free(bin_edges_sq);
  
  free_hash(grid);
  return 0;
}
