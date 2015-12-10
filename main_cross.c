#include "hash.h"

int main(int argc, char *argv[])
{
  /* input: number of points to test, number of cells per axis */
  const long seed = 42;
  const int nbins = 20;
  double Lbox = 1.0;
  int input_npoints, input_ngrid;
  float input_boxsize;

  /* check inputs */
  if(argc != 4) {
    printf("./hash ngrid npoints_on_side box_size\n");
    exit(-1);
  }

  input_ngrid = atoi(argv[1]);
  input_npoints = atoi(argv[2]);
  input_boxsize = atof(argv[3]);

  if(input_ngrid <= 0) {
    printf("ngrid must be positive!\n");
    exit(-1);
  }
  if(input_npoints <= 0) {
    printf("npoints must be positive\n");
    exit(-1);
  }
  if(input_boxsize <= 0) {
    printf("boxsize must be positive!\n");
    exit(-1);
  }

  size_t npoints1, npoints2;
  int ngrid;
  npoints1 = CUBE((size_t)input_npoints);
  npoints2 = CUBE((size_t)input_npoints);
  ngrid = (int)input_ngrid;
  Lbox = (double)input_boxsize;

  printf("computing with %ld random points in a (%lf)^3 periodic box...\n",npoints1,Lbox);

  /* generate random points (x,y,z) in unit cube */
  // separate arrays (or Fortran-style arrays) are necessary both for SIMD and cache efficiency
  FLOAT *x1 = (FLOAT*) my_malloc(npoints1*sizeof(FLOAT));
  FLOAT *y1 = (FLOAT*) my_malloc(npoints1*sizeof(FLOAT));
  FLOAT *z1 = (FLOAT*) my_malloc(npoints1*sizeof(FLOAT));

  FLOAT *x2 = (FLOAT*) my_malloc(npoints2*sizeof(FLOAT));
  FLOAT *y2 = (FLOAT*) my_malloc(npoints2*sizeof(FLOAT));
  FLOAT *z2 = (FLOAT*) my_malloc(npoints2*sizeof(FLOAT));
  
  const gsl_rng_type * T;
  gsl_rng * r;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
  gsl_rng_set(r, seed); /* Seeding random distribution */
  
  size_t n;
  for(n=0;n<npoints1;n++)
    {
      x1[n] = gsl_rng_uniform(r)*Lbox; 
      y1[n] = gsl_rng_uniform(r)*Lbox; 
      z1[n] = gsl_rng_uniform(r)*Lbox; 
    }
  for(n=0;n<npoints2;n++)
    {
      x2[n] = gsl_rng_uniform(r)*Lbox; 
      y2[n] = gsl_rng_uniform(r)*Lbox; 
      z2[n] = gsl_rng_uniform(r)*Lbox; 
    }

  gsl_rng_free(r);

  /* hash into grid cells */
  GHash *grid1 = allocate_hash(ngrid, Lbox, npoints1, x1, y1, z1);
  GHash *grid2 = allocate_hash(ngrid, Lbox, npoints2, x2, y2, z2);
  if ((int)grid1 == 0) {
    printf("allocating grid1 failed!\n");
    exit(-1);
  }
  if ((int)grid2 == 0) {
    printf("allocating grid2 failed!\n");
    exit(-1);
  }

  geometric_hash(grid1, x1, y1, z1, npoints1);
  geometric_hash(grid2, x2, y2, z2, npoints2);
  

  /* compute pair counts */
  double *bin_edges_sq = my_malloc((nbins+1)*sizeof(double));
  double *bin_edges = my_malloc((nbins+1)*sizeof(double));
  long int *pcounts = my_malloc(nbins*sizeof(long int));
  long int *pcounts_naive = my_malloc(nbins*sizeof(long int));
  int i;
  for(i=0;i<nbins;i++) {
    pcounts[i] = (long int) 0;
    pcounts_naive[i] = (long int) 0;
  }
  double maxr = grid1->Lbox/(double)(grid1->ngrid);
  double minr = 0.1;
  double dlogr = (log10(maxr)-log10(minr))/(double)nbins;
  for(i=0;i<=nbins;i++) {
    double bin_edge = pow(10.0, ((double)i)*dlogr + log10(minr));
    bin_edges[i] = bin_edge;
    bin_edges_sq[i] = SQ(bin_edge);
  }

  cross_count_pairs(grid1, grid2, pcounts, bin_edges_sq, nbins);
  if(testing)
    cross_count_pairs_naive(x1,y1,z1,npoints1, x2,y2,z2,npoints2, pcounts_naive, bin_edges_sq, nbins, Lbox);

  for(i=0;i<nbins;i++) {
    double ndens1 = npoints1/CUBE(Lbox);
    double exp_counts = (4./3.)*M_PI*(CUBE(bin_edges[i+1])-CUBE(bin_edges[i]))*ndens1*npoints2;
    printf("pair counts between (%lf, %lf] = %ld\n",bin_edges[i],bin_edges[i+1],pcounts[i]);
    if(testing)
      printf("(naive) pair counts between (%lf, %lf] = %ld\n",bin_edges[i],bin_edges[i+1],pcounts_naive[i]);
    printf("expected pair counts between (%lf, %lf] = %lf\n\n",bin_edges[i],bin_edges[i+1],exp_counts);
  }

  my_free(pcounts);
  my_free(bin_edges);
  my_free(bin_edges_sq);
  
  my_free(x1);
  my_free(y1);
  my_free(z1);
  my_free(x2);
  my_free(y2);
  my_free(z2);

  free_hash(grid1);
  free_hash(grid2);
  return 0;
}
