#include "hash.h"

int main(int argc, char *argv[])
{
  /* input: number of points to test, number of cells per axis */
  int nbins;
  double Lbox, minr, maxr;
  int input_nbins;
  float input_boxsize, input_rmin, input_rmax, input_njackknife;
  char *filename;

  /* check inputs */
  if(argc != 7) {
    printf("./auto nbins rmin rmax box_size njackknife_samples filename\n");
    exit(-1);
  }

  input_nbins = atoi(argv[1]);
  input_rmin = atof(argv[2]);
  input_rmax = atof(argv[3]);
  input_boxsize = atof(argv[4]);
  input_njackknife = atoi(argv[5]);

  filename = malloc(sizeof(char)*(strlen(argv[6])+1));
  if(filename) { /* filename is not null */
    sprintf(filename,"%s",argv[6]);
  } else {
    printf("malloc failure! cannot allocate filename array\n");
    exit(-1);
  }

  if(input_nbins <= 0) {
    printf("ngrid must be positive!\n");
    exit(-1);
  }
  if(input_rmin <= 0.) {
    printf("rmin must be positive!\n");
    exit(-1);
  }
  if(!(input_rmax > input_rmin)) {
    printf("rmax must be greater than rmin!\n");
    exit(-1);
  }
  if(input_boxsize <= 0.) {
    printf("boxsize must be positive!\n");
    exit(-1);
  }
  if(input_njackknife < 0) {
    printf("njackknife_samples must be nonnegative!\n");
    exit(-1);
  }
  if(pow(floor(pow((float)input_njackknife, 1./3.)), 3) != input_njackknife) {
    printf("njackknife_samples must be a perfect cube!\n");
    exit(-1);
  }

  size_t npoints;
  int ngrid, njack;
  nbins = input_nbins;
  njack = input_njackknife;
  minr = input_rmin;
  maxr = input_rmax;
  Lbox = (double)input_boxsize;
  /* compute ngrid from rmax */
  ngrid = (int)floor(Lbox/maxr);

  /* read from file */
  particle *points = read_particles_hdf5(filename, "particles", &npoints);

  // separate arrays (or Fortran-style arrays) are necessary both for SIMD and cache efficiency
  FLOAT *x = (FLOAT*) my_malloc(npoints*sizeof(FLOAT));
  FLOAT *y = (FLOAT*) my_malloc(npoints*sizeof(FLOAT));
  FLOAT *z = (FLOAT*) my_malloc(npoints*sizeof(FLOAT));

  size_t n;
  for(n=0;n<npoints;n++)
    {
      x[n] = points[n].x;
      y[n] = points[n].y;
      z[n] = points[n].z;
    }

  free(points);

  /* hash into grid cells */
  GHash *grid = allocate_hash(ngrid, njack, Lbox, npoints, x, y, z);
  if ((int)grid == 0) {
    printf("allocating grid failed!\n");
    exit(-1);
  }

  fprintf(stderr,"computing geometric hash...");
  geometric_hash(grid, x, y, z, npoints);
  fprintf(stderr,"done!\n");

  /* compute pair counts assuming periodic box */
  double *bin_edges_sq = my_malloc((nbins+1)*sizeof(double));
  double *bin_edges = my_malloc((nbins+1)*sizeof(double));
  long int *pcounts = my_malloc(nbins*sizeof(long int));
  long int *pcounts_jackknife = my_malloc(njack*nbins*sizeof(long int));
  long int *pcounts_naive = my_malloc(nbins*sizeof(long int));
  long int *pcounts_jackknife_naive = my_malloc(njack*nbins*sizeof(long int));
  int i;
  for(i=0;i<nbins;i++) {
    pcounts[i] = (long int) 0;
    pcounts_naive[i] = (long int) 0;
  }
  for(int i=0;i<njack;i++) {
    for(int j=0;j<nbins;j++) {
      pcounts_jackknife[i*nbins + j] = (long int) 0;
    }
  }
  double dlogr = (log10(maxr)-log10(minr))/(double)nbins;
  for(i=0;i<=nbins;i++) {
    double bin_edge = pow(10.0, ((double)i)*dlogr + log10(minr));
    bin_edges[i] = bin_edge;
    bin_edges_sq[i] = SQ(bin_edge);
  }

  fprintf(stderr,"computing pair counts...");
  count_pairs(grid, pcounts, pcounts_jackknife, bin_edges_sq, nbins);
  fprintf(stderr,"done!\n");

#ifdef TEST_ALL_PAIRS
  count_pairs_naive(x,y,z, npoints, pcounts_naive, pcounts_jackknife, bin_edges_sq, nbins, njack, Lbox);
#endif

  /* output pair counts */
  printf("min_bin\tmax_bin\tbin_counts\tnatural_estimator\n");

  for(i=0;i<nbins;i++) {
    double ndens = npoints/CUBE(Lbox);
    double exp_counts = (2./3.)*M_PI*(CUBE(bin_edges[i+1])-CUBE(bin_edges[i]))*ndens*npoints;
#ifdef JACKKNIFE_SUBSAMPLES
    double exp_counts_jackknife = exp_counts*(double)(((double)njack-1.0)/(double)njack);
#else /* bootstrap subsamples */
    double exp_counts_jackknife = exp_counts*(double)(1.0/(double)njack);
#endif
    printf("%lf\t%lf\t%ld\t%lf", bin_edges[i],bin_edges[i+1],pcounts[i],(double)pcounts[i]/exp_counts - 1.0);
    for(int j=0;j<njack;j++) {
      printf("\t%lf",(double)pcounts_jackknife[j*nbins + i]/exp_counts_jackknife - 1.0);
    }
    printf("\n");
#ifdef TEST_ALL_PAIRS
    printf("(naive) pair counts = %ld\n",pcounts_naive[i]);
#endif
  }

  /* free memory */
  my_free(pcounts);
  my_free(pcounts_jackknife);
  my_free(bin_edges);
  my_free(bin_edges_sq);

  my_free(x);
  my_free(y);
  my_free(z);
  
  free_hash(grid);
  return 0;
}
