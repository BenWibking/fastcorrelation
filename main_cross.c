#include "hash.h"

int main(int argc, char *argv[])
{
  /* input: number of points to test, number of cells per axis */
  const int nbins = 20;
  double Lbox;
  int input_ngrid;
  float input_boxsize;
  char filenameA[1000], filenameB[1000];

  /* check inputs */
  if(argc != 5) {
    printf("./cross ngrid box_size filenameA filenameB\n");
    exit(-1);
  }

  input_ngrid = atoi(argv[1]);
  input_boxsize = atof(argv[2]);

  sprintf(filenameA,"%s",argv[3]);
  sprintf(filenameB,"%s",argv[4]);

  if(input_ngrid <= 0) {
    printf("ngrid must be positive!\n");
    exit(-1);
  }
  if(input_boxsize <= 0) {
    printf("boxsize must be positive!\n");
    exit(-1);
  }

  size_t npointsA, npointsB;
  int ngrid;
  ngrid = (int)input_ngrid;
  Lbox = (double)input_boxsize;

  /* read from file */
  particle *pointsA = read_particles_hdf5(filenameA, "particles", &npointsA);
  particle *pointsB = read_particles_hdf5(filenameB, "particles", &npointsB);

  /* generate random points (x,y,z) in unit cube */
  // separate arrays (or Fortran-style arrays) are necessary both for SIMD and cache efficiency
  FLOAT *x1 = (FLOAT*) my_malloc(npointsA*sizeof(FLOAT));
  FLOAT *y1 = (FLOAT*) my_malloc(npointsA*sizeof(FLOAT));
  FLOAT *z1 = (FLOAT*) my_malloc(npointsA*sizeof(FLOAT));

  FLOAT *x2 = (FLOAT*) my_malloc(npointsB*sizeof(FLOAT));
  FLOAT *y2 = (FLOAT*) my_malloc(npointsB*sizeof(FLOAT));
  FLOAT *z2 = (FLOAT*) my_malloc(npointsB*sizeof(FLOAT));
  
  size_t n;
  for(n=0;n<npointsA;n++)
    {
      x1[n] = pointsA[n].x;
      y1[n] = pointsA[n].y;
      z1[n] = pointsA[n].z;
    }
  for(n=0;n<npointsB;n++)
    {
      //      printf("n: %ld x: %f y: %f z: %f\n",n,pointsB[n].x,pointsB[n].y,pointsB[n].z);
      x2[n] = pointsB[n].x;
      y2[n] = pointsB[n].y;
      z2[n] = pointsB[n].z;
    }

  free(pointsA);
  free(pointsB);

  /* hash into grid cells */
  GHash *grid1 = allocate_hash(ngrid, Lbox, npointsA, x1, y1, z1);
  GHash *grid2 = allocate_hash(ngrid, Lbox, npointsB, x2, y2, z2);
  if ((int)grid1 == 0) {
    printf("allocating grid1 failed!\n");
    exit(-1);
  }
  if ((int)grid2 == 0) {
    printf("allocating grid2 failed!\n");
    exit(-1);
  }

  geometric_hash(grid1, x1, y1, z1, npointsA);
  geometric_hash(grid2, x2, y2, z2, npointsB);
  

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

#ifdef TEST_ALL_PAIRS
  cross_count_pairs_naive(x1,y1,z1,npoints1, x2,y2,z2,npoints2, pcounts_naive, bin_edges_sq, nbins, Lbox);
#endif

  /* output pair counts */
  printf("min_bin\tmax_bin\tbin_counts\tnatural_estimator\n");
  for(i=0;i<nbins;i++) {
    double ndensA = npointsA/CUBE(Lbox);
    double exp_counts = (4./3.)*M_PI*(CUBE(bin_edges[i+1])-CUBE(bin_edges[i]))*ndensA*npointsB;
    printf("%lf\t%lf\t%ld\t%lf\n",bin_edges[i],bin_edges[i+1],pcounts[i],(double)pcounts[i]/exp_counts);

#ifdef TEST_ALL_PAIRS
    printf("(naive) pair counts between (%lf, %lf] = %ld\n",bin_edges[i],bin_edges[i+1],pcounts_naive[i]);
#endif
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
