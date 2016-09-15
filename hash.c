#include "hash.h"

void* my_malloc(size_t size)
{
#ifdef __INTEL_COMPILER
  void* pointer = _mm_malloc(size,32);
#else
  void* pointer = malloc(size);
#endif
  if(pointer) {
    return pointer;
  } else {
    printf("malloc failure for size: %zu\n", size);
    exit(-1);
    return 0;
  }
}

void my_free(void* block)
{
#ifdef __INTEL_COMPILER
  _mm_free(block);
#else
  free(block);
#endif
}

void* my_realloc(void* old_block, size_t new_size, size_t old_size)
{
  /* emulate the behavior of realloc, ensuring alignment (but we always have to do a memcpy) */
#ifdef __INTEL_COMPILER
  void* new_block = _mm_malloc(new_size,64);
  memcpy(new_block, old_block, old_size);
  _mm_free(old_block);
#else
  void* new_block = realloc(old_block, new_size);
#endif
  if(new_block) {
    return new_block;
  } else {
    printf("realloc failure!\n");
    exit(-1);
    return 0;
  }
}

MHash* allocate_1d_hash(int nbins, double * bin_edges, size_t npoints, halo_metadata * h)
{
  if ((nbins <= 0) || (npoints <= 0)) {
    return (void*)0;
  }

  MHash * m = my_malloc(sizeof(MHash));
  m->nbins = nbins;
  m->bin_edges = bin_edges;
  m->counts = my_malloc(nbins*sizeof(size_t));
  m->allocated = my_malloc(nbins*sizeof(size_t));
  m->h = my_malloc(nbins*sizeof(halo_metadata*));

  int i;
  for(i=0;i<nbins;i++)
    {
      m->counts[i] = 0;
      m->allocated[i] = 0;
    }

  /* compute needed allocation sizes for each bin */
  size_t n;
  for(n=0;n<npoints;n++) {
    double s = h[n].mass;
    /* linear search */
    size_t j;
    for(j=nbins-1;j>0;j--) {
      if((s < bin_edges[j]) && (s >= bin_edges[j-1])) {
	/* add space in bin[j-1] */
	m->allocated[j-1]++;
      }
      break;
    }
  }

  /* allocate each bin */
  size_t j;
  for(j=0;j<nbins;j++) {
    size_t count = m->allocated[j];
    m->h[j] = my_malloc(count*sizeof(halo_metadata));
  }

  /* add elements to each bin */
  for(n=0;n<npoints;n++) {
    double s = h[n].mass;
    /* linear search */
    size_t j;
    for(j=nbins-1;j>0;j--) {
      if((s < bin_edges[j]) && (s >= bin_edges[j-1])) {
	size_t count = m->counts[j-1];
	m->h[j-1][count] = h[n]; /* or memcpy? */
	(m->counts[j-1])++;
      }
      break;
    }
  }
  
  return m;
}

int compare_halo_metadata_by_mass(const void* a, const void* b)
{
  halo_metadata * haloA = (halo_metadata*)a;
  halo_metadata * haloB = (halo_metadata*)b;

  if( haloA->mass < haloB->mass ) return -1;
  if( haloA->mass < haloB->mass ) return 1;
  return 0;
}

int compare_halo_metadata_by_id(const void* a, const void* b)
{
  halo_metadata * haloA = (halo_metadata*)a;
  halo_metadata * haloB = (halo_metadata*)b;

  if( haloA->id < haloB->id ) return -1;
  if( haloA->id < haloB->id ) return 1;
  return 0;
}

void sort_1d_hash(MHash * m)
{
  /* sort each bin within the 1d hash */
  size_t nbins = m->nbins;

  size_t j;
  for(j=0;j<nbins;j++) {
    /* sort halos within this bin */
    size_t array_size = m->counts[j];
    halo_metadata * array_to_sort = m->h[j];
    
    /* use standard C qsort with custom sorter function */
    qsort(array_to_sort, array_size, sizeof(halo_metadata), compare_halo_metadata_by_mass);

    /* walk through array_to_sort in order, adding percentiles */
    size_t i;
    for(i=0;i<array_size;i++) {
      double percentile = ((double)(i+1))/((double)(array_size));
      array_to_sort[i].percentile = percentile;
    }
  }
}

void linearize_1d_hash(MHash * m, size_t len, halo_metadata * linear_halos)
{
  size_t nbins = m->nbins;
  size_t j;
  size_t q = 0;
  for(j=0;j<nbins;j++) {
    size_t array_size = m->counts[j];
    halo_metadata * this_halos = m->h[j];
    size_t i;
    for(i=0;i<array_size;i++) {
      linear_halos[q] = this_halos[i];
      q++;
    }
  }
}

GHash* allocate_hash(int ngrid, int njack, double Lbox, size_t npoints, FLOAT * x, FLOAT * y, FLOAT * z)
{
  if((ngrid <= 0) || (npoints <= 0) || (njack < 0)) {
    return (void*)0;
  }

  GHash * g    = my_malloc(sizeof(GHash));
  g->ngrid = ngrid;
  g->njack = (int)pow((float)njack, 1./3.);
  g->Lbox = Lbox;

  g->counts    = my_malloc(ngrid*ngrid*ngrid*sizeof(size_t));
  g->allocated = my_malloc(ngrid*ngrid*ngrid*sizeof(size_t));
  g->x         = my_malloc(ngrid*ngrid*ngrid*sizeof(FLOAT*));
  g->y         = my_malloc(ngrid*ngrid*ngrid*sizeof(FLOAT*));
  g->z         = my_malloc(ngrid*ngrid*ngrid*sizeof(FLOAT*));
  g->sample_excluded_from = my_malloc(ngrid*ngrid*ngrid*sizeof(grid_id));

  int i,j,k;
  for(i=0;i<ngrid;i++) {
    for(j=0;j<ngrid;j++) {
      for(k=0;k<ngrid;k++) {
	g->counts[INDEX(i,j,k)] = 0;
      }
    }
  }

  /* loop over particles, determine bin counts */
  size_t n;
  for(n=0;n<npoints;n++) {
    /* determine bin */
    int ix = (int)floor(x[n]/g->Lbox*((double)g->ngrid)) % g->ngrid;
    int iy = (int)floor(y[n]/g->Lbox*((double)g->ngrid)) % g->ngrid;
    int iz = (int)floor(z[n]/g->Lbox*((double)g->ngrid)) % g->ngrid;
    if((ix>=0)&&(iy>=0)&&(iz>=0)) {
      /* increment bin counter */
      g->counts[INDEX(ix,iy,iz)]++;
    } else {
      printf("WARNING: negative spatial coordinate [%ld]: %lf %lf %lf\n",\
	     n,x[n],y[n],z[n]);
      printf("Skipping!\n");
    }
  }

  /* now allocate cells */
  for(i=0;i<ngrid;i++) {
    for(j=0;j<ngrid;j++) {
      for(k=0;k<ngrid;k++) {
	size_t bin_count = g->counts[INDEX(i,j,k)];
	g->x[INDEX(i,j,k)] = my_malloc(bin_count*sizeof(FLOAT));
	g->y[INDEX(i,j,k)] = my_malloc(bin_count*sizeof(FLOAT));
	g->z[INDEX(i,j,k)] = my_malloc(bin_count*sizeof(FLOAT));
	g->sample_excluded_from[INDEX(i,j,k)] = my_malloc(bin_count*sizeof(grid_id));
	g->allocated[INDEX(i,j,k)] = bin_count;
	g->counts[INDEX(i,j,k)] = 0;
	for(n=0; n<bin_count; n++) {
	  g->sample_excluded_from[INDEX(i,j,k)][n].x = -1;
	  g->sample_excluded_from[INDEX(i,j,k)][n].y = -1;
	  g->sample_excluded_from[INDEX(i,j,k)][n].z = -1;
	}
	//	printf("allocating [%d][%d][%d] = %ld\n",i,j,k,g->allocated[INDEX(i,j,k)]);
      }
    }
  }
  
  return g;
}

GHash* allocate_hash_with_id(int ngrid, double Lbox, size_t npoints, FLOAT * x, FLOAT * y, FLOAT * z, uint64_t * id)
{
  if((ngrid <= 0) || (npoints <= 0)) {
    return (void*)0;
  }

  GHash * g    = my_malloc(sizeof(GHash));
  g->ngrid = ngrid;
  g->njack = 1;
  g->Lbox = Lbox;

  g->counts    = my_malloc(ngrid*ngrid*ngrid*sizeof(size_t));
  g->allocated = my_malloc(ngrid*ngrid*ngrid*sizeof(size_t));
  g->x         = my_malloc(ngrid*ngrid*ngrid*sizeof(FLOAT*));
  g->y         = my_malloc(ngrid*ngrid*ngrid*sizeof(FLOAT*));
  g->z         = my_malloc(ngrid*ngrid*ngrid*sizeof(FLOAT*));
  g->mass      = my_malloc(ngrid*ngrid*ngrid*sizeof(FLOAT*));
  g->id        = my_malloc(ngrid*ngrid*ngrid*sizeof(uint64_t));

  int i,j,k;
  for(i=0;i<ngrid;i++) {
    for(j=0;j<ngrid;j++) {
      for(k=0;k<ngrid;k++) {
	g->counts[INDEX(i,j,k)] = 0;
      }
    }
  }

  /* loop over particles, determine bin counts */
  size_t n;
  for(n=0;n<npoints;n++) {
    /* determine bin */
    int ix = (int)floor(x[n]/g->Lbox*((double)g->ngrid)) % g->ngrid;
    int iy = (int)floor(y[n]/g->Lbox*((double)g->ngrid)) % g->ngrid;
    int iz = (int)floor(z[n]/g->Lbox*((double)g->ngrid)) % g->ngrid;
    if((ix>=0)&&(iy>=0)&&(iz>=0)) {
      /* increment bin counter */
      g->counts[INDEX(ix,iy,iz)]++;
    } else {
      printf("WARNING: negative spatial coordinate [%ld]: %lf %lf %lf\n",\
	     n,x[n],y[n],z[n]);
      printf("Skipping!\n");
    }
  }

  /* now allocate cells */
  for(i=0;i<ngrid;i++) {
    for(j=0;j<ngrid;j++) {
      for(k=0;k<ngrid;k++) {
	size_t bin_count = g->counts[INDEX(i,j,k)];
	g->x[INDEX(i,j,k)] = my_malloc(bin_count*sizeof(FLOAT));
	g->y[INDEX(i,j,k)] = my_malloc(bin_count*sizeof(FLOAT));
	g->z[INDEX(i,j,k)] = my_malloc(bin_count*sizeof(FLOAT));
	g->mass[INDEX(i,j,k)] = my_malloc(bin_count*sizeof(FLOAT));
	g->id[INDEX(i,j,k)] = my_malloc(bin_count*sizeof(uint64_t));
	g->allocated[INDEX(i,j,k)] = bin_count;
	g->counts[INDEX(i,j,k)] = 0;
	//	printf("allocating [%d][%d][%d] = %ld\n",i,j,k,g->allocated[INDEX(i,j,k)]);
      }
    }
  }
  
  return g;
}

void free_hash(GHash * g)
{
  my_free(g->counts);
  my_free(g->allocated);
  int i,j,k;
  int ngrid = g->ngrid;
  for(i=0;i<ngrid;i++) {
    for(j=0;j<ngrid;j++) {
      for(k=0;k<ngrid;k++) {
	my_free(g->x[INDEX(i,j,k)]);
	my_free(g->y[INDEX(i,j,k)]);
	my_free(g->z[INDEX(i,j,k)]);
	my_free(g->sample_excluded_from[INDEX(i,j,k)]);
      }
    }
  }
  my_free(g->x);
  my_free(g->y);
  my_free(g->z);
  my_free(g->sample_excluded_from);
  my_free(g);
}

void free_hash_with_id(GHash * g)
{
  my_free(g->counts);
  my_free(g->allocated);
  int i,j,k;
  int ngrid = g->ngrid;
  for(i=0;i<ngrid;i++) {
    for(j=0;j<ngrid;j++) {
      for(k=0;k<ngrid;k++) {
	my_free(g->x[INDEX(i,j,k)]);
	my_free(g->y[INDEX(i,j,k)]);
	my_free(g->z[INDEX(i,j,k)]);
	my_free(g->id[INDEX(i,j,k)]);
      }
    }
  }
  my_free(g->x);
  my_free(g->y);
  my_free(g->z);
  my_free(g->id);
  my_free(g);
}

void insert_particle(GHash * grid, FLOAT x, FLOAT y, FLOAT z, size_t i)
{
  /* compute ix,iy,iz coordinates */
  int ix = (int)floor(x/grid->Lbox*((double)grid->ngrid)) % grid->ngrid;
  int iy = (int)floor(y/grid->Lbox*((double)grid->ngrid)) % grid->ngrid;
  int iz = (int)floor(z/grid->Lbox*((double)grid->ngrid)) % grid->ngrid;

  /* compute jackknife sample that is is excluded from */
  int jx = (int)floor(x/grid->Lbox*((double)grid->njack)) % grid->njack;
  int jy = (int)floor(y/grid->Lbox*((double)grid->njack)) % grid->njack;
  int jz = (int)floor(z/grid->Lbox*((double)grid->njack)) % grid->njack;

  if((ix>=0)&&(iy>=0)&&(iz>=0)) {
    int ngrid = grid->ngrid;

    size_t idx = grid->counts[INDEX(ix,iy,iz)];
    size_t allocated = grid->allocated[INDEX(ix,iy,iz)];
    grid->x[INDEX(ix,iy,iz)][idx] = x;
    grid->y[INDEX(ix,iy,iz)][idx] = y;
    grid->z[INDEX(ix,iy,iz)][idx] = z;
    grid->sample_excluded_from[INDEX(ix,iy,iz)][idx].x = jx;
    grid->sample_excluded_from[INDEX(ix,iy,iz)][idx].y = jy;
    grid->sample_excluded_from[INDEX(ix,iy,iz)][idx].z = jz;

    (grid->counts[INDEX(ix,iy,iz)])++;
  } else {
      printf("WARNING: negative spatial coordinate: %lf %lf %lf\n",\
	     x,y,z);
      printf("Skipping!\n");    
  }
}

void insert_particle_with_id(GHash * grid, FLOAT x, FLOAT y, FLOAT z, FLOAT mass, uint64_t id, size_t i)
{
  /* compute ix,iy,iz coordinates */
  int ix = (int)floor(x/grid->Lbox*((double)grid->ngrid)) % grid->ngrid;
  int iy = (int)floor(y/grid->Lbox*((double)grid->ngrid)) % grid->ngrid;
  int iz = (int)floor(z/grid->Lbox*((double)grid->ngrid)) % grid->ngrid;

  if((ix>=0)&&(iy>=0)&&(iz>=0)) {
    int ngrid = grid->ngrid;

    size_t idx = grid->counts[INDEX(ix,iy,iz)];
    size_t allocated = grid->allocated[INDEX(ix,iy,iz)];
    grid->x[INDEX(ix,iy,iz)][idx] = x;
    grid->y[INDEX(ix,iy,iz)][idx] = y;
    grid->z[INDEX(ix,iy,iz)][idx] = z;
    grid->mass[INDEX(ix,iy,iz)][idx] = mass;
    grid->id[INDEX(ix,iy,iz)][idx] = id;

    (grid->counts[INDEX(ix,iy,iz)])++;
  } else {
      printf("WARNING: negative spatial coordinate: %lf %lf %lf\n",\
	     x,y,z);
      printf("Skipping!\n");    
  }
}

void geometric_hash(GHash * grid, FLOAT *x, FLOAT *y, FLOAT *z, size_t npoints)
{
  /* compute geometric hash of points (x,y,z) */
  size_t i;
  for(i=0;i<npoints;i++)
    {
      insert_particle(grid, x[i],y[i],z[i], i);
    }
}

void geometric_hash_with_id(GHash * grid, FLOAT *x, FLOAT *y, FLOAT *z, FLOAT *mass, \
			    uint64_t * id, size_t npoints)
{
  /* compute geometric hash of points (x,y,z) */
  size_t i;
  for(i=0;i<npoints;i++)
    {
      insert_particle_with_id(grid, x[i],y[i],z[i], mass[i], id[i], i);
    }
}
