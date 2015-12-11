#include "hash.h"

void* my_malloc(size_t size)
{
#ifdef __INTEL_COMPILER
  return _mm_malloc(size,32);
#else
  return malloc(size);
#endif
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
  return new_block;
}

GHash* allocate_hash(int ngrid, double Lbox, size_t npoints, FLOAT * x, FLOAT * y, FLOAT * z)
{
  if((ngrid <= 0) || (npoints <= 0)) {
    return (void*)0;
  }

  GHash * g    = my_malloc(sizeof(GHash));
  g->ngrid = ngrid;
  g->Lbox = Lbox;

  g->counts    = my_malloc(ngrid*ngrid*ngrid*sizeof(size_t));
  g->allocated = my_malloc(ngrid*ngrid*ngrid*sizeof(size_t));
  g->x         = my_malloc(ngrid*ngrid*ngrid*sizeof(FLOAT*));
  g->y         = my_malloc(ngrid*ngrid*ngrid*sizeof(FLOAT*));
  g->z         = my_malloc(ngrid*ngrid*ngrid*sizeof(FLOAT*));

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
    int ix = floor(x[n]/g->Lbox*((double)g->ngrid));
    int iy = floor(y[n]/g->Lbox*((double)g->ngrid));
    int iz = floor(z[n]/g->Lbox*((double)g->ngrid));
    /* increment bin counter */
    g->counts[INDEX(ix,iy,iz)]++;
  }

  /* now allocate cells */
  for(i=0;i<ngrid;i++) {
    for(j=0;j<ngrid;j++) {
      for(k=0;k<ngrid;k++) {
	size_t bin_count = g->counts[INDEX(i,j,k)];
	g->x[INDEX(i,j,k)] = my_malloc(bin_count*sizeof(FLOAT));
	g->y[INDEX(i,j,k)] = my_malloc(bin_count*sizeof(FLOAT));
	g->z[INDEX(i,j,k)] = my_malloc(bin_count*sizeof(FLOAT));
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
      }
    }
  }
  my_free(g->x);
  my_free(g->y);
  my_free(g->z);
  my_free(g);
}

void insert_particle(GHash * grid, FLOAT x, FLOAT y, FLOAT z)
{
  /* compute ix,iy,iz coordinates */
  int ix = floor(x/grid->Lbox*((double)grid->ngrid));
  int iy = floor(y/grid->Lbox*((double)grid->ngrid));
  int iz = floor(z/grid->Lbox*((double)grid->ngrid));

  int ngrid = grid->ngrid;

  size_t idx = grid->counts[INDEX(ix,iy,iz)];
  grid->x[INDEX(ix,iy,iz)][idx] = x;
  grid->y[INDEX(ix,iy,iz)][idx] = y;
  grid->z[INDEX(ix,iy,iz)][idx] = z;
  (grid->counts[INDEX(ix,iy,iz)])++;
}

void geometric_hash(GHash * grid, FLOAT *x, FLOAT *y, FLOAT *z, size_t npoints)
{
  /* compute geometric hash of points (x,y,z) */
  size_t i;
  for(i=0;i<npoints;i++)
    {
      insert_particle(grid, x[i],y[i],z[i]);
    }
}
