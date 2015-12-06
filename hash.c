#include "hash.h"

GHash* allocate_hash(int ngrid, double Lbox, size_t npoints)
{
  if((ngrid <= 0) || (npoints <= 0)) {
    return (void*)0;
  }

  GHash * g    = _mm_malloc(sizeof(GHash),64);
  g->ngrid = ngrid;
  g->Lbox = Lbox;

  g->counts    = _mm_malloc(ngrid*ngrid*ngrid*sizeof(size_t),64);
  g->allocated = _mm_malloc(ngrid*ngrid*ngrid*sizeof(size_t),64);
  g->x         = _mm_malloc(ngrid*ngrid*ngrid*sizeof(FLOAT*),64);
  g->y         = _mm_malloc(ngrid*ngrid*ngrid*sizeof(FLOAT*),64);
  g->z         = _mm_malloc(ngrid*ngrid*ngrid*sizeof(FLOAT*),64);

  size_t part_per_cell = (size_t)MAX(ceil(npoints/(ngrid*ngrid*ngrid)),1);

  /* now allocate cells */
  int i,j,k;
  for(i=0;i<ngrid;i++) {
    for(j=0;j<ngrid;j++) {
      for(k=0;k<ngrid;k++) {
	g->x[INDEX(i,j,k)] = malloc(part_per_cell*sizeof(FLOAT));
	g->y[INDEX(i,j,k)] = malloc(part_per_cell*sizeof(FLOAT));
	g->z[INDEX(i,j,k)] = malloc(part_per_cell*sizeof(FLOAT));
	g->allocated[INDEX(i,j,k)] = part_per_cell;
	g->counts[INDEX(i,j,k)] = 0;
	//	printf("allocating [%ld][%ld][%ld] = %ld\n",i,j,k,g->allocated[INDEX(i,j,k)]);
      }
    }
  }
  
  return g;
}

void free_hash(GHash * g)
{
  _mm_free(g->counts);
  _mm_free(g->allocated);
  int i,j,k;
  int ngrid = g->ngrid;
  for(i=0;i<ngrid;i++) {
    for(j=0;j<ngrid;j++) {
      for(k=0;k<ngrid;k++) {
	free(g->x[INDEX(i,j,k)]);
	free(g->y[INDEX(i,j,k)]);
	free(g->z[INDEX(i,j,k)]);
      }
    }
  }
  _mm_free(g->x);
  _mm_free(g->y);
  _mm_free(g->z);
  _mm_free(g);
}

void insert_particle(GHash * grid, FLOAT x, FLOAT y, FLOAT z)
{
  /* compute ix,iy,iz coordinates */
  int ix = floor(x/grid->Lbox*((double)grid->ngrid));
  int iy = floor(y/grid->Lbox*((double)grid->ngrid));
  int iz = floor(z/grid->Lbox*((double)grid->ngrid));

  int ngrid = grid->ngrid;

  /* see if there is space in the cell */
  if(grid->counts[INDEX(ix,iy,iz)] < grid->allocated[INDEX(ix,iy,iz)])
    {
      size_t idx = grid->counts[INDEX(ix,iy,iz)];
      grid->x[INDEX(ix,iy,iz)][idx] = x;
      grid->y[INDEX(ix,iy,iz)][idx] = y;
      grid->z[INDEX(ix,iy,iz)][idx] = z;
      (grid->counts[INDEX(ix,iy,iz)])++;
    }
  else
    {
      /* reallocate cell with larger size */
      size_t old_size = grid->allocated[INDEX(ix,iy,iz)];
      size_t new_size = MAX(old_size*REALLOC_FAC, old_size+1);
      //      printf("realloc'ing for [%d][%d][%d]!\n",ix,iy,iz);
      grid->x[INDEX(ix,iy,iz)] = realloc(grid->x[INDEX(ix,iy,iz)], new_size*sizeof(FLOAT));
      grid->y[INDEX(ix,iy,iz)] = realloc(grid->y[INDEX(ix,iy,iz)], new_size*sizeof(FLOAT));
      grid->z[INDEX(ix,iy,iz)] = realloc(grid->z[INDEX(ix,iy,iz)], new_size*sizeof(FLOAT));
      grid->allocated[INDEX(ix,iy,iz)] = new_size;
      
      size_t idx = grid->counts[INDEX(ix,iy,iz)];
      if(new_size > idx) {
	grid->x[INDEX(ix,iy,iz)][idx] = x;
	grid->y[INDEX(ix,iy,iz)][idx] = y;
	grid->z[INDEX(ix,iy,iz)][idx] = z;
	(grid->counts[INDEX(ix,iy,iz)])++;
      } else {
	printf("realloc failed for [%d][%d][%d]!\n",ix,iy,iz);
	printf("old_size = %ld\n",old_size);
	printf("new_size = %ld\n",new_size);
	printf("idx = %ld\n",idx);
	exit(1);
      }
    }
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
