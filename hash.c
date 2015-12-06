#include "hash.h"

void* my_malloc(size_t size)
{
  return _mm_malloc(size,64);
}

void my_free(void* block)
{
  _mm_free(block);
}

void* my_realloc(void* old_block, size_t new_size, size_t old_size)
{
  /* emulate the behavior of realloc, ensuring alignment (but we always have to do a memcpy) */
  void* new_block = _mm_malloc(new_size,64);
  memcpy(new_block, old_block, old_size);
  _mm_free(old_block);
  return new_block;
}

GHash* allocate_hash(int ngrid, double Lbox, size_t npoints)
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

  size_t part_per_cell = (size_t)MAX(ceil(npoints/(ngrid*ngrid*ngrid)),1);

  /* now allocate cells */
  int i,j,k;
  for(i=0;i<ngrid;i++) {
    for(j=0;j<ngrid;j++) {
      for(k=0;k<ngrid;k++) {
	g->x[INDEX(i,j,k)] = my_malloc(part_per_cell*sizeof(FLOAT));
	g->y[INDEX(i,j,k)] = my_malloc(part_per_cell*sizeof(FLOAT));
	g->z[INDEX(i,j,k)] = my_malloc(part_per_cell*sizeof(FLOAT));
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
      grid->x[INDEX(ix,iy,iz)] = my_realloc(grid->x[INDEX(ix,iy,iz)], \
					    new_size*sizeof(FLOAT), old_size*sizeof(FLOAT));
      grid->y[INDEX(ix,iy,iz)] = my_realloc(grid->y[INDEX(ix,iy,iz)], \
					    new_size*sizeof(FLOAT), old_size*sizeof(FLOAT));
      grid->z[INDEX(ix,iy,iz)] = my_realloc(grid->z[INDEX(ix,iy,iz)],	\
					    new_size*sizeof(FLOAT), old_size*sizeof(FLOAT));
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
