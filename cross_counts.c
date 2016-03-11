#include "hash.h"

#define SIMD_WIDTH 4

void count_pairs_disjoint(FLOAT * x1, FLOAT * y1, FLOAT * z1, grid_id * restrict label1, size_t npoints1, FLOAT * x2, FLOAT * y2, FLOAT * z2, grid_id * restrict label2, size_t npoints2, long int * pcounts, double *  bin_edges_sq, const int nbins, const double Lbox)
{
  size_t i;
  int Nj = pow(njack, 1.0/3.0)
  for(i=0;i<npoints1;i++)
    {
      const size_t simd_size = npoints2/SIMD_WIDTH;
      size_t jj;
      for(jj=0;jj<simd_size;jj++)
	{
          double dist_sq[SIMD_WIDTH];
	  size_t k;
#ifdef __INTEL_COMPILER
	  __assume_aligned(x1, 32);
	  __assume_aligned(y1, 32);
	  __assume_aligned(z1, 32);
	  __assume_aligned(x2, 32);
	  __assume_aligned(y2, 32);
	  __assume_aligned(z2, 32);
#endif
	  //#pragma simd
	  for(k=0;k<SIMD_WIDTH;k++)
	    {
	      const size_t kk = k+jj*SIMD_WIDTH;
	      dist_sq[k] = SQ(PERIODIC(x1[i]-x2[kk])) + SQ(PERIODIC(y1[i]-y2[kk])) + SQ(PERIODIC(z1[i]-z2[kk]));
	    }

	  for(k=0;k<SIMD_WIDTH;k++) {
	    if(!(dist_sq[k] > bin_edges_sq[nbins])) {
	      int n;
	      for(n=nbins-1; n>=0; n--) {
		if(dist_sq[k] > bin_edges_sq[n]) {
		  pcounts[n]++;
		  int a = label[i].x;
		  int b = label[i].y;
		  int c = label[i].z;
		  int nsample = c + Nj*(b + Nj*a);
		  for(int p=0;p<njack;p++){
		  	if(p==nsample){
		  		pcounts_jackknife[nsample*nbins + n]++;
		  	} /* Bootstrap */
		  	/*if(p=!nsample){
		  		pcounts_jackknife[nsample*nbins + n]++;
		  	}*/ /* Jackknife */
		  }
		  break;
		}
	      }
	    }
	  }
	}

      size_t k;
      for(k=((simd_size)*SIMD_WIDTH);k<npoints2;k++)
	{
	  double dist_sq = SQ(PERIODIC(x1[i]-x2[k])) + SQ(PERIODIC(y1[i]-y2[k])) + SQ(PERIODIC(z1[i]-z2[k]));
	  if(!(dist_sq > bin_edges_sq[nbins])) {
	    int n;
	    for(n=nbins-1; n>=0; n--) {
	      if(dist_sq > bin_edges_sq[n]) {
		pcounts[n]++;
		int a = label[i].x;
		int b = label[i].y;
		int c = label[i].z;
		int nsample = c + Nj*(b + Nj*a);
		for(int p=0;p<njack;p++){
			if(p==nsample){
				pcounts_jackknife[nsample*nbins + n]++;
		  	} /* Bootstrap */
		  	/*if(p=!nsample){
		  		pcounts_jackknife[nsample*nbins + n]++;
		  	}*/ /* Jackknife */
		  }
		break;
	      }
	    }
	  }
	}
    }
}

void cross_count_pairs_naive(FLOAT * x1, FLOAT * y1, FLOAT * z1, size_t npoints1, FLOAT * x2, FLOAT * y2, FLOAT * z2, size_t npoints2, long int * pcounts, double *  bin_edges_sq, const int nbins, const double Lbox)
{
  count_pairs_disjoint(x1,y1,z1,label1,npoints1,x2,y2,z2,label2,npoints2,\
		       pcounts,bin_edges_sq,nbins,Lbox);
}

void cross_count_pairs(GHash * restrict g1, GHash * restrict g2, long int * restrict pcounts, double * restrict bin_edges_sq, int nbins)
{
  /* check that g1 and g2 have the same ngrid and Lbox */
  if(g1->ngrid != g2->ngrid) {
    printf("grids do not align!\n");
    exit(-1);
  }
  if(g1->Lbox != g2->Lbox) {
    printf("box geometries do not align!\n");
    exit(-1);
  }

  int ngrid = g1->ngrid;
  double Lbox = g1->Lbox;

  int ix,iy,iz;
  for(ix=0;ix<ngrid;ix++) {
    for(iy=0;iy<ngrid;iy++) {
      for(iz=0;iz<ngrid;iz++) {

	/* do this for each cell */
	size_t count1 = g1->counts[INDEX(ix,iy,iz)];
	size_t count2 = g2->counts[INDEX(ix,iy,iz)];
	FLOAT * restrict x1 = g1->x[INDEX(ix,iy,iz)];
	FLOAT * restrict y1 = g1->y[INDEX(ix,iy,iz)];
	FLOAT * restrict z1 = g1->z[INDEX(ix,iy,iz)];
	grid_id * restrict label1 = g1->sample_excluded_from[INDEX(ix, iy, iz)];
	
	FLOAT * restrict x2 = g2->x[INDEX(ix,iy,iz)];
	FLOAT * restrict y2 = g2->y[INDEX(ix,iy,iz)];
	FLOAT * restrict z2 = g2->z[INDEX(ix,iy,iz)];
	grid_id * restrict label2 = g2->sample_excluded_from[INDEX(ix, iy, iz)];
	
	count_pairs_disjoint(x1,y1,z1,label1,count1,x2,y2,z2,label2,count2,\
		       pcounts,bin_edges_sq,nbins,Lbox);
	
	int iix,iiy,iiz;
	for(iix=-1;iix<=1;iix++) {
	  for(iiy=-1;iiy<=1;iiy++) {
	    for(iiz=-1;iiz<=1;iiz++) {
	      if(iix==0 && iiy==0 && iiz==0)
		continue;

	      int aix = (ix+iix+ngrid) % ngrid; // careful to ensure this is nonnegative!
	      int aiy = (iy+iiy+ngrid) % ngrid;
	      int aiz = (iz+iiz+ngrid) % ngrid;

	      /* now count pairs with adjacent cells */
	      size_t adj_count = g2->counts[INDEX(aix,aiy,aiz)];
	      FLOAT * restrict adj_x = g2->x[INDEX(aix,aiy,aiz)];
	      FLOAT * restrict adj_y = g2->y[INDEX(aix,aiy,aiz)];
	      FLOAT * restrict adj_z = g2->z[INDEX(aix,aiy,aiz)];
	      grid_id * restrict adj_label = g2->sample_excluded_from[INDEX(ix, iy, iz)];

	      count_pairs_disjoint(x1,y1,z1,label1,count1,adj_x,adj_y,adj_z,adj_label,adj_count,\
				   pcounts,bin_edges_sq,nbins,Lbox);

	    }
	  }
	}
      }
    }
  }
}
