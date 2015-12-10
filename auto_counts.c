#include "hash.h"

#define SIMD_WIDTH 4

void count_pairs_naive(FLOAT * x, FLOAT * y, FLOAT * z, size_t npoints, long int * pcounts, double *  bin_edges_sq, const int nbins, const double Lbox)
{
  size_t i;
  for(i=0;i<npoints;i++)
    {
      const size_t simd_size = npoints/SIMD_WIDTH;
      size_t jj;
      for(jj=0;jj<simd_size;jj++)
	{
          double dist_sq[SIMD_WIDTH];
	  size_t k;
	  __assume_aligned(x, 32);
	  __assume_aligned(y, 32);
	  __assume_aligned(z, 32);
#pragma simd
	  for(k=0;k<SIMD_WIDTH;k++)
	    {
	      const size_t kk = k+jj*SIMD_WIDTH;
	      dist_sq[k] = SQ(PERIODIC(x[i]-x[kk])) + SQ(PERIODIC(y[i]-y[kk])) + SQ(PERIODIC(z[i]-z[kk]));
	    }

	  for(k=0;k<SIMD_WIDTH;k++) {
	    size_t kk = k+jj*SIMD_WIDTH;
	    if(!(dist_sq[k] > bin_edges_sq[nbins])) {
	      int n;
	      for(n=nbins-1; n>=0; n--) {
		if(dist_sq[k] > bin_edges_sq[n]) {
		  pcounts[n]++;
		  break;
		}
	      }
	    }
	  }
	}

      size_t k;
      for(k=((simd_size)*SIMD_WIDTH);k<npoints;k++)
	{
	  double dist_sq = SQ(PERIODIC(x[i]-x[k])) + SQ(PERIODIC(y[i]-y[k])) + SQ(PERIODIC(z[i]-z[k]));
	  if(!(dist_sq > bin_edges_sq[nbins])) {
	    int n;
	    for(n=nbins-1; n>=0; n--) {
	      if(dist_sq > bin_edges_sq[n]) {
		pcounts[n]++;
		break;
	      }
	    }
	  }
	}
    }

  for(i=0;i<nbins;i++) {
    pcounts[i] = pcounts[i]/2;
  }  
}

void count_pairs(GHash * restrict g, long int * restrict pcounts, double * restrict bin_edges_sq, int nbins)
{
  int ngrid = g->ngrid;
  double Lbox = g->Lbox;
  int ix,iy,iz;
  for(ix=0;ix<ngrid;ix++) {
    for(iy=0;iy<ngrid;iy++) {
      for(iz=0;iz<ngrid;iz++) {

	//	printf("ix = %ld; iy = %ld; iz = %ld\n",ix,iy,iz);

	/* do this for each cell */
	size_t count = g->counts[INDEX(ix,iy,iz)];
	FLOAT * restrict x = g->x[INDEX(ix,iy,iz)];
	FLOAT * restrict y = g->y[INDEX(ix,iy,iz)];
	FLOAT * restrict z = g->z[INDEX(ix,iy,iz)];
	/* // scalar version
	size_t i,j;
	for(i=0;i<count;i++)
	  {
	    for(j=0;j<count;j++)
	      {
		if(i==j)
		  continue;

		double dist_sq = SQ(PERIODIC(x[i]-x[j])) + SQ(PERIODIC(y[i]-y[j])) + SQ(PERIODIC(z[i]-z[j]));
		int n;
		if(!(dist_sq > bin_edges_sq[nbins])) {
		  for(n=nbins-1; n>=0; n--) {
		    if(dist_sq > bin_edges_sq[n]) {
		      pcounts[n]++;
		      break;
		    }
		  }
		}
	      }
	      } */

	size_t i;
	for(i=0;i<count;i++)
	  {
	    const size_t simd_size = count/SIMD_WIDTH;
	    size_t jj;
	    for(jj=0;jj<simd_size;jj++)
	      {
		double dist_sq[SIMD_WIDTH];
		size_t k;
		__assume_aligned(x, 32);
		__assume_aligned(y, 32);
		__assume_aligned(z, 32);
#pragma simd
		for(k=0;k<SIMD_WIDTH;k++)
		  {
		    const size_t kk = k+jj*SIMD_WIDTH;
		    dist_sq[k] = SQ(PERIODIC(x[i]-x[kk])) + SQ(PERIODIC(y[i]-y[kk])) + SQ(PERIODIC(z[i]-z[kk]));
		  }
		
		for(k=0;k<SIMD_WIDTH;k++) {
		  size_t kk = k+jj*SIMD_WIDTH;
		  if(!(dist_sq[k] > bin_edges_sq[nbins])) {
		    int n;
		    for(n=nbins-1; n>=0; n--) {
		      if(dist_sq[k] > bin_edges_sq[n]) {
			pcounts[n]++;
			break;
		      }
		    }
		  }
		}
	      }
	    
	    size_t k;
	    for(k=((simd_size)*SIMD_WIDTH);k<count;k++)
	      {
		double dist_sq = SQ(PERIODIC(x[i]-x[k])) + SQ(PERIODIC(y[i]-y[k])) + SQ(PERIODIC(z[i]-z[k]));
		if(!(dist_sq > bin_edges_sq[nbins])) {
		  int n;
		  for(n=nbins-1; n>=0; n--) {
		    if(dist_sq > bin_edges_sq[n]) {
		      pcounts[n]++;
		      break;
		    }
		  }
		}
	      }
	  }
	
	int iix,iiy,iiz;
	for(iix=-1;iix<=1;iix++) {
	  for(iiy=-1;iiy<=1;iiy++) {
	    for(iiz=-1;iiz<=1;iiz++) {
	      if(iix==0 && iiy==0 && iiz==0)
		continue;

	      int aix = (ix+iix+ngrid) % ngrid; // careful to ensure this is nonnegative!
	      int aiy = (iy+iiy+ngrid) % ngrid;
	      int aiz = (iz+iiz+ngrid) % ngrid;

	      //	      printf("aix = %d; aiy = %d; aiz = %d\n",aix,aiy,aiz);
	      //	      printf("aix = %d; aiy = %d; aiz = %d\n",iix,iiy,iiz);
	      
	      /* now count pairs with adjacent cells */
	      size_t adj_count = g->counts[INDEX(aix,aiy,aiz)];
	      FLOAT * restrict adj_x = g->x[INDEX(aix,aiy,aiz)];
	      FLOAT * restrict adj_y = g->y[INDEX(aix,aiy,aiz)];
	      FLOAT * restrict adj_z = g->z[INDEX(aix,aiy,aiz)];

	      /* // scalar version
	      size_t i,j;
	      for(i=0;i<count;i++) {
		for(j=0;j<adj_count;j++) {
		  double dist_sq = SQ(PERIODIC(x[i]-adj_x[j])) + SQ(PERIODIC(y[i]-adj_y[j])) + SQ(PERIODIC(z[i]-adj_z[j]));
		  int n;
		  if(!(dist_sq > bin_edges_sq[nbins])) {
		    for(n=nbins-1; n>=0; n--) {
		      if(dist_sq > bin_edges_sq[n]) {
			pcounts[n]++;
			break;
		      }
		    }
		  }
		}
		} */
	      
	      /* SIMD version */
	      size_t i;
	      for(i=0;i<count;i++)
		{
		  const size_t simd_size = adj_count/SIMD_WIDTH;
		  size_t jj;
		  for(jj=0;jj<simd_size;jj++)
		    {
		      double dist_sq[SIMD_WIDTH];
		      size_t k;
		      __assume_aligned(adj_x, 32);
		      __assume_aligned(adj_y, 32);
		      __assume_aligned(adj_z, 32);
#pragma simd
		      for(k=0;k<SIMD_WIDTH;k++)
			{
			  const size_t kk = k+jj*SIMD_WIDTH;
			  dist_sq[k] = SQ(PERIODIC(x[i]-adj_x[kk])) + SQ(PERIODIC(y[i]-adj_y[kk])) + SQ(PERIODIC(z[i]-adj_z[kk]));
			}
		      
		      for(k=0;k<SIMD_WIDTH;k++) {
			size_t kk = k+jj*SIMD_WIDTH;
			if(!(dist_sq[k] > bin_edges_sq[nbins])) {
			  int n;
			  for(n=nbins-1; n>=0; n--) {
			    if(dist_sq[k] > bin_edges_sq[n]) {
			      pcounts[n]++;
			      break;
			    }
			  }
			}
		      }
		    }
		  
		  size_t k;
		  for(k=((simd_size)*SIMD_WIDTH);k<adj_count;k++)
		    {
		      double dist_sq = SQ(PERIODIC(x[i]-adj_x[k])) + SQ(PERIODIC(y[i]-adj_y[k])) + SQ(PERIODIC(z[i]-adj_z[k]));
		      if(!(dist_sq > bin_edges_sq[nbins])) {
			int n;
			for(n=nbins-1; n>=0; n--) {
			  if(dist_sq > bin_edges_sq[n]) {
			    pcounts[n]++;
			    break;
			  }
			}
		      }
		    }
		}

	    }
	  }
	}
      }
    }
  }

  size_t i;
  for(i=0;i<nbins;i++) {
    pcounts[i] = pcounts[i]/2;
  }
}
