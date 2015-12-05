#include "hash.h"

void count_pairs_naive(double *x, double *y, double *z, size_t npoints, long int * pcounts, double * bin_edges_sq, int nbins, double Lbox)
{
  size_t i,j;
#pragma omp simd collapse(2)
  for(i=0;i<npoints;i++)
    {
      for(j=0;j<npoints;j++)
	{
	  if(i==j)
	    continue;

	  double dist_sq = SQ(PERIODIC(x[i]-x[j])) + SQ(PERIODIC(y[i]-y[j])) + SQ(PERIODIC(z[i]-z[j]));
	  int n;
	  if(!(dist_sq > bin_edges_sq[nbins])) {
	    for(n=nbins-1; n>=0; n--) {
	      if((dist_sq > bin_edges_sq[n]) && (dist_sq < bin_edges_sq[n+1])) {
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

void count_pairs(GHash * g, long int * pcounts, double * bin_edges_sq, int nbins)
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
	FLOAT * x = g->x[INDEX(ix,iy,iz)];
	FLOAT * y = g->y[INDEX(ix,iy,iz)];
	FLOAT * z = g->z[INDEX(ix,iy,iz)];
	size_t i,j;
#pragma omp simd collapse(2)
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
		    if((dist_sq > bin_edges_sq[n]) && (dist_sq < bin_edges_sq[n+1])) {
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
	      FLOAT * adj_x = g->x[INDEX(aix,aiy,aiz)];
	      FLOAT * adj_y = g->y[INDEX(aix,aiy,aiz)];
	      FLOAT * adj_z = g->z[INDEX(aix,aiy,aiz)];

	      size_t i,j;
#pragma omp simd collapse(2)
	      for(i=0;i<count;i++) {
		for(j=0;j<adj_count;j++) {
		  double dist_sq = SQ(PERIODIC(x[i]-adj_x[j])) + SQ(PERIODIC(y[i]-adj_y[j])) + SQ(PERIODIC(z[i]-adj_z[j]));
		  int n;
		  if(!(dist_sq > bin_edges_sq[nbins])) {
		    for(n=nbins-1; n>=0; n--) {
		      if((dist_sq > bin_edges_sq[n]) && (dist_sq < bin_edges_sq[n+1])) {
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
