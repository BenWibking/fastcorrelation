#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "hdf5.h"
#include "hdf5_hl.h"

typedef struct
{
  float x;
  float y;
  float z;
} particle;

void* read_halo_hdf5(char infile[],char dataset_name[],size_t *len);

