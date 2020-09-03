#include "linalg.h"
#include <stdlib.h>

struct matrix
{
    float *data;
    unsigned r;
    unsigned c;
};

struct csr_matrix
{
    float *data;
    unsigned *col_ind;
    unsigned *row_ptr;
    unsigned r;
    unsigned c;
}

struct vector 
{
    float *data;
    unsigned dim;
};

void delete(matrix_ptr *pm)
{
    
}

void delete(csr_matrix_ptr *pcm)
{
    
}

void delete(vector_ptr     *pv )
{
    
}