#include "linalg.h"
#include <stdlib.h>

struct _matrix
{
    float *data;
    unsigned r;
    unsigned c;
};

struct _csr_matrix
{
    float *data;
    unsigned *col_ind;
    unsigned *row_ptr;
    unsigned r;
    unsigned c;
};

struct _vector 
{
    float *data;
    unsigned dim;
};

static inline void reset_ptr(void **ptr)
{
    free(*ptr);
    *ptr = NULL;
}

void delete_matrix(matrix *pm)
{
    reset_ptr((void **) &((*pm)->data));
    reset_ptr((void **) pm);
}

void delete_csr_matrix(csr_matrix *pcm)
{
    reset_ptr((void **) &((*pcm)->data));
    reset_ptr((void **) &((*pcm)->col_ind));
    reset_ptr((void **) &((*pcm)->row_ptr));
    reset_ptr((void **) pcm);
}

void delete_vector(vector *pv)
{
    reset_ptr((void **) &((*pv)->data));
    reset_ptr((void **) pv);
}