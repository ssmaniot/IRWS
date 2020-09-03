#include "linalg.h"
#include <stdio.h>
#include <stdlib.h>

/* data definition */
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

/* data type new() */
matrix new_matrix(unsigned row, unsigned col) 
{
    matrix m = (matrix) malloc(sizeof(struct _matrix));
    m->data = (float *) malloc(sizeof(float) * row * col);
    m->r = row;
    m->c = col;
    return m;
}

matrix new_matrix_(float *v, unsigned row, unsigned col)
{
    matrix m = (matrix) malloc(sizeof(struct _matrix));
    m->data = v;
    m->r = row;
    m->c = col;
    return m;
}

csr_matrix new_csr_matrix(unsigned row, unsigned col)
{
    return NULL;
}

vector new_vector(unsigned dim)
{
    vector v = (vector) malloc(sizeof(struct _vector));
    v->data = (float *) malloc(sizeof(float) * dim);
    v->dim = dim;
    return v;
}

vector new_vector_(float *d, unsigned dim)
{
    vector v = (vector) malloc(sizeof(struct _vector));
    v->data = d;
    v->dim = dim;
    return v;
}

/* helper function for free/set to NULL */
static inline void reset_ptr(void **ptr)
{
    free(*ptr);
    *ptr = NULL;
}

/* delete functions */
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

/* printer functions */
void print_matrix(matrix m)
{
    unsigned i, j;
    for (i = 0; i < m->r; ++i)
    {
        for (j = 0; j < m->c; ++j)
            printf("%.1f ", m->data[i*m->r+j]);
        putchar('\n');
    }
    putchar('\n');
}

void print_csr_matrix(csr_matrix cm)
{
    
}

void print_vector(vector v)
{
    unsigned i;
    for (i = 0; i < v->dim; ++i)
        printf("%.1f\n", v->data[i]);
    putchar('\n');
}

/* matrix/vector multiplications */
vector mmul(matrix m, vector v)
{
    unsigned i, j;
    vector r;
    
    if (m->r != v->dim) 
    {
        perror("ERROR - mmul() dimension mismatch: m[%u,%u] x v[%u]");
        exit(EXIT_FAILURE);
    }
    
    r = new_vector(v->dim);
    for (i = 0; i < m->r; ++i)
    {
        r->data[i] = 0.f;
        for (j = 0; j < m->c; ++j)
            r->data[i] += m->data[i * m->r + j] * v->data[j];
    }
    
    return r;
}

vector smmul(csr_matrix m, vector v)
{
    return NULL;
}