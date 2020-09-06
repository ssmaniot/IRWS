#include "linalg.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

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

csr_matrix new_csr_matrix_(float *v, unsigned *ci, unsigned *rp, unsigned row, unsigned col)
{
    csr_matrix m = (csr_matrix) malloc(sizeof(struct _csr_matrix));
    m->data = v;
    m->col_ind = ci;
    m->row_ptr = rp;
    m->r = row;
    m->c = col;
    return m;
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
static inline void _reset_ptr(void **ptr)
{
    if (*ptr == NULL) return;
    free(*ptr);
    *ptr = NULL;
}

/* delete functions */
void delete_matrix(matrix *pm)
{
    if (*pm == NULL) return;
    _reset_ptr((void **) &((*pm)->data));
    _reset_ptr((void **) pm);
}

void delete_csr_matrix(csr_matrix *pm)
{
    if (*pm == NULL) return;
    _reset_ptr((void **) &((*pm)->data));
    _reset_ptr((void **) &((*pm)->col_ind));
    _reset_ptr((void **) &((*pm)->row_ptr));
    _reset_ptr((void **) pm);
}

void delete_vector(vector *pv)
{
    if (*pv == NULL) return;
    _reset_ptr((void **) &((*pv)->data));
    _reset_ptr((void **) pv);
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

/* TODO: implement sparse matrix printer */
void print_csr_matrix(csr_matrix m)
{
    unsigned i, j, k;
    for (i = 0; i < m->r; ++i)
    {
        k = m->row_ptr[i];
        j = 0;
        for (j = 0; j < m->c; ++j)
            if (k < m->row_ptr[i+1] && j == m->col_ind[k])
            {
                printf("%.1f ", m->data[k]);
                ++k;
            }
            else 
                printf("0.0 ");
        putchar('\n');
    }
    putchar('\n');
}

void print_csr_matrix_(csr_matrix m)
{
    unsigned i;
    
    printf("m->data    = {");
    for (i = 0; i < m->row_ptr[m->r]; ++i)
    {
        printf("%.1f", m->data[i]);
        if (i < m->row_ptr[m->r] - 1)
            printf(", ");
    }
    printf("}\n");
    printf("m->col_ind = {");
    for (i = 0; i < m->row_ptr[m->r]; ++i)
    {
        printf("%2u", m->col_ind[i]);
        if (i < m->row_ptr[m->r] - 1)
            printf(", ");
    }
    printf("}\n");
    printf("m->row_ptr = {");
    for (i = 0; i <= m->r; ++i)
    {
        printf("%2u", m->row_ptr[i]);
        if (i < m->r)
            printf(", ");
    }
    printf("}\n");
    printf("m->r = %2u\n", m->r);
    printf("m->c = %2u\n", m->c);
    putchar('\n');
}

void print_vector(vector v)
{
    unsigned i;
    for (i = 0; i < v->dim; ++i)
        printf("%.1f\n", v->data[i]);
    putchar('\n');
}

unsigned get_csr_row(csr_matrix m)
{
    return m->r;
}

/* Operations on structures */
vector vsum(vector v1, vector v2)
{
    vector r;
    unsigned i;

    r = new_vector(v1->dim);
    for (i = 0; i < v1->dim; ++i)
        r->data[i] = v1->data[i] + v2->data[i];

    return r;
}

/* scalar x vector */
vector svmul(float f, vector v)
{
    vector r;
    unsigned i;

    r = new_vector(v->dim);
    for (i = 0; i < r->dim; ++i)
        r->data[i] = v->data[i] * f;

    return r;
}

/* scalar product between vectors */
float dot(vector v1, vector v2)
{
    float res = 0.f;
    unsigned i;

    for (i = 0; i < v1->dim; ++i)
        res += v1->data[i] * v2->data[i];

    return res;
}

/* L1-distance of v1 and v2 */
float dist(vector v1, vector v2)
{
    float d = 0.0, di;
    unsigned i;

    for (i = 0; i < v1->dim; ++i)
    {
        di = v1->data - v2->data;
        d += (di < 0) ? -di : di;
    }

    return d;
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
        j = m->row_ptr[i];
        while (j < m->row_ptr[i+1])
        {
            r->data[i] += m->data[j] * v->data[m->col_ind[j]];
            ++j;
        }
    }
    
    return r;
}

/* Filler functions */
void fill_vector(vector v, float val)
{
    unsigned i;

    for (i = 0; i < v->dim; ++i)
        v->data[i] = val;
}

void fill_row(matrix m, unsigned row, float val)
{
    unsigned col;
    for (col = 0; col < m->c; ++col)
        m->data[row * m->c + col] = val;
}

void fill_csr_row(csr_matrix m, unsigned row, float val)
{
    float *data;
    unsigned *col_ind;
    unsigned pc, bc, nc;
    unsigned nsize;
    unsigned i;
    
    pc = m->row_ptr[row];
    bc = m->row_ptr[row+1] - m->row_ptr[row];
    nc = m->row_ptr[m->r] - bc;
    
    nsize = (val == 0.f) ? pc + nc : pc + m->c + nc;
    
    data = (float *) malloc(sizeof(float) * nsize);
    col_ind = (unsigned *) malloc(sizeof(unsigned) * nsize);

    for (i = 0; i < m->row_ptr[row]; ++i)
    {
        data[i] = m->data[i];
        col_ind[i] = m->col_ind[i];
    }
    
    if (val == 0.f)
    {
        for (i = m->row_ptr[row+1]; i < m->row_ptr[m->r]; ++i)
        {
            data[i-bc] = m->data[i];
            col_ind[i-bc] = m->col_ind[i];
        }
        for (i = row + 1; i <= m->r; ++i)
            m->row_ptr[i] -= bc;
    }
    else 
    {
        for (i = 0; i < m->c; ++i)
        {
            data[m->row_ptr[row]+i] = val;
            col_ind[m->row_ptr[row]+i] = i;
        }
        for (i = m->row_ptr[row+1]; i < m->row_ptr[m->r]; ++i)
        {
            data[i+m->c-bc] = m->data[i];
            col_ind[i+m->c-bc] = m->col_ind[i];
        }
        for (i = row + 1; i <= m->r; ++i)
            m->row_ptr[i] += m->c - bc;
    }
    
    _reset_ptr((void **) (&(m->data)));
    m->data = data;
    _reset_ptr((void **) (&(m->col_ind)));
    m->col_ind = col_ind;
}

void normalize_csr_by_row(csr_matrix m)
{
    float sum;
    unsigned i, j;

    for (i = 0; i < m->r; ++i)
    {
        sum = 0.f;
        for (j = m->row_ptr[i]; j < m->row_ptr[i+1]; ++i)
            sum += m->data[i];
        if (sum != 0.f)
            for (j = m->row_ptr[i]; j < m->row_ptr[i+1]; ++i)
                m->data[j] /= sum;
    }
}

vector find_dandling_nodes(csr_matrix m)
{
    vector d;
    unsigned i;

    d = new_vector(m->r);
    for (i = 0; i < m->r; ++i)
        d->data[i] = (m->row_ptr[i] == m->row_ptr[i+1]) ? 0.f : 1.f;

    return d;
}