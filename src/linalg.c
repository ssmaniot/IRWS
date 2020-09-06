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

void delete_csr_matrix(csr_matrix *pcm)
{
    if (*pcm == NULL) return;
    _reset_ptr((void **) &((*pcm)->data));
    _reset_ptr((void **) &((*pcm)->col_ind));
    _reset_ptr((void **) &((*pcm)->row_ptr));
    _reset_ptr((void **) pcm);
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
void print_csr_matrix(csr_matrix cm)
{
    unsigned i, j, k;
    for (i = 0; i < cm->r; ++i)
    {
        k = cm->row_ptr[i];
        j = 0;
        for (j = 0; j < cm->c; ++j)
            if (k < cm->row_ptr[i+1] && j == cm->col_ind[k])
            {
                printf("%.1f ", cm->data[k]);
                ++k;
            }
            else 
                printf("0.0 ");
        putchar('\n');
    }
    putchar('\n');
}

void print_csr_matrix_(csr_matrix cm)
{
    unsigned i;
    
    printf("cm->data    = {");
    for (i = 0; i < cm->row_ptr[cm->r]; ++i)
    {
        printf("%.1f", cm->data[i]);
        if (i < cm->row_ptr[cm->r] - 1)
            printf(", ");
    }
    printf("}\n");
    printf("cm->col_ind = {");
    for (i = 0; i < cm->row_ptr[cm->r]; ++i)
    {
        printf("%2u", cm->col_ind[i]);
        if (i < cm->row_ptr[cm->r] - 1)
            printf(", ");
    }
    printf("}\n");
    printf("cm->row_ptr = {");
    for (i = 0; i <= cm->r; ++i)
    {
        printf("%2u", cm->row_ptr[i]);
        if (i < cm->r)
            printf(", ");
    }
    printf("}\n");
    printf("cm->r = %2u\n", cm->r);
    printf("cm->c = %2u\n", cm->c);
    putchar('\n');
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

/* Row-filler functions */
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
    
    nsize = (val == 0.f) ? pc + m->c + nc : pc + nc;
    
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

void fill_csr_row_(csr_matrix m, unsigned row, float val)
{
    float *data = NULL;
    unsigned *col_ind = NULL;
    unsigned len, diff;
    unsigned i;
    
    printf("fill_csr_row\n");
    /* compute the number of elements in the row */
    diff = m->row_ptr[row+1] - m->row_ptr[row];
    /* compute the length of the data array after the deletion of the elements in the row */
    len = m->row_ptr[m->r] - diff;
    
    /* allocate the new arrays for data and column index */
    data = (float *) malloc(sizeof(float) * len);
    col_ind = (unsigned *) malloc(sizeof(unsigned) * len);
    
    /* fill the data before the row to be edited */
    for (i = 0; i < m->row_ptr[row]; ++i)
    {
        data[i] = m->data[i];
        col_ind[i] = m->col_ind[i];
    }
    
    /* if val = 0, we copy the elements after the edited row and shifts the row pointers */
    if (val == 0.f)
    {
        printf("case val = 0.f\n");
        for (i = m->row_ptr[row+1]; i < m->row_ptr[m->r]; ++i)
        {
            data[i] = m->data[i+diff];
            col_ind[i] = m->col_ind[i];
            m->row_ptr[i] -= diff;
        }
        for (i = row + 1; i <= m->r; ++i)
            m->row_ptr[i] -= diff;
    }
    else 
    {
        printf("case val != 0.f\n");
        /* first, we insert the new values on the row */
        printf("first, we insert the new values on the row\n");
        for (i = 0; i < m->r; ++i)
        {
            data[m->row_ptr[row]+i] = val;
            col_ind[m->row_ptr[row]+i] = i;
        }
        /* then, we copy the data from row+1 on */
        diff = m->c - diff;
        printf("then, we copy the data from row+1 on");
        for (i = m->row_ptr[row+1]; i < m->row_ptr[m->r]; ++i)
        {
            data[i+diff] = data[i];
            col_ind[i+diff] = col_ind[i] + diff;
        }
        /* finally, we shift the row pointers */
        printf("finally, we shift the row pointers");
        for (i = row + 1; i <= m->r; ++i)
            m->row_ptr[i] += diff;
    }
    
    _reset_ptr((void **) (&(m->data)));
    m->data = data;
    _reset_ptr((void **) (&(m->col_ind)));
    m->col_ind = col_ind;
}

/* TODO: fix this! When substitute with val=0.f there are issues. */
void fill_csr_row__(csr_matrix m, unsigned row, float val)
{
    float *data;
    unsigned *col_ind;
    unsigned len, diff;
    unsigned i, j;
    
    if (val == 0.f)
    {
        len = m->row_ptr[m->r] - m->row_ptr[row+1] + m->row_ptr[row];
        col_ind = (unsigned *) malloc(sizeof(unsigned) * len);
        diff = m->row_ptr[row] - m->row_ptr[row+1];
    }
    else 
    {
        len = m->row_ptr[m->r] - m->row_ptr[row+1] + m->row_ptr[row] + m->c;
        col_ind = (unsigned *) malloc(sizeof(unsigned) * len);
        diff = m->row_ptr[row] - m->row_ptr[row+1] + m->c;
    }
    data = (float *) malloc(sizeof(float) * len);
    
    /* copy data and indexes before row */
    for (j = 0; j < m->row_ptr[row]; ++j)
    {
        data[j] = m->data[j];
        col_ind[j] = m->col_ind[j];
    }
    
    /* copy data and shifted indexes after row */
    for (j = m->row_ptr[row+1]; j < m->row_ptr[m->r]; ++j)
    {
        data[j + diff] = m->data[j];
        col_ind[j + diff] = m->col_ind[j] + diff;
    }
    
    /* update row pointers after row */
    for (i = row + 1; i <= m->r; ++i)
        m->row_ptr[i] += diff;
    
    /* insert new value if new value is nonzero */
    if (val != 0.f)
        for (j = m->row_ptr[row]; j < m->row_ptr[row+1]; ++j)
        {
            data[j] = val;
            col_ind[j] = j - m->row_ptr[row];
        }
    
    _reset_ptr((void **) (&(m->data)));
    m->data = data;
    _reset_ptr((void **) (&(m->col_ind)));
    m->col_ind = col_ind;
}