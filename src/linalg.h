#ifndef LINALG_H
#define LINALG_H 

struct matrix;
struct csr_matrix;
struct vector;

typedef struct matrix     *matrix_ptr;
typedef struct csr_matrix *csr_matrix_ptr;
typedef struct vector     *vector_ptr;

void delete(matrix_ptr     *pm );
void delete(csr_matrix_ptr *pcm);
void delete(vector_ptr     *pv );

#endif