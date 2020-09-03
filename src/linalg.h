#ifndef LINALG_H
#define LINALG_H 

typedef struct _matrix     *matrix;
typedef struct _csr_matrix *csr_matrix;
typedef struct _vector     *vector;

void delete_matrix(matrix *pm);
void delete_csr_matrix(csr_matrix *pcm);
void delete_vector(vector *pv);

#endif