#ifndef LINALG_H
#define LINALG_H 

typedef struct _matrix     *matrix;
typedef struct _csr_matrix *csr_matrix;
typedef struct _vector     *vector;

matrix new_matrix(unsigned row, unsigned col);
matrix new_matrix_(float *v, unsigned row, unsigned col);
csr_matrix new_csr_matrix(unsigned row, unsigned col);
csr_matrix new_csr_matrix_(float *v, unsigned *ci, unsigned *rp, unsigned row, unsigned col);
vector new_vector(unsigned dim);
vector new_vector_(float *v, unsigned dim);

void delete_matrix(matrix *pm);
void delete_csr_matrix(csr_matrix *pcm);
void delete_vector(vector *pv);

void print_matrix(matrix m);
void print_csr_matrix(csr_matrix cm);
void print_csr_matrix_(csr_matrix cm);
void print_vector(vector v);

vector svmul(float f, vector v);
float dot(vector v1, vector v2);
vector mmul(matrix m, vector v);
vector smmul(csr_matrix m, vector v);

void fill_row(matrix m, unsigned row, float val);
void fill_csr_row(csr_matrix m, unsigned row, float val);

#endif