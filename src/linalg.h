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
void delete_csr_matrix(csr_matrix *pm);
void delete_vector(vector *pv);

void print_matrix(matrix m);
void print_csr_matrix(csr_matrix m);
void print_csr_matrix_(csr_matrix m);
void print_vector(vector v);

unsigned get_csr_row(csr_matrix m);

vector vsum(vector v1, vector v2);
vector svmul(float f, vector v);
float dot(vector v1, vector v2);
float dist(vector v1, vector v2);
vector mmul(matrix m, vector v);
vector smmul(csr_matrix m, vector v);

void fill_vector(vector v, float val);
void fill_row(matrix m, unsigned row, float val);
void fill_csr_row(csr_matrix m, unsigned row, float val);

vector find_dandling_nodes(csr_matrix m);

#endif