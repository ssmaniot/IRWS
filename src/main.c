#include "linalg.h"
#include <stdlib.h>

int main(int argc, char *argv[])
{
    float *md, *vd;
    unsigned *ci, *rp;
    csr_matrix m;
    vector v, r;
    
    md = (float *) malloc(sizeof(float) * 2);
    vd = (float *) malloc(sizeof(float) * 2);
    ci = (unsigned *) malloc(sizeof(unsigned) * 2);
    rp = (unsigned *) malloc(sizeof(unsigned) * 3);
    
    md[0] = 1.f; md[1] = 1.f; 
    ci[0] =   1; ci[1] =   0;
    rp[0] =   0; rp[1] =   1; rp[2] = 2;
    vd[0] = 1.f; vd[1] = 2.f;
    
    m = new_csr_matrix_(md, ci, rp, 2, 2);
    v = new_vector_(vd, 2);
    
    print_csr_matrix(m);
    print_vector(v);
    
    r = smmul(m, v);
    print_vector(r);
    
    fill_csr_row(m, 0, 2.f);
    print_csr_matrix(m);
    print_vector(v);
    delete_vector(&r);
    r = smmul(m, v);
    print_vector(r);
    
    delete_csr_matrix(&m);
    delete_vector(&v);
    delete_vector(&r);
    
    return EXIT_SUCCESS;
}