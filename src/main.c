#include "linalg.h"
#include <stdlib.h>

int main(int argc, char *argv[])
{
    float *md, *vd;
    matrix m;
    vector v, r;
    
    md = (float *) malloc(sizeof(float) * 4);
    vd = (float *) malloc(sizeof(float) * 2);
    
    md[0] = 0.f; md[1] = 1.f; md[2] = 1.f; md[3] = 0.f;
    vd[0] = 1.f;              vd[1] = 2.f;
    
    m = new_matrix_(md, 2, 2);
    v = new_vector_(vd, 2);
    
    print_matrix(m);
    print_vector(v);
    
    r = mmul(m, v);
    print_vector(r);
    
    delete_matrix(&m);
    delete_vector(&v);
    delete_vector(&r);
    
    return EXIT_SUCCESS;
}