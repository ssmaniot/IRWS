#include "pagerank.h"
#include <stdlib.h>

vector pagerank(matrix A, float d)
{
    return NULL;
}

vector csr_pagerank(csr_matrix At, float d)
{
    vector dandlings, coef, p, pk, Atp;
    vector lhs, rhs;
    float etp;
    unsigned n = get_csr_row(At);

    dandlings = find_dandling_nodes(At);
    coef = new_vector(n);
    fill_vector(coef, (1.f - d) / (float) n);
    pk = new_vector(n);
    fill_vector(p, 1.f / (float) n);
    normalize_csr_by_row(At);

    do 
    {
        delete_vector(&p);
        p = pk;

        Atp = smmul(At, p);
        lhs = svmul(d, Atp);
        delete_vector(&Atp);

        etp = dot(dandlings, p);
        rhs = svmul(etp, coef);

        pk = vsum(rhs, lhs);
        delete_vector(&lhs);
        delete_vector(&rhs);
    } while(dist(p, pk) >= 1.e-6f);

    delete_vector(&dandlings);
    delete_vector(&p);
    delete_vector(&coef);

    return pk;
}

vector opt_pagerank(matrix A, float d)
{
    return NULL;
}

vector opt_csr_pagerank(csr_matrix A, float d)
{
    return NULL;
}