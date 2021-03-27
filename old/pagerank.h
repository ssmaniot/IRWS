#ifndef PAGERANK_H
#define PAGERANK_H

#include "linalg.h"

vector pagerank(matrix At, float d);
vector csr_pagerank(csr_matrix At, float d);
vector opt_pagerank(matrix At, float d);
vector opt_csr_pagerank(csr_matrix At, float d);

#endif