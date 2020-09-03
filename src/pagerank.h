#ifndef PAGERANK_H
#define PAGERANK_H

#include "linalg.h"

vector pagerank(matrix A);
vector csr_pagerank(csr_matrix A);
vector opt_pagerank(matrix A);
vector opt_csr_pagerank(csr_matrix A);

#endif