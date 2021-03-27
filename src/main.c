#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define DIM 3
#define MAXITER 100

void print_matrix(double **A, int n);
void print_vector(double *v, int n);

int main(int argc, char * argv[])
{
    int row_id, col_id;
    int i, j;
    int iter;
    double d, error;
    double row_sum, err_sum;
    double **A, **T, **M;
    double *p, *p_new;
    FILE *fp;
    char *s = NULL;
    size_t slen = 0;
    ssize_t bytes;
    int no_nodes, no_edges;
    int lines;
    
    if (argc != 2)
    {
        printf("[ERROR] No input file given.\n");
        exit(1);
    }
    
    fp = fopen(argv[1], "r");
    if (fp == NULL)
    {
        printf("[ERROR] cannot open input file \"%s\".\n", argv[1]);
        exit(1);
    }
    
    getline(&s, &slen, fp);
    getline(&s, &slen, fp);
    getline(&s, &slen, fp);
    sscanf(s, "# Nodes: %d Edges: %d", &no_nodes, &no_edges);
    printf("This graph has %d nodes and %d edges\n", no_nodes, no_edges);
    A = (double **) malloc(sizeof(double *) * no_nodes);
    for (i = 0; i < no_nodes; ++i)
        A[i] = (double *) calloc(no_nodes, sizeof(double));
    bytes = getline(&s, &slen, fp);
    lines = 0;
    while ((bytes = getline(&s, &slen, fp)) != -1)
    {
        sscanf(s, "%d %d", &i, &j);
        A[i][j] = 1.;
        lines += 1;
        if (lines % 10 == 0)
            printf("Lines %d/%d\r", lines, no_edges);
    }
    
    fclose(fp);
    free(s);
    
    printf("Lines %d/%d\n", lines, no_edges);
    printf("number of nodes in graph A: %d\n", no_nodes);
    
    /* Fix dangling nodes */
    for (row_id = 0; row_id < no_nodes; ++row_id)
    {
        row_sum = 0.;
        for (col_id = 0; col_id < no_nodes; ++col_id)
            row_sum += A[row_id][col_id];
        if (row_sum == 0.)
            for (col_id = 0; col_id < no_nodes; ++col_id)
                A[row_id][col_id] = 1.;
    }
    
    /* Normalize, making the matrix stochastic */
    for (row_id = 0; row_id < no_nodes; ++row_id)
    {
        row_sum = 0.;
        for (col_id = 0; col_id < no_nodes; ++col_id)
            row_sum += A[row_id][col_id];
        for (col_id = 0; col_id < no_nodes; ++col_id)
            A[row_id][col_id] /= row_sum;
    }
    
    /* Initialize the teleportation matrix T */
    T = (double **) malloc(sizeof(double *) * no_nodes);
    for (i = 0; i < no_nodes; ++i)
        T[i] = (double *) malloc(sizeof(double) * no_nodes);
    
    for (row_id = 0; row_id < no_nodes; ++row_id)
        for (col_id = 0; col_id < no_nodes; ++col_id)
            T[row_id][col_id] = 1. / (double) no_nodes;
    
    /* Markov matrix */
    /* Final Matrix after removing dangling nodes and adding teleportation */
    M = (double **) malloc(sizeof(double *) * no_nodes);
    for (i = 0; i < no_nodes; ++i)
        M[i] = (double *) malloc(sizeof(double) * no_nodes);
    d = 0.85;
    for (row_id = 0; row_id < no_nodes; ++row_id)
        for (col_id = 0; col_id < no_nodes; ++col_id)
            M[row_id][col_id] = d * A[row_id][col_id] + (1. - d) * T[row_id][col_id];
    
    /* Initial probability distribution */
    p = (double *) malloc(sizeof(double) * no_nodes);
    for (i = 0; i < no_nodes; ++i)
        p[i] = 1. / (double) no_nodes;
    
    p_new = (double *) malloc(sizeof(double) * no_nodes);
    
    for (iter = 1; iter <= MAXITER; ++iter)
    {
        for (i = 0; i < no_nodes; ++i)
        {
            p_new[i] = 0.;
            for (j = 0; j < no_nodes; ++j)
                p_new[i] += M[j][i] * p[j];
        }
        
        if (iter % 10 == 0)
            printf("iter: %d\r", iter);
        
        err_sum = 0.;
        for (i = 0; i < DIM; ++i)
            err_sum += (p[i] - p_new[i]) * (p[i] - p_new[i]);
        error = sqrt(err_sum);
        
        for (i = 0; i < DIM; ++i)
            p[i] = p_new[i];
        
        if (error < 10.e-10)
        {
            printf("early exit at iteration: %d\n", iter);
            break;
        }
    }
    
    if (iter == MAXITER)
        printf("iter: %d\n", iter);
    
    print_vector(p, no_nodes);
    
    for (i = 0; i < no_nodes; ++i)
    {
        free(A[i]); free(T[i]); free(M[i]);
    }
    free(A); free(T); free(M);
    free(p); free(p_new);
    
    return EXIT_SUCCESS;
}

/* Auxiliary function implementations */

void print_matrix(double **A, int n)
{
    int i, j;
    printf("[");
    for (i = 0; i < n; ++i)
    {
        printf("[ ");
        for (j = 0; j < n; ++j)
            printf("%.8f ", A[i][j]);
        printf("]");
        if (i < n - 1)
            printf("\n ");
    }
    printf("]\n");
}

void print_vector(double *v, int n)
{
    int i;
    printf("[[ ");
    for (i = 0; i < n; ++i)
        printf("%.8f ", v[i]);
    printf("]]\n");
}