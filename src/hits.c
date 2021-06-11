#include <stdio.h>
#include <stdlib.h>
#define DEBUG 

/* Data for compression */
typedef struct 
{
    int no_nodes;
    int no_edges;
} CSR_data;

int main(int argc, char *argv[])
{
    /* Data to save/load CSR matrix */
    FILE *pdata;
    char dir[DNAME];
    char row_ptr_p[2][PATH];
    char col_ind_p[2][PATH];
    char csr_data_p[PATH];
    CSR_data csr_data;
    struct stat st = {0};
    
    /* Reading data from input file */
    FILE *pf;
    int no_nodes, no_edges;
    char *s = NULL;
    size_t slen = 0;
    ssize_t bytes;
    int ri, ci;
    int *from, *to;
    int f, t;
    int i, j;
    
    /* CSR matrix representation */
    double *val;
    int *col_ind;
    int *row_ptr;
    
    /* Pagerank computation data */
    double danglings_dot_product;
    double *p, *p_new;
    double d;
    double dist;
    int iter;
    
    /* Extra data */
    double sum;
    int err;
    
    if (argc != 2)
    {
        fprintf(stderr, " [ERROR] *1* argument required: ./pagerank <arg_name>\n");
        exit(EXIT_FAILURE);
    }
    
    /* Check if input data has already been compressed */
    strcpy(dir, "PR_");
    strncpy(dir + 3, argv[1] + 5, strlen(argv[1])-9);
    dir[strlen(argv[1])-6] = '/';
    dir[strlen(argv[1])-5] = '\0';
    printf("%s\n", dir);
    
    /* Create CSR file names */
    strcpy(row_ptr_p, dir);
    strcat(row_ptr_p, "row_ptr.bin");
    strcpy(col_ind_p, dir);
    strcat(col_ind_p, "col_ind.bin");
    strcpy(val_p, dir);
    strcat(val_p, "val.bin");
    strcpy(danglings_p, dir);
    strcat(danglings_p, "danglings.bin");
    strcpy(csr_data_p, dir);
    strcat(csr_data_p, "csr_data.bin");
    
    /* If data has NOT yet been compressed, then perform compression */
    if (stat(dir, &st) == -1)
    {
        printf("Input file data \"%s\" is not compressed, ready to perform compression...\n\n", argv[1]);
        mkdir(dir, 0700);
    
        if ((pf = fopen(argv[1], "r")) == NULL)
        {
            fprintf(stderr, " [ERROR] cannot open input file \"%s\"\n", argv[1]);
            exit(EXIT_FAILURE);
        }
        
        /* Parsing input file header */
        printf("Parsing input data...\n");
        bytes = getline(&s, &slen, pf);
        bytes = getline(&s, &slen, pf);
        bytes = getline(&s, &slen, pf);
        sscanf(s, "# Nodes: %d Edges: %d", &no_nodes, &no_edges);
        printf("This graph has %d nodes and %d edges\n", no_nodes, no_edges);
        bytes = getline(&s, &slen, pf);
        
        csr_data.no_nodes = no_nodes;
        csr_data.no_edges = no_edges;
        
        /* Reading data from input file */
        i = 0;
        from = (int *) malloc(sizeof(int) * no_edges);
        to = (int *) malloc(sizeof(int) * no_edges);
        out_links = (int *) calloc(no_nodes, sizeof(int));
        while ((bytes = getline(&s, &slen, pf)) != -1)
        {
            sscanf(s, "%d %d", from + i, to + i);
            out_links[from[i]] += 1;         
            if (i % 10 == 0)
                printf("\rEdge %d/%d", i, no_edges);
            ++i;
        }
        printf("\rEdge %d/%d\n", i, no_edges);
        printf("Done\n\n");
        fclose(pf);
        free(s);
        
        /* Keeping track of danglings data */
        no_danglings = 0;
        for (i = 0; i < no_nodes; ++i)
            if (out_links[i] == 0)
                ++no_danglings;
        danglings = (int *) malloc(sizeof(int) * no_danglings);
        j = 0;
        for (i = 0; i < no_nodes; ++i)
            if (out_links[i] == 0)
                danglings[j++] = i;
        
        csr_data.no_danglings = no_danglings;

        printf("Sorting edges...\n");
        sort_input_data(from, to, no_edges);
        printf("Done.\n\n");
        
        /* csr matrix initialization */
        val = (double *) malloc(sizeof(double) * no_edges);
        col_ind = (int *) malloc(sizeof(int) * no_edges);
        row_ptr = (int *) malloc(sizeof(int) * (no_nodes + 1));
        ri = 0;
        row_ptr[ri] = 0;
        ci = 0;
        
        /* Writing data in CSR matrix */
        for (ci = 0; ci < no_edges; ++ci)
        {
            f = from[ci];
            t = to[ci];
            if (t > ri)
            {
                for (i = ri + 1; i <= t; ++i)
                    row_ptr[i] = ci;
                ri = t;
            }
            val[ci] = 1. / (double) out_links[f];
            col_ind[ci] = f;
        }
        while (ri < no_nodes)
            row_ptr[++ri] = ci;

#ifdef DEBUG 
        printf("CSR Transposed matrix\n");
        printf("---------------------\n");
        printf("val:     [ ");
        for (i = 0; i < no_edges; ++i)
            printf("%.3f ", val[i]);
        printf("]\n");
        
        printf("col_ind: [ ");
        for (i = 0; i < no_edges; ++i)
            printf("%d ", col_ind[i]);
        printf("]\n");
        
        printf("row_ptr: [ ");
        for (i = 0; i < no_nodes+1; ++i)
            printf("%d ", row_ptr[i]);
        printf("]\n\n");
        printf("danglings: [ ");
        for (j = 0; j < no_danglings; ++j)
        {
            printf("%d", danglings[j]);
            if (j < no_danglings - 1)
                printf(", ");
        }
        printf(" ]\n");
        printf("Number of danglings nodes: %d\n\n", no_danglings);
#endif 
    
        /* Writing data back to memory */
        err = (write_data(row_ptr_p, (void *) row_ptr, sizeof(int), no_nodes + 1) == EXIT_FAILURE) 
           || (write_data(col_ind_p, (void *) col_ind, sizeof(int), no_edges) == EXIT_FAILURE)
           || (write_data(val_p, (void *) val, sizeof(double), no_edges) == EXIT_FAILURE)
           || (write_data(danglings_p, (void *) danglings, sizeof(int), no_danglings) == EXIT_FAILURE)
           || (write_data(csr_data_p, (void *) &csr_data, sizeof(CSR_data), 1) == EXIT_FAILURE);
        
        /* Input data */
        free(from); free(to);
        from = NULL; to = NULL;
        /* Danglings data */
        free(out_links); free(danglings);
        out_links = NULL; danglings = NULL;
        /* CSR data structure */
        free(val); free(col_ind); free(row_ptr);
        val = NULL; col_ind = NULL; row_ptr = NULL;
        
        /* Manage error from writing data to memory */
        if (err)
        {
            delete_folder(dir);
            fprintf(stderr, " [ERROR] data could not be written in memory.\n");
            exit(EXIT_FAILURE);
        }
    }
    
    printf("Reading csr matrix data...\n");
    /* Reading CSR matrix metadata info from file */
    pdata = fopen(csr_data_p, "rb");
    bytes = fread(&no_nodes, sizeof(csr_data.no_nodes), 1, pdata);
    bytes = fread(&no_edges, sizeof(csr_data.no_edges), 1, pdata);
    bytes = fread(&no_danglings, sizeof(csr_data.no_danglings), 1, pdata);
    fclose(pdata);
    printf("no_nodes: %d\nno_edges: %d\nno_danglings: %d\n", no_nodes, no_edges, no_danglings);
    
    /* mmapping the CSR matrix data from files */
    err = 0;
    if ((row_ptr = (int *) mmap_data(row_ptr_p, sizeof(int), no_nodes + 1)) == NULL)          ++i; 
    else if ((col_ind = (int *) mmap_data(col_ind_p, sizeof(int), no_edges)) == NULL)         ++i;
    else if ((val = (double *) mmap_data(val_p, sizeof(double), no_edges)) == NULL)           ++i;
    else if ((danglings = (int *) mmap_data(danglings_p, sizeof(int), no_danglings)) == NULL) ++i;
    
    if (err > 0)
    {
        fprintf(stderr, " [ERROR] data could not be mmapped from memory.\n");
        fprintf(stderr, "         Data is corrupted, the folder will be destroyed.\n");
        delete_folder(dir);
        /* Using Duff's device to manage errors */
        switch (err)
        {
            case 1: do { munmap(row_ptr, (no_nodes + 1) * sizeof(int));
            case 2:      munmap(col_ind, no_edges * sizeof(int));
            case 3:      munmap(val, no_edges * sizeof(double));
            case 4:      munmap(col_ind, no_nodes * sizeof(int));
                       } while (--err > 0);
        }
        exit(EXIT_FAILURE);
    }
    
    printf("Done.\n\n");

#ifdef DEBUG 
        printf("CSR Transposed matrix\n");
        printf("---------------------\n");
        printf("val:     [ ");
        for (i = 0; i < no_edges; ++i)
            printf("%.3f ", val[i]);
        printf("]\n");
        
        printf("col_ind: [ ");
        for (i = 0; i < no_edges; ++i)
            printf("%d ", col_ind[i]);
        printf("]\n");
        
        printf("row_ptr: [ ");
        for (i = 0; i < no_nodes+1; ++i)
            printf("%d ", row_ptr[i]);
        printf("]\n\n");
        printf("danglings: [ ");
        for (j = 0; j < no_danglings; ++j)
        {
            printf("%d", danglings[j]);
            if (j < no_danglings - 1)
                printf(", ");
        }
        printf(" ]\n");
        printf("Number of danglings nodes: %d\n\n", no_danglings);
#endif 
    
    /* Setting data up for PageRank computation */
    d = 0.85;
    p = (double *) malloc(sizeof(double) * no_nodes);
    for (i = 0; i < no_nodes; ++i)
        p[i] = 1. / (double) no_nodes;
    p_new = (double *) malloc(sizeof(double) * no_nodes);
    dist = DBL_MAX;
    iter = 0;
    
    /* Computing PageRank */
    printf("Computing PageRank...\n");
    while (dist > TOL && iter < MAX_ITER)
    {
#ifdef DEBUG 
        if (iter % 10 == 0)
        {
#endif 
            printf("\riter %d", iter);
#ifdef DEBUG 
            printf("\n");
            print_p(p, no_nodes); 
        }
#endif 
        
        /* DTp = DanglingsT @ p */
        danglings_dot_product = 0.;
        for (j = 0; j < no_danglings; ++j)
            danglings_dot_product += p[danglings[j]];
        danglings_dot_product /= (double) no_nodes;
        
        /* ATp = AT @ p + DTp */
        for (ri = 0; ri < no_nodes; ++ri)
        {
            p_new[ri] = danglings_dot_product;
            for (ci = row_ptr[ri]; ci < row_ptr[ri+1]; ++ci)
                p_new[ri] += p[col_ind[ci]] * val[ci];
        }
        
        /* d*AT @ p + (1-d)eeT @ p */
        for (i = 0; i < no_nodes; ++i)
            p_new[i] = d * p_new[i] + (1. - d) / (double) no_nodes;
        
        dist = 0.;
        for (i = 0; i < no_nodes; ++i)
            dist += (p[i] - p_new[i]) * (p[i] - p_new[i]);
        dist = sqrt(dist);
        
        for (i = 0; i < no_nodes; ++i)
            p[i] = p_new[i];
        
        ++iter;
    }
    printf("\riter %d\n", iter);
#ifdef DEBUG 
    print_p(p, no_nodes);
#endif 
    printf("Done.\n\n");
    
    sum = 0.;
    for (i = 0; i < no_nodes; ++i)
        sum += p[i];
    printf("Proof of correctness:\n");
    printf("sum(p) = %f\n\n", sum);
    
    /* un-mmapping data */
    munmap(row_ptr, (no_nodes + 1) * sizeof(int));
    munmap(col_ind, no_edges * sizeof(int));
    munmap(val, no_edges * sizeof(double));
    munmap(col_ind, no_nodes * sizeof(int));
    
    /* Vectors of probability */
    free(p); free(p_new);
    
    exit(EXIT_SUCCESS);
}