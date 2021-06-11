#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>
#include <sys/mman.h>
#include <fcntl.h>

#define TOL 1.e-10
#define MAX_ITER 3
#define MOD_ITER 1
#define FNAME 256
#define DNAME 1024
#define PATH 1024
#define MMAP 2048
#define DEBUG

/* Data for compression */
typedef struct 
{
    int no_nodes;
    int no_edges;
    int no_danglings;
} LCSR_data;

/* Helper functions */
int write_data(char path[], void *data, size_t nmemb, size_t size);
void delete_folder(char dir[]);
void * mmap_data(char path[], size_t nmemb, size_t size);
void print_vec_f(double *v, int n);
void print_vec_d(int *v, int n);
void double_merge(int *from, int *to, int lo, int mid, int hi);
void double_merge_sort(int *from, int *to, int lo, int hi);
void sort_input_data(int *from, int *to, int n);

int main(int argc, char *argv[])
{
    /* Data to save/load LCSR matrix */
    FILE *pdata;
    char dir[DNAME];
    /*   Matrix L            Matrix L^T         */
    char row_ptr_p[PATH],    row_ptr_tp[PATH];
    char col_ind_p[PATH],    col_ind_tp[PATH];
    char val_p[PATH],        val_tp[PATH];
    char danglings_p[PATH],  danglings_tp[PATH];
    char lcsr_data_p[PATH];
    LCSR_data lcsr_data;
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
    
    /* LCSR matrix representation */
    double *val;
    int *col_ind, *col_ind_t;
    int *row_ptr, *row_ptr_t;
    
    /* HITS computation data */
    double *a, *a_new;
    double *h, *h_new;
    double a_dist, h_dist;
    
    /* PageRank computation data */
    int *danglings;
    int no_danglings;
    double danglings_dot_product;
    double *p, *p_new;
    double d;
    double dist;
    int iter;
    
    /* Extra data */
    int *out_links;
    double sum;
    int err;
    
    if (argc != 2)
    {
        fprintf(stderr, " [ERROR] *1* argument required: ./pagerank <arg_name>\n");
        exit(EXIT_FAILURE);
    }
    
    /* Init data folder name */
    strcpy(dir, "HITS_");
    strncpy(dir + 5, argv[1] + 5, strlen(argv[1])-9);
    dir[strlen(argv[1])-4] = '/';
    dir[strlen(argv[1])-3] = '\0';
    printf("%s\n", dir);
    
    /* Create LCSR file names */
    strcpy(row_ptr_p,   dir); strcat(row_ptr_p,   "row_ptr.bin");
    strcpy(col_ind_p,   dir); strcat(col_ind_p,   "col_ind.bin");
    /*strcpy(val_p,       dir); strcat(val_p,       "val.bin");*/
    /*strcpy(danglings_p, dir); strcat(danglings_p, "danglings.bin");*/
    
    /* Create transposed LCSR file names */
    strcpy(row_ptr_tp,   dir); strcat(row_ptr_tp,   "row_ptr_t.bin");
    strcpy(col_ind_tp,   dir); strcat(col_ind_tp,   "col_ind_t.bin");
    /*strcpy(val_tp,       dir); strcat(val_tp,       "val_t.bin");*/
    /*strcpy(danglings_tp, dir); strcat(danglings_tp, "danglings_t.bin");*/
    
    /* Create LCSR metadata file */
    strcpy(lcsr_data_p,  dir); strcat(lcsr_data_p,  "lcsr_data.bin");
    
    /* Check if input data has already been compressed. 
     * If data has NOT yet been compressed, then perform compression */
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
        
        lcsr_data.no_nodes = no_nodes;
        lcsr_data.no_edges = no_edges;
        
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
        
        /* LCSR matrix initialization */
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
                ri = f;
            }
            col_ind[ci] = t;
        }
        while (ri < no_nodes)
            row_ptr[++ri] = ci;

#ifdef DEBUG 
        printf("CSR matrix\n");
        printf("---------------------\n");
        
        printf("col_ind: [ ");
        for (i = 0; i < no_edges; ++i)
            printf("%d ", col_ind[i]);
        printf("]\n");
        
        printf("row_ptr: [ ");
        for (i = 0; i < no_nodes+1; ++i)
            printf("%d ", row_ptr[i]);
        printf("]\n\n");
#endif 

        printf("Sorting edges for transposed matrix...\n");
        sort_input_data(from, to, no_edges);
        printf("Done.\n\n");
        
        /* Transposed LCSR matrix initialization */
        col_ind_t = (int *) malloc(sizeof(int) * no_edges);
        row_ptr_t = (int *) malloc(sizeof(int) * (no_nodes + 1));
        ri = 0;
        row_ptr_t[ri] = 0;
        ci = 0;
        
        /* Writing data in CSR matrix */
        for (ci = 0; ci < no_edges; ++ci)
        {
            f = from[ci];
            t = to[ci];
            if (t > ri)
            {
                for (i = ri + 1; i <= t; ++i)
                    row_ptr_t[i] = ci;
                ri = t;
            }
            col_ind_t[ci] = f;
        }
        while (ri < no_nodes)
            row_ptr_t[++ri] = ci;

#ifdef DEBUG 
        printf("Transposed CSR matrix\n");
        printf("---------------------\n");
        
        printf("col_ind_t: [ ");
        for (i = 0; i < no_edges; ++i)
            printf("%d ", col_ind_t[i]);
        printf("]\n");
        
        printf("row_ptr_t: [ ");
        for (i = 0; i < no_nodes+1; ++i)
            printf("%d ", row_ptr_t[i]);
        printf("]\n\n");
#endif 
    
        /* Writing data back to memory */
        err = (write_data(row_ptr_p, (void *) row_ptr, sizeof(int), no_nodes + 1) == EXIT_FAILURE) 
           || (write_data(col_ind_p, (void *) col_ind, sizeof(int), no_edges) == EXIT_FAILURE)
           || (write_data(row_ptr_tp, (void *) row_ptr_t, sizeof(int), no_nodes + 1) == EXIT_FAILURE) 
           || (write_data(col_ind_tp, (void *) col_ind_t, sizeof(int), no_edges) == EXIT_FAILURE)
         /*|| (write_data(val_p, (void *) val, sizeof(double), no_edges) == EXIT_FAILURE)
           || (write_data(danglings_p, (void *) danglings, sizeof(int), no_danglings) == EXIT_FAILURE)*/
           || (write_data(lcsr_data_p, (void *) &lcsr_data, sizeof(LCSR_data), 1) == EXIT_FAILURE);
        
        /* Input data */
        free(from); free(to);
        from = NULL; to = NULL;
        /* Danglings data */
        /*free(out_links); free(danglings);*/
        /*out_links = NULL; danglings = NULL;*/
        /* CSR data structure */
        /*free(val);*/ free(col_ind); free(row_ptr);
        free(col_ind_t); free(row_ptr_t);
        /*val = NULL;*/ col_ind = NULL; row_ptr = NULL;
        col_ind_t = NULL; row_ptr_t = NULL;
        
        /* Manage error from writing data to memory */
        if (err)
        {
            delete_folder(dir);
            fprintf(stderr, " [ERROR] data could not be written in memory.\n");
            exit(EXIT_FAILURE);
        }
    }
    
    /* Reading LCSR matrix metadata info from file */
    printf("Reading CLSR matrix data...\n");
    pdata = fopen(lcsr_data_p, "rb");
    bytes = fread(&no_nodes, sizeof(lcsr_data.no_nodes), 1, pdata);
    bytes = fread(&no_edges, sizeof(lcsr_data.no_edges), 1, pdata);
    /*bytes = fread(&no_danglings, sizeof(lcsr_data.no_danglings), 1, pdata);*/
    fclose(pdata);
    printf("no_nodes: %d\nno_edges: %d\n\n", no_nodes, no_edges);
    
    /* mmapping the CSR matrix data from files */
    err = 0;
    if      ((row_ptr   = (int *) mmap_data(row_ptr_p,  sizeof(int), no_nodes + 1)) == NULL) ++i; 
    else if ((row_ptr_t = (int *) mmap_data(row_ptr_tp, sizeof(int), no_nodes + 1)) == NULL) ++i; 
    else if ((col_ind   = (int *) mmap_data(col_ind_p,  sizeof(int), no_edges)) == NULL)     ++i;
    else if ((col_ind_t = (int *) mmap_data(col_ind_tp, sizeof(int), no_edges)) == NULL)     ++i;
    /*else if ((val = (double *) mmap_data(val_p, sizeof(double), no_edges)) == NULL)           ++i;
    else if ((danglings = (int *) mmap_data(danglings_p, sizeof(int), no_danglings)) == NULL) ++i;*/
    
    if (err > 0)
    {
        fprintf(stderr, " [ERROR] data could not be mmapped from memory.\n");
        fprintf(stderr, "         Data is corrupted, the folder will be destroyed.\n");
        delete_folder(dir);
        /* Using Duff's device to manage errors */
        switch (err)
        {
            case 1: do { munmap(row_ptr,   (no_nodes + 1) * sizeof(int));
            case 2:      munmap(row_ptr_t, (no_nodes + 1) * sizeof(int));
            case 3:      munmap(col_ind,   no_edges * sizeof(int));
            case 4:      munmap(col_ind_t, no_edges * sizeof(int));
                       } while (--err > 0);
        }
        exit(EXIT_FAILURE);
    }
    
    printf("Done.\n\n");

#ifdef DEBUG 
        printf("CSR Transposed matrix\n");
        printf("---------------------\n");
        
        printf("col_ind: [ ");
        for (i = 0; i < no_edges; ++i)
            printf("%d ", col_ind[i]);
        printf("]\n");
        
        printf("row_ptr: [ ");
        for (i = 0; i < no_nodes+1; ++i)
            printf("%d ", row_ptr[i]);
        printf("]\n\n");
        
        printf("Transposed CSR matrix\n");
        printf("---------------------\n");
        
        printf("col_ind_t: [ ");
        for (i = 0; i < no_edges; ++i)
            printf("%d ", col_ind_t[i]);
        printf("]\n");
        
        printf("row_ptr_t: [ ");
        for (i = 0; i < no_nodes+1; ++i)
            printf("%d ", row_ptr_t[i]);
        printf("]\n\n");
#endif 
    /* HITS computation data 
    double *a, *a_new;
    double *h, *h_new;
    double a_dist, h_dist;  */
    
    /* Setting data up for PageRank computation */
    a = (double *) malloc(sizeof(double) * no_nodes);
    h = (double *) malloc(sizeof(double) * no_nodes);
    for (i = 0; i < no_nodes; ++i)
    {
        a[i] = 1.;
        h[i] = 1.;
    }
    a_new = (double *) malloc(sizeof(double) * no_nodes);
    h_new = (double *) malloc(sizeof(double) * no_nodes);
    a_dist = DBL_MAX;
    h_dist = DBL_MAX;
    iter = 0;
    
    /* Computing HITS */
    printf("Computing PageRank...\n");
    while ((a_dist > TOL || h_dist > TOL) && iter < MAX_ITER)
    {
#ifdef DEBUG 
        if (iter % MOD_ITER == 0)
        {
#endif 
            printf("\riter %d", iter);
#ifdef DEBUG 
            printf("\n");
            printf("a: ");
            print_vec_f(a, no_nodes); 
            printf("h: ");
            print_vec_f(h, no_nodes); 
        }
#endif 
        
        /* a_new = Lt @ h, h_new = L @ a */
        for (ri = 0; ri < no_nodes; ++ri)
        {
            a_new[ri] = 0.;
            for (ci = row_ptr_t[ri]; ci < row_ptr_t[ri+1]; ++ci)
                a_new[ri] += h[col_ind_t[ci]];
            h_new[ri] = .0;
            for (ci = row_ptr[ri]; ci < row_ptr[ri+1]; ++ci)
                h_new[ri] += a[col_ind[ci]];
        }
        
        /* Normalization step */
        sum = 0.;
        for (i = 0; i < no_nodes; ++i) sum += a_new[i];
        for (i = 0; i < no_nodes; ++i) a_new[i] /= sum;
        sum = 0.;
        for (i = 0; i < no_nodes; ++i) sum += h_new[i];
        for (i = 0; i < no_nodes; ++i) h_new[i] /= sum;
        
        /* Computing distance between current and old a/h */
        a_dist = 0.;
        h_dist = 0.;
        for (i = 0; i < no_nodes; ++i)
        {
            a_dist += (a[i] - a_new[i]) * (a[i] - a_new[i]);
            h_dist += (h[i] - h_new[i]) * (h[i] - h_new[i]);
        }
        a_dist = sqrt(a_dist);
        h_dist = sqrt(h_dist);
        
        /* Copy new values in a/h */
        for (i = 0; i < no_nodes; ++i)
        {
            a[i] = a_new[i];
            h[i] = h_new[i];
        }
        
        ++iter;
    }
    printf("\riter %d\n", iter);
#ifdef DEBUG 
    printf("a: ");
    print_vec_f(a, no_nodes); 
    printf("h: ");
    print_vec_f(h, no_nodes); 
#endif 
    printf("Done.\n\n");
    
    printf("Proof of correctness:\n");
    sum = 0.; for (i = 0; i < no_nodes; ++i) sum += a[i];
    printf("sum(a) = %f\n", sum);
    sum = 0.; for (i = 0; i < no_nodes; ++i) sum += h[i];
    printf("sum(h) = %f\n", sum);
    
    /* un-mmapping data */
    munmap(row_ptr,   (no_nodes + 1) * sizeof(int));
    munmap(row_ptr_t, (no_nodes + 1) * sizeof(int));
    munmap(col_ind,   no_edges * sizeof(int));
    munmap(col_ind_t, no_edges * sizeof(int));
    
    /* Vectors of probability */
    free(a); free(a_new);
    free(h); free(h_new);
    
    exit(EXIT_SUCCESS);
}

/* Helper functions */

int write_data(char path[], void *data, size_t nmemb, size_t size)
{
    FILE *pdata;
    
    if ((pdata = fopen(path, "wb")) == NULL)
    {
        fprintf(stderr, " [ERROR] Cannot create file \"%s\"\n", path);
        return EXIT_FAILURE;
    }
    fwrite(data, size, nmemb, pdata);
    fclose(pdata);
    return EXIT_SUCCESS;
}

void delete_folder(char dir[])
{
    DIR *pf = opendir(dir);
	struct dirent *next_file;
	char fpath[PATH];
	
	while ((next_file = readdir(pf)) != NULL)
	{
		sprintf(fpath, "%s/%s", dir, next_file->d_name);
		remove(fpath);
	}
	closedir(pf);
	rmdir(dir);
}

void * mmap_data(char path[], size_t nmemb, size_t size)
{
    int fd;
    char mmap_p[MMAP];
    void * mp;
    sprintf(mmap_p, "./%s", path);
#ifdef DEBUG 
    printf("mmapping \"%s\"\n", mmap_p);
#endif 
    fd = open(mmap_p, O_RDONLY);
    mp = mmap(NULL, nmemb * size, PROT_READ, MAP_SHARED, fd, 0);
    if (mp == MAP_FAILED)
        mp = NULL;
    close(fd);
    return mp;
}

void print_vec_f(double *v, int n)
{
    int i;
    printf("[ ");
    for (i = 0; i < n; ++i)
        printf("%.3f ", v[i]);
    printf("]\n");
}

void print_vec_d(int *v, int n)
{
    int i;
    printf("[ ");
    for (i = 0; i < n; ++i)
        printf("%d ", v[i]);
    printf("]\n");
}

void double_merge(int *from, int *to, int lo, int mid, int hi)
{
    int *FL, *TL; 
    int *FR, *TR;
    int i, j;
    int NL, NR;
    
    NL = mid - lo;
    NR = hi - mid;
    
    FL = (int *) malloc(sizeof(int) * NL);
    TL = (int *) malloc(sizeof(int) * NL);
    FR = (int *) malloc(sizeof(int) * NR);
    TR = (int *) malloc(sizeof(int) * NR);
    
    for (i = 0; i < NL; ++i)
    {
        FL[i] = from[lo+i];
        TL[i] = to[lo+i];
    }
    
    for (j = 0; j < NR; ++j)
    {
        FR[j] = from[mid+j];
        TR[j] = to[mid+j];
    }
    
    i = 0; j = 0;
    while (i < NL && j < NR)
    {
        if (TL[i] < TR[j] || (TL[i] == TR[j] && FL[i] <= FR[j]))
        {
            from[lo+i+j] = FL[i];
            to[lo+i+j] = TL[i];
            ++i;
        }
        else 
        {
            from[lo+i+j] = FR[j];
            to[lo+i+j] = TR[j];
            ++j;
        }
    }
    
    while (i < NL)
    {
        from[lo+i+j] = FL[i];
        to[lo+i+j] = TL[i];
        ++i;
    }
    
    while (j < NR)
    {
        from[lo+i+j] = FR[j];
        to[lo+i+j] = TR[j];
        ++j;
    }
    
    free(FL); free(TL);
    free(FR); free(TR);
}

void double_merge_sort(int *from, int *to, int lo, int hi)
{
    if (hi - lo > 1)
    {
        int mid = (lo + hi) / 2;
        double_merge_sort(from, to, lo, mid);
        double_merge_sort(from, to, mid, hi);
        double_merge(from, to, lo, mid, hi);
    }
}

void sort_input_data(int *from, int *to, int n)
{
    double_merge_sort(from, to, 0, n);
}