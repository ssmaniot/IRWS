#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
    unsigned *row_ptr, *col_ptr;
    unsigned from, to;
    double *val;
    unsigned n, e;
    FILE *pf;
    
    if (argc != 2)
    {
        printf("no file was given!\n");
        exit(EXIT_FAILURE);
    }
    
    pf = fopen(argv[1], "r");
    if (pf == NULL)
    {
        printf("cannot open file `%s'\n", argv[1]);
        exit(EXIT_FAILURE);
    }
    
    /* Skip the first two rows */
    fscanf(pf, "%*[^\n]\n");
    fscanf(pf, "%*[^\n]\n");
    /* Reading graph order and size */
    fscanf(pf, "# Nodes: %u Edges: %u", &n, &e);
    /* Skip one more line */
    fscanf(pf, "%*[^\n]\n");
    
    while (feof(pf))
    {
        fscanf(pf, "%u\t%u", &from, &to);
    }
    
    fclose(pf);
    
    return EXIT_SUCCESS;
}