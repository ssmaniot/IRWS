#include "linalg.h"
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
    csr_matrix At;
    float m[] = { 0.f, 1.f, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f };

    At = new_csr_matrix__(m, 3, 3);

    print_csr_matrix(At);

    delete_csr_matrix(&At);
    
    return EXIT_SUCCESS;
}