#include <stdio.h>
#include <stdlib.h>
#include "persistent_path_homology.h"

int main() {

    unsigned int i, j, k;
    collection_of_basis B;
    base teste;

    Basis_of_the_vector_spaces_spanned_by_regular_paths (&B, 1, 3);

    for (i = 0; i <= 2; ++i) {
        teste = (B.basis)[i];
        for (j = 0; j < teste.dimension_of_the_vector_space_spanned_by_base; ++j) {
            printf ("Dimension = %u,    ", teste.dimension_of_the_regular_path );

            printf ("[");
            for (k = 0; k <= teste.dimension_of_the_regular_path; ++k) {
                printf ("%u ", (teste.base_matrix)[j][k] );
            }
            printf ("]");
            printf("\n");
        }
    }

    return 0;
}
