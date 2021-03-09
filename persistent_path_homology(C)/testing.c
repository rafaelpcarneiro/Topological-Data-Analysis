#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include "persistent_path_homology.h"

int main() {

    unsigned int i, j, k;
    collection_of_basis B;
    base teste;

    unsigned int network_set_size = 3;
    unsigned int pph = 1;

    double **network_weight = malloc (network_set_size * sizeof (double*));

    T_p Tp;
    T_p_tuple_collection Tp_dim_i;

    gsl_rng *random_number = gsl_rng_alloc (gsl_rng_mt19937);
    gsl_rng_set (random_number, 10);

    printf ("Network_weight:\n");
    for (i = 0; i < network_set_size; ++i) {
        network_weight[i] = malloc (network_set_size * sizeof (double));
        for (j = 0; j < network_set_size; ++j) {
            if (i == j) network_weight[i][j] = 0.0;
            else network_weight[i][j] = gsl_rng_uniform (random_number);
            printf ("%6.2f ", network_weight[i][j]);
        }
        printf ("\n");
    }

    Basis_of_the_vector_spaces_spanned_by_regular_paths (&B, pph, network_set_size);

    /*  printing the basis */
    for (i = 0; i <= pph + 1; ++i) {
        teste = (B.basis)[i];
        for (j = 0; j < teste.dimension_of_the_vector_space_spanned_by_base; ++j) {
            printf ("Dimension = %u, size = %u   ", teste.dimension_of_the_regular_path,
                    teste.dimension_of_the_vector_space_spanned_by_base);

            printf ("[");
            for (k = 0; k <= teste.dimension_of_the_regular_path; ++k) {
                printf ("%u ", (teste.base_matrix)[j][k] );
            }
            printf ("]");
            printf("\n\n");
        }
    }


    initialize_Marking_basis_vectors (&B);

    /*  printing the marks */
   
    for (i = 0; i <= pph + 1; ++i) {
        teste = (B.basis)[i];
        printf ("Dimension = %u\n", teste.dimension_of_the_regular_path );
        for (j = 0; j < teste.dimension_of_the_vector_space_spanned_by_base; ++j) {
            printf ("%u  ", (teste.marks)[j] );
        }
            printf("\n\n");
    }



    sorting_the_basis_by_their_allow_times (&B, network_weight);

    /*  printing the basis now */
    printf ("\nChecking\n");
    for (i = 0; i <= pph + 1; ++i) {
        teste = (B.basis)[i];
        for (j = 0; j < teste.dimension_of_the_vector_space_spanned_by_base; ++j) {
            printf ("Dimension = %u,    ", teste.dimension_of_the_regular_path );

            printf ("[");
            for (k = 0; k <= teste.dimension_of_the_regular_path; ++k) {
                printf ("%u ", (teste.base_matrix)[j][k] );
            }
            printf ("], allow time = %.2f", allow_time_auxiliary (network_weight, teste.base_matrix[j], i) );
            printf("\n\n");
        }
    }

    generating_T_p (&Tp, &B, network_weight);


    return 0;
}
