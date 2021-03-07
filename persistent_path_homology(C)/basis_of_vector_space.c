#include <stdlib.h>
#include "basis_of_vector_space.h"

dim_vector_space get_dimVS_of_ith_base (collection_of_basis *B, dim_path dim_p) {
    base *B_dim_p = (B->basis) + dim_p;
    return (B_dim_p->dimension_of_the_vs_spanned_by_base);
}

void set_dim_path_of_ith_base (collection_of_basis *B, dim_path dim_p) {
    base *B_dim_p                          = (B->basis) + dim_p;
    B_dim_p->dimension_of_the_regular_path =  dim_p;
}

void set_dimVS_of_ith_base (collection_of_basis *B, dim_path dim_p, dim_vector_space dimVS) {
    base *B_dim_p                                = (B->basis) + dim_p;
    B_dim_p->dimension_of_the_vs_spanned_by_base =  dimVS;
}

collection_of_basis *alloc_all_basis (unsigned int number_of_bases_to_allocate,
                                      unsigned int network_set_size) {

    collection_of_basis       *B = malloc (sizeof (collection_of_basis));
    base_index                i;
    tuple_regular_path_double ith_tuple;

    B->basis        = malloc( (number_of_bases_to_allocate + 1) * sizeof(base) );
    B->max_of_basis = number_of_bases_to_allocate;

    (B->basis)->base_matrix = malloc ( network_set_size * sizeof (tuple_regular_path_double) );

    (B->basis)->dimension_of_the_regular_path       = 0;
    (B->basis)->dimension_of_the_vs_spanned_by_base = network_set_size;

    for (i = 0; i < network_set_size; ++i) {
        /*Base referent to regular paths of dimension 0*/
        ith_tuple             = ((B->basis)->base_matrix)[i];
        ith_tuple.ith_base    = malloc ( sizeof (vertex_index) );
        ith_tuple.ith_base[0] = i;
        ith_tuple.allow_time  = 0.0;
    }

    for (i = 1; i <= number_of_bases_to_allocate; ++i) {
        /*Base referent to regular paths of dimension > 0*/
        generating_all_regular_paths_dim_p (B, i, network_set_size);
    }

    return B;
} /*  Tested Ok */

void generating_all_regular_paths_dim_p (collection_of_basis *B,
                                         dim_path dim_p,
                                         unsigned int network_set_size){

    base *B_dim_p           = (B->basis) + dim_p;
    base *B_dim_p_minus_one = (B->basis) + (dim_p - 1);

    tuple_regular_path_double temp_dim_p, temp_dim_p_minus_one;

    unsigned int i, j, k, l;

    set_dim_path_of_ith_base (B, dim_p);
    set_dimVS_of_ith_base (B, dim_p, get_dimVS_of_ith_base (B, dim_p - 1) * (network_set_size - 1));

    B_dim_p->base_matrix = malloc ( get_dimVS_of_ith_base (B, dim_p) * sizeof (tuple_regular_path_double) );

    l = 0;
    for (i = 0; i < get_dimVS_of_ith_base (B, dim_p - 1); ++i) {
        /*  continue here */
        temp_dim_p_minus_one = (B_dim_p_minus_one -> base_matrix) [i];

        for (j = 0; j < network_set_size; ++j) {
            temp_dim_p = malloc( (dim_p + 1) * sizeof(vertex_index) );

            if ( temp_dim_p_minus_one [dim_p - 1] != j ) {
                for (k = 0; k <= dim_p - 1; ++k) temp_dim_p[k] = temp_dim_p_minus_one[k];
                temp_dim_p[k] = j;
                (B_dim_p -> base_matrix) [l] = temp_dim_p;
                ++l;
            }
        }
    }
} /*  Tested Ok */


void initialize_Marking_basis_vectors (collection_of_basis *B);


void sorting_the_basis_by_their_allow_times (collection_of_basis *B, double **network_weight);

void marking_vector_basis (collection_of_basis *B,
                           unsigned int vector_path_dim,
                           unsigned int vector_index);

double allow_time_auxiliary (double **network_weight, regular_path path, unsigned int path_dim);

int compareTuple (const void *p1, const void *p2);
