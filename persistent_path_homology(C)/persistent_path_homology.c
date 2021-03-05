#include <stdio.h>
#include <stdlib.h>
#include "persistent_path_homology.h"

void generating_all_regular_paths_dim_p (collection_of_basis *B,
                                         unsigned int dim_p,
                                         unsigned int network_set_size){

    base *B_dim_p           = (B -> basis) + dim_p;
    base *B_dim_p_minus_one = (B -> basis) + (dim_p - 1);

    regular_path temp_dim_p, temp_dim_p_minus_one;

    unsigned int i, j, k, l;


    B_dim_p -> dimension_of_the_regular_path                 = dim_p;
    B_dim_p -> dimension_of_the_vector_space_spanned_by_base =  (B_dim_p_minus_one -> dimension_of_the_vector_space_spanned_by_base) *
        (network_set_size - 1);

    (B_dim_p -> base_matrix) = malloc ( (B_dim_p -> dimension_of_the_vector_space_spanned_by_base) * sizeof (regular_path) );

    l = 0;
    for (i = 0; i < B_dim_p_minus_one -> dimension_of_the_vector_space_spanned_by_base; ++i) {

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

void Basis_of_the_vector_spaces_spanned_by_regular_paths (collection_of_basis *B,
                                                          unsigned int pph_dim,
                                                          unsigned int network_set_size) {

    unsigned int i;

    B -> basis        = malloc( (pph_dim + 2) * sizeof(base) );
    B -> max_of_basis = pph_dim + 1;

    (B -> basis) -> base_matrix                                   = malloc ( network_set_size * sizeof (regular_path) );
    (B -> basis) -> dimension_of_the_regular_path                 = 0;
    (B -> basis) -> dimension_of_the_vector_space_spanned_by_base = network_set_size;


    for (i = 0; i < network_set_size; ++i) {

        ((B -> basis) -> base_matrix)[i]    = malloc ( sizeof (vertex_index) );
        ((B -> basis) -> base_matrix)[i][0] = i;

    }

    for (i = 1; i <= pph_dim + 1; ++i) {
        generating_all_regular_paths_dim_p (B, i, network_set_size);
    }
} /*  Tested Ok */

void initialize_Marking_basis_vectors (collection_of_basis *B) {

    unsigned int i, j;

    for (i = 0; i <= B->max_of_basis + 1; ++i) {
        (B->basis + i)->marks = malloc ( (B->basis + i)->dimension_of_the_vector_space_spanned_by_base * sizeof (boolean) );

        for (j = 0; j < (B->basis + i)->dimension_of_the_vector_space_spanned_by_base; ++j)
            ((B->basis + i)->marks) [j] = NOT_MARKED;
    }

    /*  Marking all regular paths of dimension 0 */
    for (j = 0; j < (B->basis)->dimension_of_the_vector_space_spanned_by_base; ++j)
        ((B->basis)->marks) [j] = MARKED;
}

void marking_vector_basis (collection_of_basis *B,
                           unsigned int vector_path_dim,
                           unsigned int vector_index){

    ((B->basis + vector_path_dim)->marks) [vector_index] = MARKED;
}

void generating_T_p ( T_p *Tp, collection_of_basis *B, double **network_weight) {

    unsigned int i, j, k;

    Tp->all_Tp    = malloc (B->max_of_basis * sizeof (T_p_tuple_collection));
    Tp->max_of_Tp = B->max_of_basis;

    for (i = 0; i <= B->max_of_basis; ++i) {

        (Tp->all_Tp + i)->array_of_T_p_tuple = malloc ((B->basis + i)->dimension_of_the_vector_space_spanned_by_base * sizeof (T_p_tuple));
        (Tp->all_Tp +i)->size                = (B->basis + i)->dimension_of_the_vector_space_spanned_by_base;

        for (j = 0; j < (Tp->all_Tp +i)->size; ++j) {
            ((Tp->all_Tp + i)->array_of_T_p_tuple + j)->path_vector = malloc ((Tp->all_Tp +i)->size * sizeof (boolean));

            for (k = 0; k < (Tp->all_Tp +i)->size; ++k) {
                if ( j == k ) (((Tp->all_Tp + i)->array_of_T_p_tuple + j)->path_vector) [k] = TRUE;
                else          (((Tp->all_Tp + i)->array_of_T_p_tuple + j)->path_vector) [k] = FALSE;
            }

            ((Tp->all_Tp + i)->array_of_T_p_tuple + j)->entry_time = entry_time (network_weight,
                                                                                 ((Tp->all_Tp + i)->array_of_T_p_tuple + j)->path_vector,
                                                                                 i, (Tp->all_Tp +i)->size);
            ((Tp->all_Tp + i)->array_of_T_p_tuple + j)->is_empty   = EMPTY;
        }
    }
}


double allow_time (double **network_weight, collection_of_basis *B,
                   vector path_vector, unsigned int path_dim, unsigned int base_dim) {

    unsigned int i, j;
    double distance = 0.0;
    regular_path temp_path;

    if (path_dim == 0) return 0;

    for (i = 0; i <= base_dim; ++i) {
        if (path_vector[i] == TRUE) {
            temp_path = ((B->basis + path_dim)->base_matrix) [i];

            for (j = 0; j < path_dim; ++j) {
                distance = distance < network_weight[temp_path[j]] [temp_path[j+1]]
                           ? network_weight [temp_path[j]] [temp_path[j+1]] : distance;
            }
        }
    }
    return distance;
}


