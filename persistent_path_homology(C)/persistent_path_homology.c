#include <stdio.h>
#include <stdlib.h>
#include "persistent_path_homology.h"

/*  Auxiliary functions */

/*  Functions used to help the qsort function to sort the array representing the base  */
typedef struct {
    vertex_index index;
    double       allow_time_of_index;
} tuple;

double allow_time_auxiliary (double **network_weight, regular_path path, unsigned int path_dim) {
    /* Calculates the allow time of a regular path. It will be used to sort our basis */

    unsigned int j;
    double distance = 0.0;

    for (j = 0; j < path_dim; ++j)
        distance = distance < network_weight[path[j]][path[j+1]] ? network_weight[path[j]][path[j+1]] : distance;

    return distance;
} /*  Tested Ok */

int compareTuple (const void *p1, const void *p2) {
    if ( ((tuple*) p1)->allow_time_of_index <  ((tuple*) p2)->allow_time_of_index )       return -1;
    else if ( ((tuple*) p1)->allow_time_of_index ==  ((tuple*) p2)->allow_time_of_index ) return 0;
    else                                                                                  return 1;
} /*  Tested Ok */

/*  Auxiliary functions to help the algorithms involved with ComputePPH and BasisChange
 *  to work
 */
boolean are_these_regular_paths_the_same (regular_path path1, regular_path path2, unsigned int path_dim) {

    unsigned int i;
    for (i = 0; i <= path_dim; ++i)
        if ( path1[i] != path2[i] ) return FALSE;

    return TRUE;
}

boolean is_this_path_a_regular_path (regular_path path, unsigned int path_dim) {
    unsigned int i;

    for (i = 0; i < path_dim; ++i)
        if( path[i] == path[i+1] ) return FALSE;

    return TRUE;
}

boolean is_this_vector_zero (vector x, collection_of_basis *B, unsigned int path_dim) {
    unsigned int i;

    for (i = 0; i < (B->basis + path_dim)->dimension_of_the_vector_space_spanned_by_base; ++i)
        if (x[i] != 0) return FALSE;
    return TRUE;
}

/*  Main Functions */
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


    for (i = 0; i <= B->max_of_basis; ++i) {

        (B->basis + i)->marks = malloc ( (B->basis + i)->dimension_of_the_vector_space_spanned_by_base * sizeof (boolean) );

        for (j = 0; j < (B->basis + i)->dimension_of_the_vector_space_spanned_by_base; ++j){
            ((B->basis + i)->marks) [j] = NOT_MARKED;
        }
    }

    /*  Marking all regular paths of dimension 0 */
  
    for (j = 0; j < (B->basis)->dimension_of_the_vector_space_spanned_by_base; ++j)
        ((B->basis)->marks) [j] = MARKED;

} /*  Teste ok */


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
                                                                                 B,
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

    if (path_dim == 0) return 0.0;

    for (i = 0; i < base_dim; ++i) {
        if (path_vector[i] == TRUE) {
            temp_path = ((B->basis + path_dim)->base_matrix) [i];

            for (j = 0; j < path_dim; ++j) {
                distance = distance < network_weight[temp_path[j]] [temp_path[j+1]]
                           ? network_weight [temp_path[j]] [temp_path[j+1]] : distance;
            }
        }
    }
    return distance;
} /*  Tested Ok */


double entry_time (double **network_weight, collection_of_basis *B,
                   vector path_vector, unsigned int path_dim, unsigned int base_dim) {

    double distance;
    unsigned int i, j, k, l;
    regular_path boundary, temp;

    if (path_dim == 0) return 0.0;

    else if (path_dim == 1) return allow_time (network_weight, B, path_vector, path_dim, base_dim);

    else {
        distance = allow_time (network_weight, B, path_vector, path_dim, base_dim);

        /*Now we will have to calculate the boudary operator of the path_vector
         * then we will take its allow times */
        for (i = 0; i < base_dim; ++i) {

            if (path_vector[i] == TRUE) {
                temp = ((B->basis + path_dim)->base_matrix)[i];

                for (j = 0; j <= path_dim; ++j) {
                    boundary = malloc ((path_dim) * sizeof (vertex_index));
                    l = 0;
                    for (k = 0; k <= path_dim; ++k) {
                        if (k != j) {
                            boundary[l] = temp[k];
                            ++l;
                        }
                    }
                    if (is_this_path_a_regular_path (boundary, path_dim - 1) == FALSE) break;

                    for (l = 0; l < path_dim - 1; ++l) {
                        distance = distance < network_weight[boundary[l]] [boundary[l+1]]
                                ? network_weight [boundary[l]] [boundary[l+1]] : distance;
                    }
                    free (boundary);
                }
            }
        }
        return distance;
    }
}


boolean is_T_p_dim_i_vector_j_empty (T_p * Tp,
                                     unsigned int dim,
                                     unsigned int index) {
    return ((Tp->all_Tp + dim)->array_of_T_p_tuple + index)->is_empty; /*returns EMPTY or NOT_EMPTY*/
}


void fill_T_p_dim_i_vector_j (T_p *Tp,
                              unsigned int dim,
                              unsigned int index,
                              vector u,
                              double et) {

    ((Tp->all_Tp + dim)->array_of_T_p_tuple + index)->path_vector = u;
    ((Tp->all_Tp + dim)->array_of_T_p_tuple + index)->entry_time  = et;
    ((Tp->all_Tp + dim)->array_of_T_p_tuple + index)->is_empty    = NOT_EMPTY;

}


void sorting_the_basis_by_their_allow_times (collection_of_basis *B, double **network_weight) {

    unsigned int i, j, k, l, ii, iter_loop;
    tuple        *sort_this_array;
    regular_path copy1_regular_path, copy2_regular_path;
    boolean      *change_indexes;

    for (i = 1; i <= B->max_of_basis; ++i) {
        sort_this_array = malloc ((B->basis + i)->dimension_of_the_vector_space_spanned_by_base * sizeof (tuple) );

        for (j = 0; j < (B->basis + i)->dimension_of_the_vector_space_spanned_by_base; ++j ) {
            (sort_this_array + j)->index               = j;
            (sort_this_array + j)->allow_time_of_index = allow_time_auxiliary (network_weight, ((B->basis + i)->base_matrix)[j], i);
        }

        qsort (sort_this_array, (B->basis + i)->dimension_of_the_vector_space_spanned_by_base,
               sizeof (tuple), compareTuple);

        l = 0;
        j = 0;
        iter_loop = j;
        copy1_regular_path = malloc ((i + 1) * sizeof (vertex_index));
        copy2_regular_path = malloc ((i + 1) * sizeof (vertex_index));
        change_indexes     = malloc ((B->basis + i)->dimension_of_the_vector_space_spanned_by_base * sizeof (boolean));

        for (k = 0; k < i + 1; ++k) copy1_regular_path[k] = ((B->basis + i)->base_matrix)[0][k];

        for (k = 0; k < (B->basis + i)->dimension_of_the_vector_space_spanned_by_base; ++k ) change_indexes[k] = NOT_SORTED;

        while( l < (B->basis + i)->dimension_of_the_vector_space_spanned_by_base) {
            for (k = 0; k < (B->basis + i)->dimension_of_the_vector_space_spanned_by_base; ++k ) {

                if ((sort_this_array + k)->index == j) {
                    if (change_indexes[k] == NOT_SORTED) {
                        for (ii = 0; ii < i + 1; ++ii) copy2_regular_path[ii]               = ((B->basis + i)->base_matrix)[k][ii];
                        for (ii = 0; ii < i + 1; ++ii) ((B->basis + i)->base_matrix)[k][ii] = copy1_regular_path[ii];
                        for (ii = 0; ii < i + 1; ++ii) copy1_regular_path[ii]               = copy2_regular_path[ii];

                        change_indexes[k] = SORTED;
                        ++l;
                        j = k;
                        break;
                    }
                    else {
                        j = iter_loop + 1;
                        for (ii = 0; ii < i + 1; ++ii) copy1_regular_path[ii] = ((B->basis + i)->base_matrix)[j][ii];

                        iter_loop += 1;
                        break;
                    }
                }
            }
        }
        free (copy1_regular_path);
        free (copy2_regular_path);
        free (sort_this_array);
        free (change_indexes);
    } /*walking through the dimension of the regular paths*/
} /*  Tested Ok */


vector BasisChange (collection_of_basis *B, T_p *Tp, double **network_weight, vector path_vector, unsigned int path_dim,
                    double *return_et, unsigned int *return_max_index) {

    vector u;

    regular_path boundary_of_path_vector, temp;

    unsigned int i, j, l, k, max_index;

    double et = 0.0;

    boundary_of_path_vector = malloc (path_dim * sizeof (vertex_index));

    u = malloc ((B->basis + path_dim - 1)->dimension_of_the_vector_space_spanned_by_base * sizeof (boolean));

    max_index = 0;

    for (i = 0; i < (B->basis + path_dim)->dimension_of_the_vector_space_spanned_by_base; ++i) {

        if (path_vector[i] == TRUE) {
            /*apply the boundary operator*/
            temp = ((B->basis + path_dim)->base_matrix)[i];

            for (j = 0; j <= path_dim; ++j) {
                l = 0;
                for (k = 0; k <= path_dim; ++k) {
                    if (k != j) {
                        boundary_of_path_vector[l] = temp[k];
                        ++l;
                    }
                }
                if (is_this_path_a_regular_path (boundary_of_path_vector, path_dim - 1) == FALSE) break;

                for (k = 0; k < (B->basis + path_dim - 1)->dimension_of_the_vector_space_spanned_by_base; ++k) {
                    if (are_these_regular_paths_the_same ((B->basis + path_dim - 1)->base_matrix[k],
                                                          boundary_of_path_vector, path_dim - 1) == TRUE
                        && ((B->basis + path_dim - 1)->marks)[k] == TRUE) {

                        u[k] = (u[k] + 1) % 2;
                        if (u[k] != 0) max_index = max_index < k ? k : max_index;
                    }

                }
            }
        }
    } /*finished calcualting the border*/

    while (max_index > 0) {
        temp = ((B->basis + path_dim - 1)->base_matrix)[max_index];

        et = allow_time (network_weight, B, u, path_dim - 1, (B->basis + path_dim - 1)->dimension_of_the_vector_space_spanned_by_base)
            <
            allow_time_auxiliary (network_weight, temp, path_dim - 1 )
            ?
            allow_time_auxiliary (network_weight, temp, path_dim - 1 ) :
            allow_time (network_weight, B, u, path_dim - 1, (B->basis + path_dim - 1)->dimension_of_the_vector_space_spanned_by_base);

        if (((Tp->all_Tp + path_dim - 1)->array_of_T_p_tuple + max_index)->is_empty == EMPTY) break;

        for (k = 0; k < (B->basis + path_dim - 1)->dimension_of_the_vector_space_spanned_by_base; ++k) {
            u[k] = (u[k] + (((Tp->all_Tp + path_dim - 1)->array_of_T_p_tuple + max_index)->path_vector)[k] ) % 2;
        }

        /*now check again max_index*/
        max_index = 0;
        for (k = 0; k < (B->basis + path_dim - 1)->dimension_of_the_vector_space_spanned_by_base; ++k) {
            if (u[k] != 0) max_index = k;
        }
    }

    /*returning the values*/
    /*return_vector = u;*/
    *return_et = et;
    *return_max_index = max_index;
    return u;
}


Pers *ComputePPH(unsigned int pph_dim, double **network_weight, unsigned int network_set_size) {

    unsigned int        i, j, k, p, max_index;
    Pers                *PPH = malloc (sizeof (Pers));
    Pers_interval_p     *copy, *new;
    collection_of_basis B;
    T_p                 Tp;
    vector              v_j, u, v_i;
    double              et;



    PPH->pph = pph_dim;
    PPH->Dgm_all_dimensions = malloc ((pph_dim + 1) * sizeof (Pers_interval_p));

    for (i = 0; i <= pph_dim; ++i) {
        copy = malloc (sizeof (Pers_interval_p));
        copy->next = NULL;
        (PPH->Dgm_all_dimensions + i)->next = copy; /*initializing the stacks*/
    }

    /*Setting the environment*/
    Basis_of_the_vector_spaces_spanned_by_regular_paths (&B,  pph_dim, network_set_size);
    initialize_Marking_basis_vectors (&B);
    sorting_the_basis_by_their_allow_times (&B, network_weight);
    generating_T_p (&Tp, &B, network_weight);

    /*Now lets start the algorithm*/
    for (p = 0; p <= pph_dim; ++p) {

        u   = malloc (((B.basis) + p)->dimension_of_the_vector_space_spanned_by_base * sizeof (boolean));
        v_j = malloc (((B.basis) + p + 1)->dimension_of_the_vector_space_spanned_by_base * sizeof (boolean));
        v_i = malloc (((B.basis) + p)->dimension_of_the_vector_space_spanned_by_base * sizeof (boolean));

        for (j = 0; j < ((B.basis) + p + 1)->dimension_of_the_vector_space_spanned_by_base; ++j) {

            for (k = 0; k < ((B.basis) + p + 1)->dimension_of_the_vector_space_spanned_by_base; ++k) {
                if (k == j) v_j[k] = TRUE;
                else v_j[k] = FALSE;
            }

            u = BasisChange (&B, &Tp, network_weight, v_j, p + 1, &et, &max_index);

            if (is_this_vector_zero (u, &B, p) == TRUE)
                marking_vector_basis (&B, p + 1, j);
            else {
                fill_T_p_dim_i_vector_j (&Tp, p, max_index, u, et);

                /*vector v_i*/
                for (k = 0; k < ((B.basis) + p)->dimension_of_the_vector_space_spanned_by_base; ++k) {
                    if (k == max_index) v_i[k] = TRUE;
                    else v_i[k] = FALSE;
                }
                ((PPH->Dgm_all_dimensions + p)->PPH_interval_dim_p)[0] = entry_time (network_weight,
                                                                                     &B,
                                                                                     v_i,
                                                                                     p,
                                                                                     ((B.basis) + p)->dimension_of_the_vector_space_spanned_by_base);
                ((PPH->Dgm_all_dimensions + p)->PPH_interval_dim_p)[1] = et;

                new       = malloc (sizeof (Pers_interval_p));
                new->next = (PPH->Dgm_all_dimensions + p)->next;
                (PPH->Dgm_all_dimensions + p)->next = new;
            }
        }
        for (j = 0; j < ((B.basis) + p)->dimension_of_the_vector_space_spanned_by_base; ++j) {

            if (is_T_p_dim_i_vector_j_empty (&Tp, p, j) == EMPTY &&
                (((B.basis) + p)->marks)[j] == MARKED ) {

                /*vector v_i*/
                for (k = 0; k < ((B.basis) + p)->dimension_of_the_vector_space_spanned_by_base; ++k) {
                    if (k == max_index) v_i[k] = TRUE;
                    else v_i[k] = FALSE;
                }

                ((PPH->Dgm_all_dimensions + p)->PPH_interval_dim_p)[0] = entry_time (network_weight,
                                                                                     &B,
                                                                                     v_i,
                                                                                     p,
                                                                                     ((B.basis) + p)->dimension_of_the_vector_space_spanned_by_base);
                ((PPH->Dgm_all_dimensions + p)->PPH_interval_dim_p)[1] = -1.0;

                new       = malloc (sizeof (Pers_interval_p));
                new->next = (PPH->Dgm_all_dimensions + p)->next;
                (PPH->Dgm_all_dimensions + p)->next = new;

            }
        }

        free (u);
        free (v_i);
        free (v_j);
    }

    return PPH;
}
