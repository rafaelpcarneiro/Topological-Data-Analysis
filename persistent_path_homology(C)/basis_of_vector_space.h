#ifndef __BASIS_OF_VECTOR_SPACE_H_
#define __BASIS_OF_VECTOR_SPACE_H_

#include "definitions.h"

/*================================================================
 * SETTING THE DATA TYPE TO DEAL WITH THE VECTOR SPACES,
 * REPRESENTED BY THEIR BASIS
 * ---------------------------------------------------------------
 *
 * Explanation of each of the data types.
 *
 *================================================================
 * tuple_regular_path_double := a tuple wich coordinates given by
 *   |                          the variables below.
 *   |
 *   --> ith_base   :=
 *   |       a regular regular path of the ith element of
 *   |       the base
 *   |
 *   --> allow_time :=
 *           the allow time of the ith element of the
 *           base
 *
 * base := a struct storing all info relative to a base of a vector
 *   |      space spanned by regular paths of dimension p
 *   |
 *   |
 *   --> base_matrix :=
 *   |       An array storing all regular paths of same dimension
 *   |       whose function is to span the vector space of all
 *   |       regular paths of same dimension
 *   |
 *   |
 *   --> dimension_of_the_regular_path :=
 *   |       the dimension of the regular paths generating the
 *   |       vector space
 *   |
 *   --> dimension_of_the_vector_space_spanned_by_base :=
 *   |       the number of regular paths stored at basis_matrix
 *   |
 *   --> marks :=
 *           an array which informs wheter or not an element of
 *           base is marked.
 *
 * collection_of_base := a struct containing all base structures.
 *   |
 *   |
 *   --> basis :=
 *   |        an array of bases
 *   |
 *   --> max_of_basis :=
 *           an int saying how many elements basis has. That is,
 *           how many base structures we have created
 *================================================================*/

typedef struct {
    regular_path ith_base;
    double       allow_time;

} tuple_regular_path_double;

/*  Here 'vs' stands for vector space */
typedef struct{
    tuple_regular_path_double  *base_matrix;
    unsigned int               dimension_of_the_regular_path;
    unsigned int               dimension_of_the_vs_spanned_by_base;
    boolean                    *marks;

} base;

typedef struct{
    base         *basis;
    unsigned int max_of_basis;

} collection_of_basis;

/*  FUNCTIONS OPERATING ON THESE STRUCTS  */
collection_of_basis *alloc_all_basis (dim_vector_space, unsigned int);

void generating_all_regular_paths_dim_p (collection_of_basis*,
                                         dim_path,
                                         unsigned int);


void Basis_of_the_vector_spaces_spanned_by_regular_paths (collection_of_basis*,
                                                          unsigned int,
                                                          unsigned int);


void initialize_Marking_basis_vectors (collection_of_basis*);


void sorting_the_basis_by_their_allow_times (collection_of_basis*, double**);

void marking_vector_basis (collection_of_basis*, unsigned int, unsigned int);

double allow_time_auxiliary (double**, regular_path, unsigned int);

int compareTuple (const void*, const void*);

/*  setters and getters */
dim_vector_space get_dimVS_of_ith_base (collection_of_basis*, dim_path);

void set_dim_path_of_ith_base (collection_of_basis*, dim_path);

void set_dimVS_of_ith_base (collection_of_basis*, dim_path, dim_vector_space);
#endif // __BASIS_OF_VECTOR_SPACE_H_
