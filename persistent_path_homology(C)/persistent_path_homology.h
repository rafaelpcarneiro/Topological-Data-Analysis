#ifndef __PERSISTENT_PATH_HOMOLOGY_H_
#define __PERSISTENT_PATH_HOMOLOGY_H_


/*================================================================
 * BASIC DATA TYPES
 * ----------------
 *
 * Explanation of each data type.
 *
 *================================================================
 * boolean      := a type of data that stores TRUE or FALSE
 *
 * vertex_index := indexes of the vertices of ours graphs
 *
 * regular_path := an array where each element is a vertex_index
 *
 * vector       := an array of booleans representing an
 *                 element of a Z/2Z-vector space of dimension p.
 *                 For instance, fixed a base B of such vector space,
 *                 then
 *                 vector w = {0, 1, 0, 0, 1};
 *                 means that w = a_1 + a_4, with
 *                 a_1, a_4 elements of the base B
 *================================================================*/

#define TRUE 1
#define FALSE 0

typedef char boolean;

typedef unsigned int vertex_index;

typedef vertex_index *regular_path;

typedef boolean *vector;


/*================================================================
 * DEALING WITH THE BASIS
 * ----------------------
 *
 * Explanation of each structs data type.
 *
 *================================================================
 * boolean      := a type of data that stores TRUE or FALSE
 *
 * vertex_index := indexes of the vertices of ours graphs
 *
 * regular_path := an array where each element is a vertex_index
 *
 * vector       := an array of booleans representing an
 *                 element of a Z/2Z-vector space of dimension p.
 *                 For instance, fixed a base B of such vector space,
 *                 then
 *                 vector w = {0, 1, 0, 0, 1};
 *                 means that w = a_1 + a_4, with
 *                 a_1, a_4 elements of the base B
 *================================================================*/

typedef struct{
    vertex_index  ***basis_matrix;
    unsigned int  number_of_vector_spaces;
    vertex_index  *dim_of_each_vector_spaces;
} base_dim_p;

typedef struct{
    vector       path_vector; /* an element of a vector space spanned
                               * by regular paths of dimension dim */
    unsigned int entry_time;
    boolean      is_empty;
    unsigned     int dim;
} T_p_tuple;

typedef struct{
}

#endif // __PERSISTENT_PATH_HOMOLOGY_H_
