#ifndef __DEFINITIONS_H_
#define __DEFINITIONS_H_

/*================================================================
 * BASIC DATA TYPES
 * ---------------------------------------------------------------
 *
 *================================================================
 * boolean      := a type of data that stores TRUE or FALSE
 *
 * vertex_index := indexes of the vertices of ours graphs
 *
 * base_index   := indexes of the elements of a given oredered
 *                 base
 *
 * dim_path     := a type making clear that the integer is relative
 *                 to the dimension of a regular path
 *
 * dim_base     := a type making clear that the integer is relative
 *                 to the dimension of a vector space spanned by
 *                 a base
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

#define TRUE       1
#define FALSE      0

#define MARKED     1
#define NOT_MARKED 0

#define EMPTY      1
#define NOT_EMPTY  0

#define SORTED     1
#define NOT_SORTED 0

typedef char         boolean;

typedef unsigned int vertex_index;

typedef unsigned int base_index;

typedef unsigned int dim_path;

typedef unsigned int dim_vector_space;

typedef vertex_index *regular_path;

typedef boolean      *vector;


#endif // __DEFINITIONS_H_
