#ifndef __PERSISTENT_PATH_HOMOLOGY_H_
#define __PERSISTENT_PATH_HOMOLOGY_H_


/*================================================================
 * BASIC DATA TYPES
 * ---------------------------------------------------------------
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

#define TRUE       1
#define FALSE      0

#define MARKED     1
#define NOT_MARKED 0

#define EMPTY      1
#define NOT_EMPTY  0

#define SORTED     1
#define NOT_SORTED 0

typedef char boolean;

typedef unsigned int vertex_index;

typedef vertex_index *regular_path;

typedef boolean *vector;


/*================================================================
 * SETTING THE DATA TYPE TO DEAL WITH THE VECTOR SPACES,
 * REPRESENTED BY THEIR BASIS
 * ---------------------------------------------------------------
 *
 * Explanation of each of the data types.
 *
 *================================================================
 * base := a struct storing all info relative to a base of a vector
 *   |      space spanned by regular paths of dimension p
 *   |
 *   --> base_matrix := an array of regular regular paths with
 *   |                   the same dimension;
 *   |
 *   --> dimension_of_the_regular_path := the dimension of the
 *   |                   regular paths generating the vector space
 *   |
 *   --> dimension_of_the_vector_space_spanned_by_base := the number
 *   |                   of regular paths stored at basis_matrix
 *   |
 *   --> marks        := an array which informs wheter or not an
 *                      element of base is marked.
 *
 * collection_of_base := a struct containing all base structures.
 *   |
 *   |
 *   --> basis := an array of bases
 *   |
 *   |
 *   --> max_of_basis := an int saying how many elements basis
 *                has. That is, how many base structures we have
 *                created
 *================================================================*/

typedef struct{
    regular_path  *base_matrix;
    unsigned int  dimension_of_the_regular_path;
    unsigned int  dimension_of_the_vector_space_spanned_by_base;
    boolean       *marks;

} base;

typedef struct{
    base         *basis;
    unsigned int max_of_basis;

} collection_of_basis;

/*================================================================
 * SETTING THE DATA TYPE TO DEAL WITH THE T_P STRUCTURE,
 * ---------------------------------------------------------------
 *
 * Explanation of each of these data types.
 *
 *================================================================
 * T_P_tuple := a struct storing the tuple (path_vector,
 *   |           entry_time, is_empty, dim).
 *   |
 *   --> path_vector:= an array of booleans representing a vactor
 *   |                  in a given base.
 *   |
 *   --> entry_time := the moment when a regular path w is feasible
 *   |              at the filtration.
 *   |
 *   --> is_empty   := checking if this tuple is empty or not;
 *   |
 *   |
 *   --> dim        := the dimension of the regular paths that span
 *                     the vector spach which path_vector is element
 *
 * T_p_tuple_collection := a struct with all T_p_tuple of same dim.
 *   |
 *   |
 *   --> array_of_T_p_tuple := an array with all tuples  sharing
 *   |                         the condition that path_vectors have
 *   |                         same dimension
 *   |
 *   --> size               := the size of the array above
 *
 * T_p := a struct which finnaly assembly all T_p ranging from all
 *   |    dimensions
 *   |
 *   --> all_Tp := an array containing all T_p_tuple_collection
 *   |
 *   |
 *   --> dim    := maximum number of T_p_tuple_collection to
 *                 store;
 *================================================================*/
typedef struct{
    vector       path_vector;
    double       entry_time;
    boolean      is_empty;

} T_p_tuple;

typedef struct{
    T_p_tuple    *array_of_T_p_tuple;
    unsigned int size;

} T_p_tuple_collection;

typedef struct{
    T_p_tuple_collection *all_Tp;
    unsigned int max_of_Tp;

} T_p;

/*================================================================
 * SETTING THE DATA TYPE TO STORE THE PERSISTENT PATH HOMOLOYG
 * DIAGRAMS
 * ---------------------------------------------------------------
 *
 *================================================================
 *================================================================*/

/*  a stack is introduced down bellow  */
typedef struct persistent_interval{
    double                     PPH_interval_dim_p[2];
    struct persistent_interval *next;
} Pers_interval_p;

typedef struct {
    Pers_interval_p *Dgm_all_dimensions;
    unsigned int pph;
} Pers;




/*================================================================
 * ALL FUNCTIONS OPERATING ON OUR DATA
 * ---------------------------------------------------------------
 *
 * Explanation of each of these functions,
 *
 *================================================================
 * generating_all_regular_paths_dim_p( ... )
 *  |
 *  |
 *  --> Given a base of dimension dim - 1, representing
 *      all regular paths of dimension dim - 1, this function
 *      will create a base of dimension dim
 *
 * Basis_of_the_vector_spaces_spanned_by_regular_paths( ... );
 *  |
 *  |
 *  --> Here we will be storing the regular paths of
 *      dimension 0, 1, ..., pph_dim + 1 into their repective basis
 *
 * allow_time (...)
 *  |
 *  |
 *  --> Given the vector path_vector, an element of the vector
 *      space spanned by the regular paths of dimension path_dim
 *      That is
 *          path_vector := sum_i alpha_i * v_i
 *      with alpha_i in Z/2Z and v_i base vectors.
 *
 *      Then this function will calculate the time of birth of
 *      such vector. In another words, this method
 *      returns the value
 *          min { k in A_k; path_vector in A_k, k >= 0 },
 *      where A_k is the set of the allowed regular paths of
 *      dimension k.
 *
 * entry_time (...)
 *  |
 *  |
 *  --> The same as allow_time but it returns
 *          min { k in A_k; path_vector in A_k and
 *                d(path_vector) in A_k, k >= 0 },
 *
 * sorting_the_basis_by_their_allow_times (...)
 *  |
 *  |
 *  --> Sorting the structures T_p and basis in agreement
 *      with the allow times.
 *
 * initialize_Marking_basis_vectors (...);
 *
 * marking_vector_basis (...)
 *
 * generating_T_p
 *
 * is_T_p_dim_i_vector_j_empty (...);
 *
 * fill_T_p_dim_i_vector_j (...);
 *
 * BasisChange (...);
 *  |
 *  |
 *  --> Implementing the function BasisChange as in the paper
 *
 * ComputePPH (...);
 *================================================================*/

void generating_all_regular_paths_dim_p (collection_of_basis *B,
                                         unsigned int dim_p,
                                         unsigned int network_set_size);


void Basis_of_the_vector_spaces_spanned_by_regular_paths (collection_of_basis *B,
                                                          unsigned int pph_dim,
                                                          unsigned int network_set_size);


void initialize_Marking_basis_vectors (collection_of_basis *B);


void sorting_the_basis_by_their_allow_times (collection_of_basis *B, double **network_weight);


double allow_time (double **network_weight, collection_of_basis *B,
                   vector path_vector, unsigned int path_dim, unsigned int base_dim);


double entry_time (double **network_weight, collection_of_basis *B,
                   vector path_vector, unsigned int path_dim, unsigned int base_dim);


void marking_vector_basis (collection_of_basis *B,
                           unsigned int vector_path_dim,
                           unsigned int vector_index);


void generating_T_p (T_p *Tp, collection_of_basis *B, double **network_weight);


boolean is_T_p_dim_i_vector_j_empty (T_p *Tp,
                                     unsigned int dim,
                                     unsigned int index);


void fill_T_p_dim_i_vector_j (T_p *Tp,
                              unsigned int dim,
                              unsigned int index,
                              vector u,
                              double et);


vector BasisChange (collection_of_basis *B, T_p *Tp, double **network_weight, vector path_vector, unsigned int path_dim,
                    double *return_et, unsigned int *return_max_index);


Pers *ComputePPH(unsigned int pph_dim,
                 double **network_weight,
                 unsigned int network_set_size);


double allow_time_auxiliary (double **network_weight, regular_path path, unsigned int path_dim);

#endif /* __PERSISTENT_PATH_HOMOLOGY_H_ */
