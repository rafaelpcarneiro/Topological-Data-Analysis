#!/usr/bin/env python3

import numpy as np
from math import *

#############################################
######## Auxiliary Functions
#############################################

def generating_all_regular_paths_dim_p( B_i, network_set_size ):
    """
    Given a numpy array B_i of shape (a,b), representing
    all regular paths of dimension b, this function
    will return a numpy array B_{i+1} of shape (a * (network_set_size - 1), b+1),
    representing all regular paths of dimension b+1.

    This auxiliary function will be used at the method:
       PPH.Basis_of_the_vector_spaces_spanned_by_regular_paths().
    """

    # Notation ->  B_next := B_{i+1}
    size_B_next = B_i.shape[0] * (network_set_size - 1 )
    B_next = np.zeros( (size_B_next, B_i.shape[1] + 1) )

    for i in range( B_i.shape[0] ):
        for j in range( network_set_size ):
            if B_i[ i, -1 ] != j:
                B_next[i, 0:B_i.shape[0] ] = B_i[i,:]
                B_next[i, -1 ] = j

    return B_next






#############################################
######## Main Class: PPH
#############################################
class PPH:
    """
    This class will calculate the persitent path homology
    of a network (network_set, network_weight) which satisfies
      (i)  network_set is a finite set represented by a
           numpy array;

      (ii) network_weight is a function definied mathematically
           as:
           network_weight: network_set x network_set -> R_+

           The machine will store this function as a numpy array
           of dimension (network_set, network_set)


    After storing the network the class PPH will calculate its
    persistent path homology of a given dimension d. The method
    responsible for such calculation is: diagram

    All the algorithm deployed here is based at the theory developed
    at the following paper:

       + PERSISTENT PATH HOMOLOGY OF DIRECTED NETWORKS, from the
         authors Samir Chowdhury and Facundo Mémoli

    This paper can be easily found at the website:
       -> https://arxiv.org/abs/1701.00565

    The same authors have a paper with an algorithm implemented
    to calculate the persistent path homology. The algorithm
    resembles the one to calculate the persitent diagram when
    considering the field Z/2Z over the vector spaces from
    the homology groups.

    The paper with the algorithm can be found here:
       -> https://epubs.siam.org/doi/10.1137/1.9781611975031.75
    exectly at the end of the file.
    """

    def __init__(self, network_set, network_weight, pph_dim):
        """
        atributes:
          + network_set: a numpy array storing the sets of
                         the network

          + network_weight: a nupy array storing the weight
                            function which do not need to be symmetric

          + network_set_size: integer representing the size of our
                              data

          + pph_dim: dimension of the persistent path diagram

          + basis: a list looking like

                 basis := [ B0, B_1, B_2, ..., B_(self.dim + 1) ].

            where each element is a numpy array:
                B0: stores regular paths of dimension 0;
                B1: stores regular paths of dimension 1;
                B2: stores regular paths of dimension 2;
                                 .
                                 .
                                 .
          + basis_dim: a list containing as elements the values
                       of the dimensions of each vector space
                       spaned by the basis B0, B1, B2, above.
                       In other words:
                       basis_dim := [ basis[0].shape[0],
                                      basis[1].shape[0], ... ]
        """

        self.pph_dim   = pph_dim
        self.basis =     []
        self.basis_dim = []

        if network_set.size != 0:
            self.network_set = network_set
            self.network_set_size = network_set.size

        else:
            print( "Please the network_set cannot be the empty set\n." )
            return 0


        if network_weight.shape[0] == network_weight.shape[1] and \
           network_weight.shape[0] > 0:

            self.network_weight = network_weight

        else:
            print( "Please the network_weight must be a square matrix and must not be empty\n." )
            return 0


    def Basis_of_the_vector_spaces_spanned_by_regular_paths(self):
        """
        Here we will be storing the regular paths of
        dimension 0, 1, ..., self.pph_dim + 1 into a list
        called basis:

             basis := [ B0, B_1, B_2, ..., B_(self.pph_dim + 1) ].

        The elements of the list basis are numpy arrays such
        that:
            B0: stores regular paths of dimension 0;
            B1: stores regular paths of dimension 1;
            B2: stores regular paths of dimension 2;
                             .
                             .
                             .

        The list Basis is an atribute of our class PPH.

        Important, the paths of dimension p will be represented
        by a sequence a_0, a_1, a_2, ..., a_p satisfying

             * a_0, a_1, a_2, ..., a_p in {0, 1, 2, ..., p}
             * and a_i != a_{i+1}, for all i in {0, 1, 2, ..., p-1}.

        Therefore, if a_i = 7, then a_i is merely an index pointer
        of the element 7 of the numpy.array network_set.

        With that in mind, suppose that we have
        network_set_size = 3, then

        B0 = numpy.array( [ [0], [1], [2] ] )
        B1 = numpy.array( [ [0,1], [0,2], [1,2], [1,0], [2,0], [2,1]] )
        B2 = numpy.array( [ [0,1,2], [1,2,0], [2,0,1], ... )

        This is merely a combinatorial task!
        """

        self.basis.append( np.zeros(self.network_set_size, 1) )
        for i in range( self.network_set ):
            self.basis[0][i,0] = i

        i = 1
        while i <= self.pph_dim + 1:
            self.basis.append( generating_all_regular_paths_dim_p( self.basis[i], self.network_set_size ) )
            i += 1


    def dimensions_of_each_vector_space_spanned_by_regular_paths(self):
        """
        Given the list self.basis of numpy arrays this method will
        store at the list self.basis_dim the dimensions of the
        vector spaces spanned by regular paths of dimension 0,1,2,...,
        self.dim + 1.
        """
        i = 0
        while i <= self.pph_dim + 1:
           self.basis_dim.append( self.basis.shape[0] )


    def allow_time (self, path_vector, path_dim):
        """
        Given the vector path_vector, an element of the vector
        space spanned by the regular paths of dimension path_dim.
        That is
            path_vector := sum_{i=0}^{self.basis_dim} alpha_i * basis[path_dim][i]
        with alpha_i in Z/2Z.

        Then this function will calculate the time of birth of
        such vector. In another words, this method
        returns the value
            min { k in A_k; path_vector in A_k, k >= = },
        where A_k is the set of the allowed regular paths of dimension
        k.

        Note tha path_vector is a numpy array that will eventually
        look like this, for instance
           path_vector = np.array([ 0, 0, 1, 1, 1, 0, 0, 1]),
        meaning that
           path_vector = self.basis[path_dim][2] + self.basis[path_dim][3] +
                         self.basis[path_dim][4] + self.basis[path_dim][7]
        """

        if path_dim == 0:
            return 0

        distance = []

        # Considering that
        # path_vector := sum_{i=0}^{self.basis_dim[path_dim]} alpha_i * self.basis[path_dim][i]
        find_indexes = np.arange( path_vector.size )[ path_vector == 1 ]

        # Let's find the elements of the basis that generate path_vector.
        # We will write vector_path = sum_i alpha_i * sigma_i,
        # with alpha_i in Z/2Z and sigma_i in self.basis[path.dim]

        for i in find_indexes:
            j = 0
            sigma_i = self.basis[path_dim][i] # sigma_i will look like: np.array( [ 0,2, 5,1, 3 ] ), for  instance

            # now we will run through the the vertices of the path sigma_i and we will
            # store the time needed for the edges to appear.
            while j < path_dim:
                distance.append( self.network_weight( sigma_i[j], sigma_i[j+1] ) )
                j += 1

        return max( distance )


    def entry_time(self, path_vector, path_dim, index_base ):

        if path_dim == 0:
            return 0

        elif path_dim == 1:
            return self.allow_time( path_vector, path_dim )

        else:
            distance = [ self.allow_time( path_vector, path_dim ) ]

            # Finding the basis vectors that generate
            # the vector path_vector.
            basis_that_generate_path_vector = (self.basis[ path_dim ])[ path_vector == 1 ]

            # Now we will apply the boundary transformation
            # d = sum_{j=0}^{path_dim} (-1)^j [a_0, a_1, ..., â_j, ...]

            # Taking
            #   path_vector = sum_i sigma_i
            # we will take the allow times of each d(sigma_i),
            # wth d() the boundary function

            for sigma_i in basis_that_generate_path_vector:
                i = 0
                while i <= path_dim:
                    aux_index = [ x != i for x in range(path_dim + 1) ]

                    # Now we will write sigma_i[ aux_index ] as
                    # a linear combination of the basis vectors of
                    # dimensionpath_dim - 1
                    # NOTATION:
                    #      aux_path := sigma_i[ aux_index ]
                    aux_path = np.zeros( path_dim - 1 )
                    for j in range( self.basis_dim[ path_dim - 1 ] ):
                        if np.all( self.basis[path_dim-1][j] == sigma_i[ aux_index ]):
                            aux_path[j] = 1


                    distance.append( self.allow_time( aux_path, dim -1 ) )
                    i += 1

            return max( distance )


       
#### constraints for the T_p structure (need to be improved)
###array_index = 0
###entry_index = 1
###allow_index = 2
###mark_index  = 3
###basis_index = 4



###################################
###### ---------> Tp structure
###### (Obs: it could be defined as a class here)
###################################
### 
###  T_p = []
###  i = 0
###
###  while i <= max_dimension_studied + 1:
###      j   = 0
###
###      # T_i = [ v_i, et(v_i), at(v_i), mark, basis_index]
###      T_i = []
###
###      while j < dimensionBasis[ i ]:
###          aux = np.zeros( dimensionBasis[i] )
###          aux[j] = 1
###
###          T_i.append( [ aux, entry_time( aux, i, j ), allow_time( aux, i ), False, j] )
###          j += 1
###
###      i += 1
###      T_p.append( T_i )
###
###
###  #############################################################
###  ### -----> Sorting the structures T_p and basis in agreement
###  ###        with the allow times.
###  #############################################################
###
###  for i in range( len( T_p ) ):
###      T_p[i].sort( key = lambda x: x[ allow_index ])
###
###      aux_index = [x[basis_index] for x in T_p[i]]
###      basis_i_copy = basis[i].copy()
###
###      for j in range( dimensionBasis[ i ] ):
###          basis[i][j] = basis_i_copy[ aux_index[j] ]
###
###
###  ###########################################################
###  ####### --------> Computing the persistence path diagram
###  ###########################################################
###
###  def BasisChange( an_array, dim ):
###
###      print('e aqui?')
###      p = dim # is dim really necessary or is implicit?
###
###      #u = np.zeros( dimensionBasis[p - 1] )
###      #for i in range( an_array.size ):
###      aux_basis = basis[ dim ][ an_array == 1  ][0]
###      i = 0
###      #aux_path = np.zeros( dimensionBasis[ dim - 1 ] )
###      u = np.zeros( dimensionBasis[ dim - 1 ] )
###      while i <= dim:
###          aux_index = [ x != i for x in range(dim+1) ]
###
###
###          for j in range( dimensionBasis[ dim - 1 ] ):
###              if np.all( basis[dim-1][j] == aux_basis[ aux_index ]) and T_p[dim-1][j][mark_index] == False:
###                  u[j] = 1
###          i += 1
###      #aux_index = [ x != i for x in range( an_array.size ) ]
###      #aux_array =  an_array[ aux_index ]
###      #print('dim, aux_array')
###      #print(dim, aux_array)
###
###      #for j in range( dimensionBasis[ p - 1 ] ):
###      #    if np.all( aux_array == basis[p][j] ) and  T_p[p - 1][j][mark_index] == False:
###      #        u[j] = 1
###
###
###
###
###      print(u)
###
###      while np.all( u != 0 ):
###          print('chegou aqui?\n')
###          aux_sigma     = np.arange( u.size )
###          aux_index_eq1 = (aux_sigma[ u == 1 ])
###
###
###          aux_index_max = aux_index_eq1.max() # equivalent to i in the paper
###          sigma = basis[ p -1 ][aux_index_max]
###
###          et = max( [allow_time(an_array, p), allow_time(sigma, p-1)] )
###
###          if  T_p[p - 1][ aux_index_max ][ array_index ][aux_index_max] == 0 :
###              break
###
###          u_next = u ^ T_p[ p-1 ][ aux_index_max ][ array_index ]
###
###
###      return [u_next, aux_index_map, et]
###
###
###  Pers = [ [], [], [] ]
###
###
###  for p in range( max_dimension_studied + 1): # max_dimension_studied + 1
###                                              # because range returns
###                                              # a interval like [a,b)
###
###      j = 0
###      print('aui000')
###      while j < dimensionBasis[ p + 1 ]:
###          return_BasisChange = BasisChange( basis[p+1][j], p+1 )
###          print('aqui')
###          print(return_BasisChange)
###
###          u = return_BasisChange[0]
###          i = return_BasisChange[1]
###          et = return_BasisChange[2]
###
###          if np.all( u == 0 ):
###              T_p[ p + 1 ][j][mark_index] = True
###
###          else:
###              T_p[p][i][ array_index ] = u
###              T_p[p][i][ entry_index ] = et
###
###              Pers[p].append( [T_p[p][i][entry_index], et ] )
###
###          j += 1
###
###      j = 0
###      while j < dimensionBasis[ p ]:
###          if T_p[ p ][j][ mark_index ] == True and \
###             np.all( T_p_[ p ][j][ array_index ] == 0):
###              Pers[p].append( [T_p[p][j][ entry_index ], np.inf] )
###
###          j += 1
###
###  print( Pers )
