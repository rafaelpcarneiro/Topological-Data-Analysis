#!/usr/bin/env python3

import numpy as np
from persistent_path_homology import *

##### First Test
def euclidean_distance( x, y ):
    #x, y numpy arrays

    return sum(x-y)**2

network_set = np.random.uniform( 0,1, (5,2) )

network_weight = np.zeros( (5,5) )
for i in range( 5 ):
    for j in range( 5 ):
        network_weight[i,j] = euclidean_distance( network_set[i], \
                                                  network_set[j])

test1 = PPH( network_set, network_weight, 1 )

test1.Basis_of_the_vector_spaces_spanned_by_regular_paths()
test1.dimensions_of_each_vector_space_spanned_by_regular_paths()
test1.generating_T_p()
