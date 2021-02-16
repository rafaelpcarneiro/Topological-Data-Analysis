#!/usr/bin/env python3

import numpy as np
from persistent_path_homology import *

##### First Test
def euclidean_distance( x, y ):
    #x, y numpy arrays

    return sum(x-y)**2

network_set = np.random.uniform( 0,1, (10,1) )

network_weight = np.zeros( (10,10) )
for i in range( 5 ):
    for j in range( 5 ):
        network_weight[i,j] = euclidean_distance( network_set[i], \
                                                  network_set[j])

test1 = PPH( network_set, network_weight, 1 )

test1.ComputePPH()
