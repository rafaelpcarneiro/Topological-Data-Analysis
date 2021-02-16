#!/usr/bin/env python3

import numpy as np
from persistent_path_homology import *

DATA_SIZE = 3
PPH_DIM   = 1

##### First Test
def euclidean_distance( x, y ):
    #x, y numpy arrays

    return sum(x-y)**2

np.random.seed(500)
network_set = np.random.uniform( 0,1, (DATA_SIZE, 1) )

network_weight = np.zeros( (DATA_SIZE, DATA_SIZE) )
for i in range( DATA_SIZE ):
    for j in range( DATA_SIZE ):
        network_weight[i,j] = euclidean_distance( network_set[i], \
                                                  network_set[j])

test1 = PPH( network_set, network_weight, PPH_DIM )

test1.ComputePPH_printing_step_by_step()
