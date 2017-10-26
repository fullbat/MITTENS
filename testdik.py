import numpy as np
from scipy.sparse.csgraph import dijkstra
import scipy.io
from scipy.sparse import csr_matrix
import scipy.sparse
import multiprocessing as mp
import os.path
import sys

matdict = scipy.io.loadmat('mgh_voxel_graph.mat')

graph = csr_matrix(matdict['graph'])
#graph.nonzero()
	


dist = []
path = []
t = 0
i = matdict['nvoxels'][0][0]-1
nvoxel = matdict['nvoxels'][0][0]




dist = distances
print(dist)

while t != nvoxel:
 distances, predecessors = dijkstra(graph, indices=(0,nvoxel-t), return_predecessors=True, directed=False)
   while i != 0:
            path.append(predecessors[i])
            i = predecessors[i]

            #path.append(value[0])
 print("path=",path[::-1])
 t+=1


