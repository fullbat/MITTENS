#!/usr/bin/python3

from networkit.graph import Graph
from tqdm import tqdm
from time import time
import networkx as nx
import networkit
from scipy.io.matlab import loadmat, savemat
from fib_io import load_fibgz
from importlib import reload
import logging
import numpy as np
from networkit.distance import Dijkstra
from spatial import Spatial, hdr
#from scipy.sparse.csgraph import dijkstra
import scipy.io
from scipy.sparse import csr_matrix
import scipy.sparse
import multiprocessing as mp
import os.path
import sys

DISCONNECTED=999999999

"""
matdict = scipy.io.loadmat('mgh_voxel_graph.mat')

graph = csr_matrix(matdict['graph'])
#graph.nonzero()


distances, predecessors = dijkstra(graph, indices=(0,10), return_predecessors=True, directed=False)	


#value = graph.data
#column_index = graph.indices
#row_pointers = graph.indptr



dist = []
path = []

i = matdict['nvoxels'][0][0]-1

dist = distances
print(dist)


while i != 0:
    
    path.append(predecessors[i])
    i = predecessors[i]

#path.append(value[0])
print("path=",path[::-1])
"""


class SS:
    
    
    def __init__(self, matfile="mgh_voxel_graph.mat",fibgz_file="mgh_1001.src.gz.odf8.f5rec.bal.gqi.1.25.fib.gz", nifti_prefix="",
            real_affine_image="", mask_image="",
            step_size=np.sqrt(3)/2. , angle_max=35, odf_resolution="odf8", 
            weighting_scheme=None,
            # Spatial mapping data
            flat_mask=None, nvoxels=None,
            real_affine=None, ras_affine=None, voxel_size=None,volume_grid=None,
            voxel_coords=None, coordinate_lut=None, graph=None, angle_weights="flat", angle_weighting_power=1.,normalize_doubleODF=True):
           
               
            
        
        # These will get filled out from loading a fibgz or niftis
        self.flat_mask = None
        self.nvoxels = None
        self.voxel_size = voxel_size
        self.voxel_coords = None
        self.coordinate_lut = None
        self.label_lut = None
        self.atlas_labels = None
        self.mask_image = mask_image
        # From args
        self.step_size = step_size
        self.odf_resolution = odf_resolution
        self.angle_max = angle_max
        self.orientation = None
        self.angle_weights = angle_weights
        self.normalize_doubleODF = normalize_doubleODF
        self.angle_weighting_power = angle_weighting_power
        self.graph = graph
        
        #presi da voxel_graph
        self.step_size = step_size
        self.angle_max = angle_max
        self.odf_resolution = odf_resolution
        self.weighting_scheme = weighting_scheme
        self.angle_weights = angle_weights
        self.normalize_doubleODF = normalize_doubleODF
        self.angle_weighting_power = angle_weighting_power

        self.flat_mask = flat_mask
        self.nvoxels = nvoxels
        self.real_affine = real_affine
        self.ras_affine = ras_affine
        self.volume_grid = volume_grid

        self.voxel_coords = voxel_coords
        self.coordinate_lut = coordinate_lut
        self.graph = graph

        # These are computed later
        self.atlas_labels = None
        self.label_lut = None
        self.background_image = None
        self.null_graph = None
        
        # Nodes that do not correspond to a voxel (eg source or sink)
        self.nonvoxel_nodes = set()
        self.undirected_component_ids = None
        print("\nUsing\n------\n  Step Size:\t\t%.4f Voxels \n  ODF Resolution:\t"
                "%s\n  Max Angle:\t\t%.2f Degrees\n"
                "  Angle Weights:\t%s\n  Angle weight power:\t%.1f",
                self.step_size, self.odf_resolution, self.angle_max, 
                self.angle_weights,self.angle_weighting_power)
        print("Ao")
        #m = loadmat(matfile)
        
        
        #graph = csr_matrix(m['graph'])
        
        #sparse_graph = nx.from_scipy_sparse_matrix(m['graph'],create_using=nx.DiGraph(),
                                #edge_attribute= 'weight')
        
        #self.graph = networkit.nxadapter.nx2nk(sparse_graph, weightAttr='weight')
               
        #self._load_fibgz(fibgz_file)
        
        if matfile:
            print("Loading voxel graph from matfile")
            self._load_matfile(matfile)
            
        
    def _load_matfile(self,matfile):
        if not os.path.exists(matfile):
            print("No such file: %s",matfile)
            return

        #try:
        m = loadmat(matfile)
        #except Exception, e:
        #    logger.critical("Unable to load %s:\n %s", matfile, e)
        #    return

        # Transition prob details
        self.step_size = float(m['step_size'].squeeze())
        self.angle_max = float(m['angle_max'].squeeze())
        self.odf_resolution = str(m['odf_resolution'].squeeze())
        self.weighting_scheme = str(m['weighting_scheme'].squeeze())
        self.angle_weights = str(m['angle_weights'].squeeze())
        self.normalize_doubleODF = bool(m['normalize_doubleODF'].squeeze())
        self.angle_weighting_power = float(m['angle_weighting_power'].squeeze())
        
        # Spatial mappings
        self.flat_mask = m['flat_mask'].squeeze().astype(np.bool)
        masked_voxels = self.flat_mask.sum()
        assert masked_voxels == m['nvoxels']
        self.ras_affine = m['ras_affine']
        self.real_affine = m.get('real_affine', self.ras_affine).squeeze()
        self.voxel_size =  m['voxel_size'].squeeze()
        self.volume_grid = m['volume_grid'].squeeze().astype(np.int)

        # Guaranteed to have a flat mask by now
        self.nvoxels = masked_voxels
        # These coordinates are LPS+ voxels
        self.voxel_coords = np.array(np.unravel_index(
            np.flatnonzero(self.flat_mask), self.volume_grid, order="F")).T
        self.coordinate_lut = dict(
            [(tuple(coord), n) for n,coord in enumerate(self.voxel_coords)])

        print("Loading graph from matfile")
        sparse_graph = nx.from_scipy_sparse_matrix(m['graph'],create_using=nx.DiGraph(),
                                edge_attribute="w")
        self.graph = networkit.nxadapter.nx2nk(sparse_graph, weightAttr="w")
    
    def _load_fibgz(self, path):
        print("Loading %s", path)
        f = load_fibgz(path)
        print("Loaded %s", path)

        
    def voxelvalue(self):
        vox = self.graph.randomNode()
        print(vox)
        
        
        return vox, ""
    
    
  
    def region_voxels_to_region_query(self,  write_trk="", 
            write_prob="", write_nifti="", write_wm_maxprob_map=""):
        
        
        """
        Query paths from one region to another. Parameters ``from_region``
        and ``to_region`` can be a path to a nifti file or a region ID
        if an atlas has been added through the ``add_atlas`` function.
        The probability of being connected to the ``to_region`` is calculated
        for every voxel in the ``from_region``.

        Parameters:
        ===========

        from_region:str,int
          Region from which the probability will be calculated for
          every voxel along the shortest path to any voxel in the
          ``to_region``.
        
        to_region:str,int
          The "sink" region. 

        write_trk:str
          Suffix that will be added to the trackvis output. If empty a 
          trackvis file won't be written.

        write_prob:str
          suffix for a txt file that contains the probability for each 
          path. These files can be loaded into DSI Studio to color the
          streamlines.

        write_nifti:str
          Write the probabilities to a NIfTI file. Each voxel in ``from_region``
          will contain its probability.

        write_wm_maxprob_map:str
          Writes a map where each voxel contains the probability of the maximum
          probability path that goes through it.
        """
        # Get lists of nodes for the query
        from_nodes, from_name = self.voxelvalue()
        
        to_nodes, to_name = self.voxelvalue()
        
        # Loop over all the voxels in the from_region
        sink_label_node = to_nodes


        # Find the connected components
        undirected_version = self.graph.toUndirected()
        components = networkit.components.ConnectedComponents(undirected_version)
        components.run()
        print("Found %d components in the graph", components.numberOfComponents())
        target_component = components.componentOfNode(sink_label_node)

        n = networkit.distance.Dijkstra(self.graph, sink_label_node)
        t0 = time()
        n.run()
        t1 = time()
        print("Computed shortest paths in %.05f sec", t1-t0)

        trk_paths = []
        probs = np.zeros(self.nvoxels,dtype=np.float)
        maxprobs = np.zeros(self.nvoxels,dtype=np.float)
        #for node in tqdm(from_nodes):
        node = from_nodes
        if components.componentOfNode(node) == target_component:
            path = n.getPath(node)
            print("path: ", path)
            if not len(path):  #if len(path) is not None 
                score = self.get_path_probability(self.graph, path)
                print("score: ", score)
                probs[node] = score
                
            
            if write_trk:
                trk_paths.append(path[1:-1])
                if write_wm_maxprob_map:
                        path = np.array(path)
                        path = path[path < self.nvoxels]
                        maxprobs[path] = np.maximum(maxprobs[path], score)

        # Write outputs
        if write_prob:
            g = open("%s_to_%s_%s.txt"%(from_name, to_name, write_trk), "w")
            for prob in probs:
                #print("probssss: ", prob)
                g.write("%.9f\n"%prob)
            g.close()
        if write_trk:
            if write_trk.endswith(".trk") or write_trk.endswith(".trk.gz"):
                self._write_trk(trk_paths, write_trk)
            elif write_trk.endswith("txt"):
                self._write_trk_txt(trk_paths, write_trk)
        if write_nifti:
            self.save_nifti(probs, write_nifti)
        if write_wm_maxprob_map:
            self.save_nifti(maxprobs, write_wm_maxprob_map)

        return trk_paths, probs
    
    def Dijkstra(self, g, source, sink):
        d = networkit.distance.Dijkstra(g, source, target = sink)
        d.run()
        path = d.getPath(sink)
        return path 
    
    def get_path_probability(self, g, path):
        prob = self.get_prob(g, path)
        prob = prob ** (1./len(path))
        print("get_path_probability: ", prob)
        return prob

    
    def get_prob(self, g, path):
        prob = 1 
        for step in range(len(path) - 1):
            prob*=np.e**(-g.weight(path[step], path[step+1]))
            print("get_prob: ", prob)
        return prob
    
    
    def _write_trk(self,paths, trk_file_name):
        if not len(paths): 
            print("Empty paths, not writing %s",trk_file_name)
            return
        header = hdr.copy()
        print("voxel_size", self.voxel_size)
        header['voxel_size'] = self.voxel_size.astype("<f4")
        header['dim'] = self.volume_grid.astype('<i2')
        trk_paths = [
            (self.voxel_coords[np.array(path[1:-1])]*self.voxel_size, 
                None, None) for path in paths ]
        nib.trackvis.write(trk_file_name, trk_paths, header )





