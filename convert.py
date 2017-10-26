import numpy as np
from scipy.sparse.csgraph import dijkstra
import scipy.io
from scipy.sparse import csr_matrix
import scipy.sparse
import multiprocessing as mp
import os.path
import sys

	

def mm2csr(ifile):

    M=scipy.io.mmread(ifile)
    theShape=M.shape
    neqns=theShape[0]
    nnz=M.nnz
    if theShape[0]!=theShape[1]:
        raise Exception('Matrix not square: {}'.format(ifile))

    print("Loaded: {0} (size {1:,}, {2:,} non-zeros)".format(ifile, neqns, nnz))
    fname, ext = os.path.splitext(ifile)

    ofile = fname + ".csr" 
    M_csr = scipy.sparse.csr_matrix(M)

    M_csr.sort_indices()

    values = M_csr.data
    column_indices = M_csr.indices
    row_pointers = M_csr.indptr

    fhdl=open(ofile, 'w')

    fhdl.write(str(theShape[0]))
    fhdl.write('\n')

    fhdl.write(str(M_csr.nnz))
    fhdl.write('\n')

    for i in row_pointers:
        fhdl.write(str(i+1))
        fhdl.write('\n')

    for j in column_indices:
        fhdl.write(str(j+1))
        fhdl.write('\n')

    for m in values:
        fhdl.write(str(m))
        fhdl.write('\n')

    fhdl.close()
    print('Written: {}'.format(ofile))

if __name__ == "__main__":

    inames=sys.argv[1:]
    print("Converting the following to CSR:")
    print('\n'.join(inames))
    print("--------------\n")

    myPool = mp.Pool(processes=2)
    myPool.map(mm2csr, inames)
