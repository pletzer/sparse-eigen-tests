import numpy
import scipy.sparse
import scipy.sparse.linalg

filename = '../data/10x10.dat'
numEigenvalues = 5

data = numpy.loadtxt(filename)
ivec = numpy.array(data[:,0], numpy.int)
jvec = numpy.array(data[:,1], numpy.int)
matrix = scipy.sparse.coo_matrix( (data[:,2], (ivec, jvec)) )
eigenvalues, eigenvectors = scipy.sparse.linalg.eigsh(matrix, numEigenvalues, which='LA')

print eigenvalues



