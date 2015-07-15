import numpy
import scipy.sparse
import scipy.sparse.linalg
import argparse
import sys

parser = argparse.ArgumentParser(description='Compute eignevalues and eigenvectors')
parser.add_argument('--input', dest='input', default='', type=str,
                    help='Sparse matrix coordinate file in i, j, value format')
parser.add_argument('--numEigen', dest='numEigen', default=1, type=int, 
                    help='Number of eigenvalues and eigenvectors')

args = parser.parse_args()

# check
if args.input == '':
    print 'ERROR: must specify input file (--input INPUT)'
    sys.exit(1)

data = numpy.loadtxt(args.input)
ivec = numpy.array(data[:,0], numpy.int)
jvec = numpy.array(data[:,1], numpy.int)
matrix = scipy.sparse.coo_matrix( (data[:,2], (ivec, jvec)) )
eigenvalues, eigenvectors = scipy.sparse.linalg.eigsh(matrix, args.numEigen, which='LA')

print 'eigenvalues: '
print eigenvalues



