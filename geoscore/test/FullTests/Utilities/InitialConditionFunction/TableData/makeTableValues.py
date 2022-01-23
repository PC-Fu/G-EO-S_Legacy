import numpy as np
import pylab as pl

data = np.loadtxt('Table_values.txt')
[nx,ny] = data.shape

X,Y = np.mgrid[0.:nx,0.:ny]

data = np.sin( 2*np.pi*X/float(nx) )*np.cos( 4*np.pi*Y/float(ny) )

np.savetxt('Table_values.txt',data)
