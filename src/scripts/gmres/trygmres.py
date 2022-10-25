import numpy as np
from scipy.sparse.linalg import gmres

print('start')

ndim = 9

#---------- making a Vandermonde matrix, just for fun -------------

A = np.zeros((ndim,ndim), dtype=float)

row = np.linspace(0.,1,ndim+1)
row = row[1:ndim+1]
print(row)

A[0,:] = row
for j in range(ndim-1):
    A[j+1,:] = A[j,:]*row

print('--------- The matrix: ')
print(A)


xgood = np.array(range(ndim), dtype = float)

xgood = xgood[::-1]

print('--------- Exact solution: ')
print(xgood)

b = A.dot(xgood)

print('--------- Right hand side: ')
print (b)

#------------------------------------------------------------------ call GMRES ---------------------
x, exitCode = gmres(A, b,maxiter=1,restart=25,tol = 1.e-15)


print(' ')
print('Exit code = ',exitCode)
print(' ')

print('--------- Found solution:')
print(x)


aa = A.dot(x)
print('--------- Discrepancy in the RHS')
print (aa-b)

print(' ')
print('--------- Error in the solution')
print (x-xgood)

