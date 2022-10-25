import numpy as np
from scipy.sparse.linalg import gmres
from scipy.sparse.linalg import LinearOperator

iter = 0

def myvandermonde():
   A = np.zeros((ndim,ndim), dtype=float)

   row = np.linspace(0.,1,ndim+1)
   row = row[1:ndim+1]

   A[0,:] = row
   for j in range(ndim-1):
      A[j+1,:] = A[j,:]*row

   return A

def my_lin_op(v):
   global iter

   iter = iter + 1
   print('Iteration # ',iter)
   return A.dot(v)

print('start')

ndim = 9
iter = 0

#---------- making a Vandermonde matrix, just for fun -------------


A = myvandermonde()

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


AA = LinearOperator((ndim,ndim), matvec = my_lin_op)


x, exitCode = gmres(AA, b,maxiter=1,restart=40,tol = 1.e-15)


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
print('Max error in x :',np.max(np.abs(x-xgood)))

