import numpy as np
import matplotlib.pyplot as plt
import spectral_poisson as spec_poi
import utils

#-----------------------------------------------------------------------
# This file outputs 5 d files and 5 surface plits
# d files- exact u on big square, numerical u on big square
#          exact u on small square, numerical u on small square, numerical v on small square
# plots - exact u on big square, numerical u on big square,
#         numerical uz on big square, numerical lapuz on big square, numerical v on big square
#------------------------------------------------------------------------




#Solving lap u = gaussian on "big" square
n=7
a=5.
#dimension of smaller square
small_dim = 2**n+1

#set up big square
big_dim = 2*small_dim-1
x1_bg = np.linspace(-a,a,big_dim)
x2_bg = np.linspace(-a,a,big_dim)
X1_bg, X2_bg = np.meshgrid(x1_bg,x2_bg)


sigma = 0.25


#comute laplacian
def lap_gauss(x1,x2,sigma):
    scale = 4*(x1**2+x2**2-sigma*sigma)/(sigma**4)
    exponent = -(x1**2+x2**2)/(sigma*sigma)
    return scale*np.exp(exponent)


#compute known solution
def gaussian(x1,x2,sigma):
    exponent = -(x1**2+x2**2)/(sigma*sigma)
    return np.exp(exponent)


#--solve u via Fourier given its laplacian
lapu = lap_gauss(X1_bg,X2_bg,sigma)
u = spec_poi.return_u(lapu,2*a)
exact_sol = gaussian(X1_bg,X2_bg,sigma)

#constructing functions for solution v

#z is smooth transition func- 1 on small,transition to 0 on big
z = spec_poi.zeta(X1_bg,X2_bg,a/2,a)


#uz is numerically solved u * z
uz = u*z

# laplacian of uz
lapuz = spec_poi.return_lapu(uz,2*a)

#solution v
v = spec_poi.v(X1_bg,X2_bg,a/2,a,u)


#extract small square from big square
mid1 = int(.5*(big_dim-small_dim))
mid2 = int(mid1 + small_dim)

u_small= u[mid1:mid2,mid1:mid2]
v_small= v[mid1:mid2,mid1:mid2]

utils.ndwrite(exact_sol[mid1:mid2,mid1:mid2],"gauss_exact_small.d")
utils.ndwrite(v_small,"gauss_v_small.d")
