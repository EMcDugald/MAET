import numpy as np
import spectral_poisson as spec_poi
import mymathradial as mmr
import utils

import time
start = time.time()

#Solving lap u = lap radial on square
n=9
a=3

#dimension of smaller square
small_dim = 2**n+1

#set up big square
big_dim = 2*small_dim-1
x1_bg = np.linspace(-a,a,big_dim)
x2_bg = np.linspace(-a,a,big_dim)
X1_bg, X2_bg = np.meshgrid(x1_bg,x2_bg)


#--solve u via Fourier given its laplacian
exact_sol, lapu = mmr.radial(big_dim,0.,0.,.9,.3,a)
u = spec_poi.return_u(lapu,2*a)

#--solutions on big square
utils.ndwrite(exact_sol,"exact_radial1.d")
utils.ndwrite(lapu,"lap_exact_radial1.d")


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

exact_small = exact_sol[mid1:mid2,mid1:mid2]
v_small = v[mid1:mid2,mid1:mid2]

utils.ndwrite(exact_small,"radial_exact_small_fft1.d")
utils.ndwrite(v_small,"radial_v_small_fft1.d")

end = time.time()
print("Elapsed Time for FFT Based Computation:",end-start)
