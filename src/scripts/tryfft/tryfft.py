import numpy as np
import scipy as sp
import utils


print('start')

def gaussian(x1,x2,sigma):
    exponent = -(x1**2+x2**2)/(sigma*sigma)
    return np.exp(exponent)


a = 3.
n = 8

sigma = 0.25
nsize = 2**n+1
nm1 = nsize - 1

x1 = np.linspace(-a/2,a/2,nsize)
x2 = np.linspace(-a/2,a/2,nsize)
X1, X2 = np.meshgrid(x1,x2)
X1m1 = X1[0:nm1,0:nm1]
X2m1 = X2[0:nm1,0:nm1]

#------------- make a function to play with --------------------

func = gaussian(X1m1,X2m1,sigma)
utils.ndwrite(func,"func.d")

#------------- shifted version of a function -------------------

shifted_func = sp.fft.fftshift(func,[0,1])
utils.ndwrite(shifted_func,"sh_func.d")

#------------- Fourier transform of the function ---------------
ftransform = sp.fft.fft2(shifted_func)

utils.ndwrite(np.real(ftransform),"ftr_r.d")
utils.ndwrite(np.imag(ftransform),"ftr_i.d")

#------------- shifted Fourier transform of the function --------

ftr_sh = sp.fft.fftshift(ftransform,[0,1])

utils.ndwrite(np.real(ftr_sh),"ftr_sh_r.d")
utils.ndwrite(np.imag(ftr_sh),"ftr_sh_i.d")

#------------- inverse Fourier transform of the unshifted transform --------

recovered_func = sp.fft.ifft2(ftransform)

utils.ndwrite(np.real(recovered_func),"func_r.d")
utils.ndwrite(np.imag(recovered_func),"func_i.d")

