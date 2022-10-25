import numpy as np
import scipy as sp
import utils

#Compute partial derivative of Gaussian

n=8
a=5.
ndim = 2**n+1
x1 = np.linspace(-a,a,ndim,dtype='complex_')
x2 = np.linspace(-a,a,ndim,dtype='complex_')
X1, X2 = np.meshgrid(x1,x2)
nm1 =ndim -1
X1m1 = X1[0:nm1,0:nm1]
X2m1 = X2[0:nm1,0:nm1]

sigma = 0.25
#comute laplacian
def gaussian(x1,x2,sigma):
    exponent = -(x1**2+x2**2)/(sigma*sigma)
    return np.exp(exponent)

#partial derivative
# i=0 -> x, i=1 -> y
def partial_gauss(x1,x2,sigma,i):
    exponent = -(x1**2+x2**2)/(sigma*sigma)
    if i==0:
        return -2*x1*np.exp(exponent)/sigma**2
    else:
        return -2*x2*np.exp(exponent)/sigma**2


######################## Basic error- Compute the gaussian via FFT ##############################
#f
f = gaussian(X1m1,X2m1,sigma)

#shifted f
shifted_f = sp.fft.fftshift(f)

#fft of shifted f
shifted_fhat = sp.fft.fft2(shifted_f)

#inverse transform of fft of shifted f
f_tst_shifted = sp.fft.ifft2(shifted_fhat)

#shifted inverse transform of fft of shifted f (same as f)
f_tst = sp.fft.fftshift(f_tst_shifted)

utils.ndwrite(f,'f_exact.d')
utils.ndwrite(f_tst,'f_approx.d')


######################## Basic error- compute f_x and f_y via FFT ##############################

#exact partials
f_x_exact = partial_gauss(X1m1,X2m1,sigma,0)
f_y_exact = partial_gauss(X1m1,X2m1,sigma,1)

#frequency space meshgrid
kx = (np.pi/a)*sp.fft.fftfreq(nm1,1./nm1)
ky = (np.pi/a)*sp.fft.fftfreq(nm1,1./nm1)
KX,KY = np.meshgrid(kx,ky)

#shifted f
fshift = sp.fft.fftshift(f)
#fft of shifted f
fshift_hat = sp.fft.fft2(fshift)

#partials in frequency space
fshift_hat_x = 1j*KX*fshift_hat
fshift_hat_y = 1j*KY*fshift_hat

#partials of shifted f in spatial domain
fshift_x = sp.fft.ifft2(fshift_hat_x)
fshift_y = sp.fft.ifft2(fshift_hat_y)

#rearrange the partials
f_x = sp.fft.fftshift(np.real(fshift_x))
f_y = sp.fft.fftshift(np.real(fshift_y))

utils.ndwrite(f_x,'f_x.d')
utils.ndwrite(f_y,'f_y.d')
utils.ndwrite(f_x_exact,'f_x_exact.d')
utils.ndwrite(f_y_exact,'f_y_exact.d')