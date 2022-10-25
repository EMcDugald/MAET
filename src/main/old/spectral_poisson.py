import numpy as np
import mymath
from scipy import signal
import scipy as sp
import utils

#solves for coefficients in fourier transform
# see https://en.wikipedia.org/wiki/Spectral_method for an example
def get_bks(fft_mat):
    k1_ind = np.shape(fft_mat)[0]
    k2_ind = np.shape(fft_mat)[1]
    U = np.zeros(shape=(k1_ind,k2_ind))
    for k1 in range(1,k1_ind):
        for k2 in range(1,k2_ind):
            U[k1,k2] = -fft_mat[k1,k2]/((k1)**2+(k2)**2)
    return U

def get_bks_lap(fft_mat):
    k1_ind = np.shape(fft_mat)[0]
    k2_ind = np.shape(fft_mat)[1]
    U = np.zeros(shape=(k1_ind,k2_ind))
    for k1 in range(1,k1_ind):
        for k2 in range(1,k2_ind):
            U[k1,k2] = -fft_mat[k1,k2]*((k1)**2+(k2)**2)
    return U


#returns un-normalized solution from laplacian
#in_dim = inout dimensions- defines normalization
def return_u(lap,sq_len):
    #normalization for inverse dst
    norm = (sq_len/np.pi)**2
    # get coefficients of dft of laplacian u
    ak_mat = mymath.mydst2(lap)

    # solve for the coefficients of u in frequency domain
    bk_mat = get_bks(ak_mat)

    # apply the dft to get the time/space domain values of u
    u = mymath.mydst2(bk_mat)
    return u*norm

#in_dim = inout dimensions- defines normalization
def return_lapu(u,sq_len):
    #normalization for inverse dst
    norm = (sq_len/np.pi)**2
    # get coefficients of dft of laplacian u
    ak_mat = mymath.mydst2(u)

    # solve for the coefficients of u in frequency domain
    bk_mat = get_bks_lap(ak_mat)

    # apply the dft to get the time/space domain values of u
    u = mymath.mydst2(bk_mat)
    return u/norm


############# WIP ###############


#smooth transition function
def zeta(x1,x2,rmin,rmax):
    return mymath.t6hat(rmax,rmin,x1)*mymath.t6hat(rmax,rmin,x2)

#step function
def chi(x1, x2, rmin):
    t1 = np.abs(x1) <= rmin
    t2 = np.abs(x2) <= rmin
    t3 = np.logical_and(t1, t2)
    p = np.where(t3 == True, 0, 1)
    return p

#finite difference approx to laplacian
def laplacian(dx, dy, w):
    laplacian_xy = np.zeros(w.shape)
    for y in range(w.shape[1]-1):
        laplacian_xy[:, y] = (1/dy)**2 * ( w[:, y+1] - 2*w[:, y] + w[:, y-1] )
    for x in range(w.shape[0]-1):
        laplacian_xy[x, :] = laplacian_xy[x, :] + (1/dx)**2 * ( w[x+1,:] - 2*w[x,:] + w[x-1,:] )
    return laplacian_xy


#fundametnal solution- small shift in arg to bypass division by zero
def phi(x1,x2):
    norm = np.sqrt(x1**2+x2**2)
    return np.log(norm+10e-50)/(2*np.pi)



def G(u,x1,x2,rmin,rmax):
    sq_len = np.abs(np.max(x1[0])-np.min(x1[0]))
    w = u*zeta(x1,x2,rmin,rmax)
    return return_lapu(w,sq_len)


# returns convolution xG*phi on big square- based on quadrature
def convol_quad(rmin,rmax,X1,X2,u):
    x1_dim = np.shape(X1)[0]
    x2_dim = np.shape(X2)[0]
    conv_arr = np.zeros((x1_dim,x2_dim))
    dx = X1[0][1]-X1[0][0]
    f = chi(X1,X2,rmin)*G(u,X1,X2,rmin,rmax)
    for i in range(1,x1_dim):
        for j in range(1,x2_dim):
            g = phi(X1[i,j]-X1,X2[i,j]-X2)
            conv_arr[i,j] = np.sum(f*g)*dx**2
    return conv_arr


def v_quad(X1,X2,rmin,rmax,u):
    xgp = convol_quad(rmin,rmax,X1,X2,u)
    return u*zeta(X1,X2,rmin,rmax) - xgp

# returns convolution xG*phi on big square- based on FFT
def convol(rmin,rmax,X1,X2,u):
    ker = chi(X1,X2,rmin)*G(u,X1,X2,rmin,rmax)
    k_start = X1[0][0]
    k_end = -k_start
    k_size = np.shape(X1)[0]
    im_start = int(2*k_start)
    im_end = int(2*k_end)
    im_len = int(2*k_size-1)
    im_sp = np.linspace(im_start,im_end,im_len)
    IM1, IM2 = np.meshgrid(im_sp,im_sp)
    im = phi(IM1,IM2)
    conv = signal.fftconvolve(im,ker,mode='valid')
    dx = X1[0][1]-X1[0][0]
    return conv*dx**2

def v(X1,X2,rmin,rmax,u):
    xgp = convol(rmin,rmax,X1,X2,u)
    return u*zeta(X1,X2,rmin,rmax) - xgp


# returns convolution xG*phi on big square- based on FFT
# does not rely on scipy method
def convol2(rmin,rmax,X1,X2,u):
    # set scaling for continuous conversion
    dx = X1[0][1] - X1[0][0]
    # get kernel on [-a,a]^2
    ker = chi(X1, X2, rmin) * G(u, X1, X2, rmin, rmax)
    utils.ndwrite(ker,"XG.d")
    # make meshgrid on [-2a,2a]^2
    k_start = X1[0][0]
    k_end = -k_start
    k_size = np.shape(X1)[0]
    im_start = 2 * k_start
    im_end = 2 * k_end
    im_len = 2 * k_size - 1
    im_sp = np.linspace(im_start, im_end, im_len)
    IM1, IM2 = np.meshgrid(im_sp, im_sp)
    # het image on [-2a,2a]^2
    im = phi(IM1, IM2)
    utils.ndwrite(im,"phi.d")
    # shape of kernel
    m1, m2 = np.shape(ker)
    # shape of image
    n1, n2 = np.shape(im)
    # total size of lienar convolution
    num_rows = m1 + n1 - 1
    num_cols = m2 + n2 - 1
    k_pad = np.zeros((num_rows, num_cols))
    i_pad = np.zeros((num_rows, num_cols))
    # place kernel in center of pad
    k_row_st = round((num_rows - m1) / 2)
    k_col_st = round((num_cols - m2) / 2)
    i_row_st = round((num_rows - n1) / 2)
    i_col_st = round((num_cols - n2) / 2)
    k_pad[k_row_st:k_row_st + m1, k_col_st:k_col_st + m2] = ker
    utils.ndwrite(k_pad,"XG_pad.d")
    i_pad[i_row_st:i_row_st + n1, i_col_st:i_col_st + n2] = im
    utils.ndwrite(i_pad,"phi_pad.d")
    #take inverse fft of product of ffts
    ker_pad_shift_fft = sp.fft.fft2(sp.fft.fftshift(k_pad))
    utils.ndwrite(ker_pad_shift_fft, "XG_shift_fft.d")
    im_pad_fft = sp.fft.fft2(i_pad)
    utils.ndwrite(im_pad_fft, "phi_fft.d")
    conv = sp.fft.ifft2(ker_pad_shift_fft * im_pad_fft)
    utils.ndwrite(np.real(conv), "ifft_fft_prod.d")
    #extract middle of array of size kernel
    conv = np.real(conv)[k_row_st-1:2*k_row_st, k_col_st-1:2*k_col_st]
    utils.ndwrite(conv, "ifft_fft_prod_mid.d")
    return conv * dx ** 2

def v2(X1,X2,rmin,rmax,u):
    xgp = convol2(rmin,rmax,X1,X2,u)
    return u*zeta(X1,X2,rmin,rmax) - xgp

