import numpy as np
import mymath
import scipy as sp
import time
from scipy import ndimage

#cutting function
def chi(x1, x2,a):
    t1 = np.abs(x1) <= a/2
    t2 = np.abs(x2) <= a/2
    t3 = np.logical_and(t1, t2)
    p = np.where(t3 == True, 0, 1)
    return p

#fundamental solution to poisson
def phi(X1,X2):
    return np.log(np.sqrt(X1**2+X2**2)+1.0e-50)/(2*np.pi)

#fundamental solution derivative
def dphi(X1,X2):
    d1 = X1/(X1**2+X2**2+1.0e-50)
    d2 = X2/(X1**2+X2**2+1.0e-50)
    return np.array([d1/(2*np.pi),d2/(2*np.pi)])

#inverts the laplacian using sine series
def sseries_u(lap,sq_len,printflag=False):
    if printflag:
        print("starting invert laplacian")
    start = time.time()
    norm = (sq_len/np.pi)**2
    ak_mat = mymath.mydst2(lap)
    nrow, ncol = np.shape(ak_mat)
    k1, k2 = np.indices((nrow-1,ncol-1))
    bk_mat = np.zeros((nrow,ncol))
    bk_mat[1:,1:] = (-ak_mat[k1+1,k2+1]/((k1+1)**2+(k2+1)**2))
    u = mymath.mydst2(bk_mat)
    end = time.time()
    if printflag:
        print("ending invert laplacian, time=",end-start)
    return u*norm

#computes the laplacian using sine series
def sseries_lapu(u,sq_len,printflag=False):
    if printflag:
        print("starting compute laplacian")
    start = time.time()
    norm = (sq_len/np.pi)**2
    ak_mat = mymath.mydst2(u)
    nrow, ncol = np.shape(ak_mat)
    k1, k2 = np.indices((nrow - 1, ncol - 1))
    bk_mat = np.zeros((nrow, ncol))
    bk_mat[1:, 1:] = (-ak_mat[k1+1,k2+1]*((k1+1)**2+(k2+1)**2))
    u = mymath.mydst2(bk_mat)
    end = time.time()
    if printflag:
        print("ending compute laplacian, time=",end-start)
    return u/norm

#computes convolution of chi and phi on [-a/2,a/2]
#source is padded to match size of ker
def conv(src,ker,nsmdim,dx,printflag=False):
    if printflag:
        print("starting compute convolution")
        print("source size =",np.shape(src))
        print("kernel size=",np.shape(ker))
    start = time.time()
    nk = np.shape(ker)[0]
    ns = np.shape(src)[0]
    s = np.zeros((nk-1,nk-1))
    nsbeg = round((nk-ns)/2)
    nsend = nsbeg + ns
    s[nsbeg:nsend,nsbeg:nsend]=src
    sshft = sp.fft.fftshift(s)
    if printflag:
        print("FFT size =",np.shape(s))
    sshftfft = sp.fft.fft2(sshft)
    kfft = sp.fft.fft2(ker[1:,1:])
    c = sp.fft.ifft2(sshftfft*kfft)
    ncbeg = round((nk-nsmdim)/2)
    ncend = ncbeg + nsmdim
    c = np.real(c)[ncbeg-1:ncend-1,ncbeg-1:ncend-1]*dx**2
    end = time.time()
    if printflag:
        print("computed convolution, time=",end-start)
    return c

#fft based partials
def fft_partial(f,nsmdim,a):
    f_x1 = np.zeros(np.shape(f))
    f_x2 = np.zeros(np.shape(f))
    kx1 = (2*np.pi / a) * sp.fft.fftfreq(nsmdim - 1, 1. / (nsmdim - 1))
    KX1, KX2 = np.meshgrid(kx1, kx1)
    func = f[0:nsmdim-1,0:nsmdim-1]
    fshft = sp.fft.fftshift(func)
    fshift_hat = sp.fft.fft2(fshft)
    fshift_hat_x1 = 1j * KX1 * fshift_hat
    fshift_hat_x2 = 1j * KX2 * fshift_hat
    fshift_x1 = sp.fft.ifft2(fshift_hat_x1)
    fshift_x2 = sp.fft.ifft2(fshift_hat_x2)
    f_x1[0:nsmdim-1,0:nsmdim-1] = sp.fft.fftshift(np.real(fshift_x1))
    f_x2[0:nsmdim-1,0:nsmdim-1] = sp.fft.fftshift(np.real(fshift_x2))
    return f_x1, f_x2

#finite difference based partials
def fin_diff_partial(f,dx):
    nrow, ncol = np.shape(f)
    fx1 = np.zeros((nrow,ncol))
    fx2 = np.zeros((nrow,ncol))
    k1, k2 = np.indices((nrow-2, ncol-2))
    fx1[k1+1,k2+1] = (f[k1+1,k2+2]-f[k1+1,k2])/(2*dx)
    fx2[k1+1,k2+1] = (f[k1+2,k2+1]-f[k1,k2+1])/(2*dx)
    return fx1, fx2

#finite difference laplacian
def fin_diff_laplacian(f,dx):
    nrow, ncol = np.shape(f)
    lap_f = np.zeros((nrow-2,ncol-2))
    k1, k2 = np.indices((nrow - 2, ncol - 2))
    lap_f[k1,k2] = f[k1,k2+1]+f[k1+2,k2+1]+f[k1+1,k2]+f[k1+1,k2+2]-4*f[k1+1,k2+1]
    return lap_f/(dx**2)

#point source function
def W1(X1,X2,wts,shifts,const=1.):
    s = np.zeros(np.shape(X1))
    i=0
    for wt in wts:
        s += wt*phi(X1-shifts[i][0],X2-shifts[i][1])
        i += 1
    return -s/const

#gradient of point source function
def dW1(X1,X2,wts,shifts,const=1.):
    s = np.array([np.zeros(np.shape(X1)),np.zeros(np.shape(X2))])
    i=0
    for wt in wts:
        s += wt*dphi(X1-shifts[i][0],X2-shifts[i][1])
        i += 1
    return -s/const


def sigma(a,X1,X2,theta,wts,rot=False):
    q1center = np.array([(a / 6), (a / 6)])
    q2center = np.array([(-a / 6), (a / 6)])
    q3center = np.array([(-a / 6), (-a / 6)])
    q4center = np.array([(a / 6), (-a / 6)])
    z1 = mymath.t6hat(a / 6., a / 50., X1 - q1center[0]) * mymath.t6hat(a / 6., a / 50., X2 - q1center[1])
    z2 = mymath.t6hat(a / 6., a / 50., X1 - q2center[0]) * mymath.t6hat(a / 6., a / 50., X2 - q2center[1])
    z3 = mymath.t6hat(a / 6., a / 50., X1 - q3center[0]) * mymath.t6hat(a / 6., a / 50., X2 - q3center[1])
    z4 = mymath.t6hat(a / 6., a / 50., X1 - q4center[0]) * mymath.t6hat(a / 6., a / 50., X2 - q4center[1])
    sigma = wts[0]*z1 + wts[1]*z2 + wts[2]*z3 + wts[3]*z4
    if rot:
        return ndimage.rotate(np.exp(sigma),theta,reshape=False,mode='reflect')
    else:
        return np.exp(sigma)

#computes convolution for fundamental solution kernel
def poisson_convolution(X1bg,X2bg,X1md,X2md,source,a,m2sbeg,nsmdim,dx):
    z = mymath.t6hat(a, a / 2, X1md) * mymath.t6hat(a, a / 2, X2md)
    x = chi(X1md, X2md, a)
    k = phi(X1bg,X2bg)
    u0 = sseries_u(source,2*a)
    u0z = u0*z
    g = sseries_lapu(u0z, 2*a)
    src = x*g
    return u0z[m2sbeg:m2sbeg+nsmdim,m2sbeg:m2sbeg+nsmdim] - conv(src,k,nsmdim,dx)

#for an array of sources, outputs array of their convolutions with fundamental solution
#takes kernel, smoothing, and indicator functions as input for iteration schemes
def poisson_convolution2(z,x,k,sources,a,m2sbeg,nsmdim,dx,printflag=False):
    print("starting convolution with fundamental solution")
    out = np.array([np.zeros((nsmdim,nsmdim)),np.zeros((nsmdim,nsmdim))])
    i = 0
    for s in sources:
        u0 = sseries_u(s,2*a,printflag)
        u0z = u0*z
        g = sseries_lapu(u0z, 2*a,printflag)
        src = x*g
        out[i] = u0z[m2sbeg:m2sbeg+nsmdim,m2sbeg:m2sbeg+nsmdim] - conv(src,k,nsmdim,dx,printflag)
        i += 1
    return out

#this method generates a neighborhood of points around the point source location
def set_init_currents_at_source(location,gap,X1,X2,num_lines,min_theta,max_theta):
    print("generating initial values near source")
    x1_inits = np.array([])
    x2_inits = np.array([])
    loc_x1 = location[0]
    loc_x2 = location[1]
    theta = np.linspace(min_theta, max_theta, num_lines)
    for thet in theta:
        x1 = loc_x1+gap*np.cos(thet)
        x2 = loc_x2+gap*np.sin(thet)
        if x1>=np.min(X1) and x1<=np.max(X1) and x2>=np.min(X2) and x2<=np.max(X2):
            x1_inits = np.append(x1_inits,x1)
            x2_inits = np.append(x2_inits,x2)
    return np.array((x1_inits,x2_inits)).T

#this method sets n inits along the line perpendicular to segment connecting source locations
def set_init_currents_midline(location1,location2,num_inits,a):
    print("generating initial values between sources")
    x1_inits = np.linspace(-a/2.,a/2.,num_inits)
    x11 = location1[0]
    x21 = location1[1]
    x12 = location2[0]
    x22 = location2[1]
    m = (x22-x21)/(x12-x11)
    slope = -1./m
    x2_inits = slope*x1_inits
    return np.array((x1_inits, x2_inits)).T

def set_inits(X1,X2,src1,src2,currx1,currx2,n_lines,a,n_disc):
    print("getting initial values")
    init_arr = []
    partial_sums,inputs = total_flux(X1,X2,src1,src2,currx1,currx2,a,n_disc)
    C = partial_sums[-1]
    print("total flux = ", C)
    partition = C*np.linspace(0.,1.,n_lines)
    for part in partition:
        idx = np.abs(partial_sums-part).argmin()
        init_arr.append(inputs[idx])
    print("initial values = ", np.array(init_arr))
    return np.array(init_arr)

def total_flux(X1,X2,src1,src2,currx1,currx2,a,n_disc):
    x = np.linspace(-a / 2., a / 2., n_disc)
    x11 = src1[0]
    x21 = src1[1]
    x12 = src2[0]
    x22 = src2[1]
    m = (x22 - x21) / (x12 - x11)
    slope = -1. / m
    y = slope * x
    inputs = np.array((x,y)).T
    normal_vec = np.array([1./np.sqrt(2.),1./np.sqrt(2.)])
    line_int = 0.
    partial_sum_arr = []
    for input in inputs:
        I11, I12, I21, I22, x11, x12, x21, x22 = get_corners(X1, X2, input[0], input[1])
        currx = bilinear_interp(I11,I12,I21,I22,x11,x12,x21,x22,currx1,input)
        curry = bilinear_interp(I11,I12,I21,I22,x11,x12,x21,x22,currx2,input)
        J = np.array([currx,curry])
        line_int += np.dot(J,normal_vec)
        partial_sum_arr.append(line_int)
    return np.array(partial_sum_arr), inputs

# bilinear interp for smooth potential, analytic for source potential
def get_grad_potentail(init,X1,X2,potX1,potX2,wts,shifts,a,theta,sig_wts):
    grad_source_pot_x1, grad_source_pot_x2 = dW1(init[0],init[1],wts,shifts)
    I11,I12,I21,I22,x11,x12,x21,x22 = get_corners(X1,X2,init[0],init[1])
    grad_smooth_pot_x1 = bilinear_interp(I11,I12,I21,I22,x11,x12,x21,x22,potX1,init)
    grad_smooth_pot_x2 = bilinear_interp(I11,I12,I21,I22,x11,x12,x21,x22,potX2,init)
    curr_x1 = sigma(a,init[0],init[1],theta,sig_wts)*(grad_source_pot_x1 + grad_smooth_pot_x1)
    curr_x2 = sigma(a,init[0],init[1],theta,sig_wts)*(grad_source_pot_x2 + grad_smooth_pot_x2)
    return np.array([curr_x1,curr_x2])


def bilinear_interp(I11,I12,I21,I22,x11,x12,x21,x22,func,init):
    fq11 = func[I11[0],I11[1]]
    fq21 = func[I21[0],I21[1]]
    fq12 = func[I12[0],I12[1]]
    fq22 = func[I22[0],I22[1]]
    x1 = init[0]
    x2 = init[1]
    fx1_x21 = ((x12-x1)/(x12-x11))*fq11 + ((x1-x11)/(x12-x11))*fq21
    fx1_x22 = ((x12-x1)/(x12-x11))*fq12 + ((x1-x11)/(x12-x11))*fq22
    return ((x22-x2)/(x22-x21))*fx1_x21 + ((x2-x21)/(x22-x21))*fx1_x22


def get_corners(X1,X2,x1,x2):
    X1proj1 = X1[0,:]
    X2proj1 = X2[:,0]
    distarrX1 = np.sqrt((X1proj1-x1)**2)
    distarrX2 = np.sqrt((X2proj1-x2)**2)
    X1inds = np.sort(np.argpartition(distarrX1,2,axis=None)[:2])
    X2inds = np.sort(np.argpartition(distarrX2,2,axis=None)[:2])
    I11 = np.array([X1inds[0],X2inds[0]])
    I12 = np.array([X1inds[0],X2inds[1]])
    I21 = np.array([X1inds[1],X2inds[0]])
    I22 = np.array([X1inds[1],X2inds[1]])
    x11 = X1proj1[X1inds[0]]
    x12 = X1proj1[X1inds[1]]
    x21 = X2proj1[X2inds[0]]
    x22 = X2proj1[X2inds[1]]
    return I11,I12,I21,I22,x11,x12,x21,x22



def md_ptv1(g,y,X1,X2,potX1,potX2,wts,src,a,src_tol,shifts,theta,sig_wts,neg=False,max_its=10000.):
    print("Starting Mid point")
    ctr = 0
    x1min = X1[0][0]
    x1max = X1[0][-1]
    x2min = X2[0][0]
    x2max = X2[-1][0]
    x1_arr = [y[0]]
    x2_arr = [y[1]]
    while y[0] >= x1min and y[0] <= x1max and y[1] >= x2min and y[1] <= x2max and np.linalg.norm(y-src) >= src_tol and ctr <= max_its:
        if np.linalg.norm(y-src)<=10.*src_tol:
            h = src_tol/100.
        else:
            h = src_tol
        k1 = h*g(y,X1,X2,potX1,potX2,wts,shifts,a,theta,sig_wts)
        nexty = h*g(y+.5*k1,X1,X2,potX1,potX2,wts,shifts,a,theta,sig_wts)
        if neg:
            y -= nexty
        else:
            y += nexty

        if np.linalg.norm(y)<1.e-8:
            print(ctr)
            return x1_arr, x2_arr
        elif np.linalg.norm(y-src)<=src_tol:
            print(ctr)
            return x1_arr, x2_arr
        else:
            ctr += 1
            x1_arr.append(y[0])
            x2_arr.append(y[1])
    print(ctr)
    return x1_arr, x2_arr


def md_ptv2(g,y,X1,X2,potX1,potX2,wts,src,a,src_tol,shifts,theta,sig_wts,neg=False,max_its=5000.):
    print("Starting Mid point")
    ctr = 0
    x1min = X1[0][0]
    x1max = X1[0][-1]
    x2min = X2[0][0]
    x2max = X2[-1][0]
    x1_arr = [y[0]]
    x2_arr = [y[1]]
    while y[0] >= x1min and y[0] <= x1max and y[1] >= x2min and y[1] <= x2max and ctr <= max_its:
        if np.linalg.norm(y-src)<=10.*src_tol:
            h = src_tol/100.
        else:
            h = src_tol
        k1 = h*g(y,X1,X2,potX1,potX2,wts,shifts,a,theta,sig_wts)
        nexty = h*g(y+.5*k1,X1,X2,potX1,potX2,wts,shifts,a,theta,sig_wts)
        if neg:
            y -= nexty
        else:
            y += nexty

        if np.linalg.norm(y)<1.e-8:
            print(ctr)
            return x1_arr, x2_arr
        elif np.linalg.norm(y-src)<=src_tol:
            print(ctr)
            return x1_arr, x2_arr
        else:
            ctr += 1
            x1_arr.append(y[0])
            x2_arr.append(y[1])
    print(ctr)
    return x1_arr, x2_arr










