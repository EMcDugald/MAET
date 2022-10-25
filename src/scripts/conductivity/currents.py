import numpy as np
import my_conductivity as mc
import mymath
import utils
import time
from scipy.sparse.linalg import gmres
from scipy.sparse.linalg import LinearOperator
import matplotlib.pyplot as plt
import configparser


#--- This script generates current lines from conductivity ---#
#--- A conductivity function is called from main module    ---#
#--- and potential is obtained as a direct solution to     ---#
#--- conductivity equation. Currents are derived from      ---#
#--- the potential, and drawn. This is repeated subject to ---#
#--- rotation of the phantom. Solution of conductivity     ---#
#--- is obtained via gmres.                                ---#


#-----------Functions for GMRES -----------------------------------------------#
def my_fredholm(rho,printflag): #generates LHS of fredholm operator#
    print("generating fredholm operator")
    Rho = np.zeros((nmddim,nmddim))
    Rho[m2sbeg:m2sbeg+nsmdim,m2sbeg:m2sbeg+nsmdim] = rho #--place rho in medium square for convolution--#
    rho_x1,rho_x2 = mc.fft_partial(Rho,nmddim,4.)  #--gradient of rho--#
    conv_x1,conv_x2 = mc.poisson_convolution2(z,x,k,np.array([rho_x1,rho_x2]),2.,m2sbeg,nsmdim,dx,printflag)  #--convolve grad(rho) with fundamental solution--#
    conv_prod = log_cond_x1*conv_x1 + log_cond_x2*conv_x2  #--compute product of convolution with gradient log conductivity--#
    mat = Rho[m2sbeg:m2sbeg+nsmdim,m2sbeg:m2sbeg+nsmdim] + conv_prod   #--operator for gmres--#
    return mat


def my_lin_op(rho_vec): #reshape output of fredholm LHS#
    global i
    printflag = False
    if i == 1:
        printflag = True
    rho = rho_vec.reshape((nsmdim,nsmdim)).T #--reshape rho_vec before calling methods--#
    out = np.ravel(my_fredholm(rho,printflag),order='F') #--LinearOperator must output a vector--#
    i+=1
    return out
#-----------------------------------------------------------------------------#


config = configparser.RawConfigParser()
config.read('currents.cfg')


# ------Parameters to set dimension of grids ------------------------#
n = config.getint('grid_dims','ndim')
nsmdim = 2**n+1
nmddim = 2*nsmdim-1
nbgdim = round( (nmddim-1)*(3/2))+1
print('nbgdim '+str(nbgdim))
print('nsmdim '+str(nsmdim))
print('nmddim '+str(nmddim))
m2sbeg = round((nmddim-nsmdim)/2)
# ------------------------------------------------------------------#

# --------Meshgrids for small, medium, and large squares ----------#
a= config.getfloat('grid_dims','sq_len')
xsm1 = np.linspace(-a/2.,a/2.,nsmdim)
xsm2 = np.linspace(-a/2.,a/2.,nsmdim)
Xsm1, Xsm2 = np.meshgrid(xsm1,xsm2) #small square#
xmd1 = np.linspace(-a,a,nmddim)
xmd2 = np.linspace(-a,a,nmddim)
Xmd1, Xmd2 = np.meshgrid(xmd1,xmd2) #domain#
xbg1 = np.linspace(-1.5*a,1.5*a,nbgdim)
xbg2 = np.linspace(-1.5*a,1.5*a,nbgdim)
Xbg1, Xbg2 = np.meshgrid(xbg1,xbg2) #large square#
dx = Xbg1[0][1] - Xbg1[0][0]


#---- Fixed Functions ---------------------------------------------#
wts = np.array([config.getfloat('src_wts', 'w1'),
                config.getfloat('src_wts', 'w2')])  # get weights for point sources
shifts = np.array([[config.getfloat('src_shifts', 'x11'),
                    config.getfloat('src_shifts', 'x21')],
                   [config.getfloat('src_shifts', 'x12'),
                    config.getfloat('src_shifts', 'x22')]])  # get shifts for point sources
point_source_potential = mc.W1(Xsm1, Xsm2, wts, shifts)  # generate point source potential#
utils.ndwrite(point_source_potential,"misc/point_source_potential.d")
z = mymath.t6hat(a, a / 2., Xmd1) * mymath.t6hat(a, a / 2., Xmd2)  # generate smooth transition function#
x = mc.chi(Xmd1, Xmd2, a)  # generate curring function, 0 in domain, 1 outside domain
k = mc.phi(Xbg1, Xbg2)  # generate fundamental solution of poisson's equation#
w_x1, w_x2 = mc.dW1(Xsm1, Xsm2, wts, shifts)  # generate gradient of point sources#
#------------------------------------------------------------------#


rot_ct = 1
nrots = config.getint('rotation','nrot')
sigma_wts = np.array([config.getfloat('sigma_wts','w1'),
                      config.getfloat('sigma_wts','w2'),
                      config.getfloat('sigma_wts','w3'),
                      config.getfloat('sigma_wts','w4')])
# thetas = np.linspace(0.,360.,nrots)
thetas = np.arange(0.,360.,360/nrots)
for theta in thetas:
    print("Rotation #: ",rot_ct)
    i = 1
    cond = mc.sigma(a,Xsm1,Xsm2,theta,sigma_wts,rot=True) #generate conductivity given theta from main module
    utils.ndwrite(cond,"sigma/sigma_"+f"{rot_ct:0>4}" + ".d")
    log_cond = np.log(cond+1.0e-50) #generate logarithm of conductivity#
    utils.ndwrite(log_cond,"sigma/log_sigma_"+f"{rot_ct:0>4}" + ".d")
    log_cond_x1, log_cond_x2 = mc.fft_partial(log_cond,nsmdim,2.) #generate gradient of conductivity#


    start = time.time()
    b = np.ravel(-(log_cond_x1*w_x1+log_cond_x2*w_x2),order='F')
    rho0 = b #initial guess for gmres is RHS of equation 5#
    utils.ndwrite(b.reshape((nsmdim,nsmdim)).T,"error/RHS_"+f"{rot_ct:0>4}"+".d") #write RHS for testing purposes#
    AA = LinearOperator((nsmdim**2,nsmdim**2),matvec=my_lin_op) #define operator for gmres#
    lap_pot, exitCode = gmres(AA,b,x0=rho0,maxiter=1,restart=40,tol = 1.e-8) #call gmres#
    print('Exit code = ',exitCode)
    lap_potential = lap_pot.reshape((nsmdim,nsmdim)).T #gmres returns laplacian of potential#


    #convolve laplacian of potential with fundamental solution to recover potential#
    print("Recovering Potential from Laplacian")
    lap_pot_pad = np.zeros((nmddim,nmddim))
    lap_pot_pad[m2sbeg:m2sbeg+nsmdim,m2sbeg:m2sbeg+nsmdim] = lap_potential
    utils.ndwrite(lap_potential,"misc/lap_potential_"f"{rot_ct:0>4}"+".d")
    smooth_potential = mc.poisson_convolution(Xbg1,Xbg2,Xmd1,Xmd2,lap_pot_pad,2.,m2sbeg,nsmdim,dx)
    utils.ndwrite(smooth_potential,"smooth_potential/potential_"+f"{rot_ct:0>4}"+".d")
    end = time.time()
    print("gmres ended, total time=",end-start)


    print("computing potential gradient ",rot_ct)
    pot_x1,pot_x2 = mc.fin_diff_partial(smooth_potential,dx) #compute potential of gradient#
    LHS = lap_potential + log_cond_x1*pot_x1+log_cond_x2*pot_x2 #compute LHS#
    utils.ndwrite(LHS,"error/LHS_"+f"{rot_ct:0>4}"+".d") #write LHS for testing purposes#

    print("computing current ",rot_ct)
    curr_x1 = cond*(pot_x1+w_x1) #current = conductivity*(grad(smooth potential + point source potential))
    curr_x2 = cond*(pot_x2+w_x2)
    utils.ndwrite(curr_x1,"currents/currx1_"+f"{rot_ct:0>4}"+".d")
    utils.ndwrite(curr_x2,"currents/currx2_"+f"{rot_ct:0>4}"+".d")
    print("computing curl ",rot_ct)
    K1, K2 = np.indices(np.shape(curr_x1))
    curl = curr_x1[K1,K2]*log_cond_x2[K1,K2]-curr_x2[K1,K2]*log_cond_x1[K1,K2] #compute curl
    utils.ndwrite(curl,"curl/curl_"+f"{rot_ct:0>4}"+".d")

    start = time.time()
    src_dist = config.getfloat('curr_lines','src_dist') #defines where to terminate current line based on proximity to source#
    num_lines = config.getint('curr_lines','nlines')
    source1 = shifts[0]  # location of point source
    source2 = shifts[1]  # location of point source

    n_disc = 1000 #number of discretization points to use for line integral
    inits = mc.set_inits(Xsm1,Xsm2,source1,source2,curr_x1,curr_x2,num_lines,a,n_disc)
    print("drawing current lines")
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5, 5))

    # plot currents from middle to source 1
    j1 = 1
    print("plotting currents from middle to source 1")
    for init in inits:
        print(j1)
        x1_vals1, x2_vals1 = mc.md_ptv2(mc.get_grad_potentail, init, Xsm1, Xsm2, pot_x1, pot_x2, wts, source1, a,
                                      src_dist, shifts, theta, sigma_wts,neg=True)
        ax.plot(x1_vals1, x2_vals1, color='k', linewidth=.5)
        j1 += 1

    print("**************************************")

    inits = mc.set_inits(Xsm1, Xsm2, source1, source2, curr_x1, curr_x2, num_lines, a, n_disc)
    #plot currents from middle to sources 2
    j2 = 1
    print("plotting currents from middle to source 2")
    for init in inits:
        print(j2)
        x1_vals2, x2_vals2 = mc.md_ptv2(mc.get_grad_potentail, init, Xsm1, Xsm2, pot_x1, pot_x2, wts, source2, a,
                                      src_dist, shifts, theta, sigma_wts,neg=False)
        ax.plot(x1_vals2, x2_vals2, color='k', linewidth=.5)
        j2 += 1

    ax.scatter(source1[0], source1[1], color='r')
    ax.scatter(source2[0], source2[1], color='b')

    ax.set_xlim((-a / 2., a / 2.))
    ax.set_ylim((-a / 2., a / 2.))
    plt.savefig("current_lines/center_midpoint_"+f"{rot_ct:0>4}"+".png", bbox_inches='tight')
    end = time.time()
    print("time to draw ", num_lines, " current lines:", end - start)
    rot_ct += 1