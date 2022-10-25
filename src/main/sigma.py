import numpy as np
import mymath
from scipy import ndimage


#generate an initial conductivity
def gen_sigma_1(a,X1,X2,wts):
    q1center = np.array([(a / 6), (a / 6)])
    q2center = np.array([(-a / 6), (a / 6)])
    q3center = np.array([(-a / 6), (-a / 6)])
    q4center = np.array([(a / 6), (-a / 6)])
    z1 = mymath.t6hat(a / 6., a / 50., X1 - q1center[0]) * mymath.t6hat(a / 6., a / 50., X2 - q1center[1])
    z2 = mymath.t6hat(a / 6., a / 50., X1 - q2center[0]) * mymath.t6hat(a / 6., a / 50., X2 - q2center[1])
    z3 = mymath.t6hat(a / 6., a / 50., X1 - q3center[0]) * mymath.t6hat(a / 6., a / 50., X2 - q3center[1])
    z4 = mymath.t6hat(a / 6., a / 50., X1 - q4center[0]) * mymath.t6hat(a / 6., a / 50., X2 - q4center[1])
    sigma = wts[0]*z1 + wts[1]*z2 + wts[2]*z3 + wts[3]*z4
    return np.exp(sigma)

#rotate conductibity
def rot_sigma(theta,sigma):
    return ndimage.rotate(np.exp(sigma),theta,reshape=False,mode='reflect')
