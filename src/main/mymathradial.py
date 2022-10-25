
import numpy as np


#-----------------------------------------------------------------------
# This function returns a radially symmetric function centered at (rcx,rcy)
# with radiii rmax,rmin,   and its exactly computed laplacian.
# These functions are computed on a grid (ndim x ndim) in the square
# [-halfsize,halfsize] x [-halfsize,halfsize]
#------------------------------------------------------------------------
def radial(ndim,rcx,rcy,rmax,rmin,halfsize):
   step = 2.0*halfsize/(ndim-1)

   #func   = np.zeros((ndim,ndim))
   #exlap  = np.zeros((ndim,ndim))

   xx = np.linspace(-halfsize,halfsize,ndim)
   yy = xx
   [x,y] = np.meshgrid(xx,yy)

   rad = np.sqrt((x-rcx)*(x-rcx)+(y-rcy)*(y-rcy)) + 1.e-50

   [hat,der,der2] = hat_rad_derivs(rmax,rmin,rad)
   der = -der
   der2 = -der2

   rloga = np.log(rad)
   exlapla = - (rloga*der2 + (2.0+rloga)*der/rad)
   func = -(1.-hat)*rloga

   return [func,exlapla]



#----------------------------------------------------------------------
# This function returns a hat function and its first and second derivatives
# The hat function equals 1 on the interval [-rmin,rmin] and it vanishes
# outside of the interval [-rmax,rmax]. Smooth C6 transition in between.
#------------------------------------------------------------------------
def hat_rad_derivs(rmax,rmin,rad):

    arg = (rmax-rad)/(rmax-rmin)*np.pi
    rkoef =  -1.0/(20.0*np.pi)

    hat = rkoef*(  np.sin(6.0*arg)/3.0 - 3.0*np.sin(4.0*arg)  + 15.0*np.sin(2.0*arg) - 20.0*arg )

    der = - (np.sin(arg)**6) /(rmax-rmin)*3.20

    der2 =6*(np.sin(arg)**5)*np.cos(arg) /(rmax-rmin)/(rmax-rmin)*3.20*np.pi

    hat[rad > rmax] = 0.
    der[rad > rmax] = 0.
    der2[rad > rmax] = 0.

    hat[rad < rmin] = 1.
    der[rad < rmin] = 0.
    der2[rad < rmin] = 0.

    return [hat,der,der2]

#-----------------------------------------------------------------------
# This function returns a radially symmetric function centered at (rcx,rcy)
# with radiii rmax,rmin,   and its exactly computed laplacian.
# These functions are computed on a grid (ndim x ndim) in the square
# [-halfsize,halfsize] x [-halfsize,halfsize]
#------------------------------------------------------------------------
def radial_old(ndim,rcx,rcy,rmax,rmin,halfsize):
   step = 2.0*halfsize/(ndim-1)

   func   = np.zeros((ndim,ndim))
   exlap  = np.zeros((ndim,ndim))

   for i in range(0,ndim):
      y = -halfsize + i*step
      for j in range(0,ndim):
         x = -halfsize + j*step
         rad = np.sqrt((x-rcx)*(x-rcx)+(y-rcy)*(y-rcy))

         if(rad > 1.e-30):

           [hat,der,der2] = hatderivs(rmax,rmin,rad)
           der = -der
           der2 = -der2

           rloga = np.log(rad)
           exlapla = - (rloga*der2 + (2.0+rloga)*der/rad)
           exlap[i,j] = exlapla
           func[i,j] = -(1.-hat)*rloga

   return [func,exlap]

#----------------------------------------------------------------------
# This function returns a hat function and its first and second derivatives
# The hat function equals 1 on the interval [-rmin,rmin] and it vanishes
# outside of the interval [-rmax,rmax]. Smooth C6 transition in between.
#------------------------------------------------------------------------
def hatderivs(rmax,rmin,xsign):
    signum = 1.0
    x = np.abs(xsign)
    if(xsign < 0.0):
       signum = -1.0

    hat = 0.0
    der = 0.0
    der2 = 0.0

    if(x > rmax):
        return hat,der,der2

    if(x < rmin):
        hat = 1.0
        return hat,der,der2

    arg = (rmax-x)/(rmax-rmin)*np.pi
    rkoef =  -1.0/(20.0*np.pi)

    hat = rkoef*(  np.sin(6.0*arg)/3.0 - 3.0*np.sin(4.0*arg)  + 15.0*np.sin(2.0*arg) - 20.0*arg )

    der = - signum*(np.sin(arg)**6) /(rmax-rmin)*3.20

    der2 =6*(np.sin(arg)**5)*np.cos(arg) /(rmax-rmin)/(rmax-rmin)*3.20*np.pi

    return [hat,der,der2]



#-----------------------------------------------------------------