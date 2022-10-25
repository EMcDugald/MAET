import numpy as np
import math
from scipy.fftpack import dst

def t6hat(rmax,rmin,x):

   xx   = np.abs(x)
   ones = np.ones_like(xx)
   xmax = rmax*ones
   xmin = rmin*ones

   xx1 = np.minimum(xx,xmax)
   xx2 = np.maximum(xx1,xmin)
   xx3 = (ones - (xx2 - rmin)/(rmax-rmin))*math.pi
   xx4 = sin6hat(xx3)
   return xx4

def sin6hat(x):
    res = np.sin(6.0*x)/3.0 - 3.0*np.sin(4.0*x)   + 15.0*np.sin(2.0*x) - 20.0*x
    res = -res/(20.0*math.pi)
    return res

def doublefilt(zerocut,onecut,lowzerocut,lowonecut,input):
    ntotal = len(input)

    z =  np.fft.fftshift(np.fft.fft(input))
    fr = np.zeros_like(input)

    for jj in range(ntotal):
        fr[jj] = (jj - ntotal/2 )/ (ntotal/2)

    hat  = t6hat(zerocut,onecut,fr)

    # z = hat*z


    if(lowonecut > 0.0001):
        hat1 = np.ones_like(fr)
        hat1 = hat1 - t6hat(lowonecut,lowzerocut,fr)
        hat = hat*hat1

    z = hat*z

    #plt.plot(hat1*hat)
    #plt.show()

    z =  np.fft.fftshift(z)
    out = np.fft.ifft(z)
    out = np.real(out)

    return out, hat

def triplefilt(zerocut,onecut,lowzerocut,lowonecut,input,shapefft):
    ntotal = len(input)

    z =  np.fft.fftshift(np.fft.fft(input))

    z = z*shapefft

    fr = np.zeros_like(input)

    for jj in range(ntotal):
        fr[jj] = (jj - ntotal/2 )/ (ntotal/2)

    hat  = t6hat(zerocut,onecut,fr)

    # z = hat*z


    if(lowonecut > 0.0001):
        hat1 = np.ones_like(fr)
        hat1 = hat1 - t6hat(lowonecut,lowzerocut,fr)
        hat = hat*hat1

    z = hat*z

    #plt.plot(hat1*hat)
    #plt.show()

    z =  np.fft.fftshift(z)
    out = np.fft.ifft(z)
    out = np.real(out)

    return out, hat


def mydst2(input):

    nd1, nd2 = np.shape(input)
    temp   = dst(input[1:nd1-1,1:nd1-1], type=1, axis = 0, norm='ortho')

    ndd1, ndd2 = np.shape(temp)
    temp1  = dst(temp,  type=1, axis = 1, norm='ortho')

    print('MyDST2: ',ndd1,ndd2)

    output = np.zeros((nd1,nd2))
    output[1:nd1-1,1:nd2-1] = temp1

    return output


