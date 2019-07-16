import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import sys
import os


def filt1(X, yvals, xvals, ny, nx):

    """ Spatially filters the first two dimensions of array X into windows of size nx*ny. I.E. nx+1 and ny+1 are the number of grid points contained in each window, in each direction."""
    
    ylen = X.shape[0]
    xlen = X.shape[1]

    yflen = (ylen-1)//ny
    xflen = (xlen-1)//nx

    Y = np.zeros((X.shape))

    #Y = Y[0:yflen,0:xflen,]

    ymax = ny*yflen+1
    xmax = nx*xflen+1

    Y = Y[0:ymax,0:xmax,]
    Xnew = X[0:ymax,0:xmax,]
    yvals = yvals[0:ymax,0:xmax,]
    xvals = xvals[0:ymax,0:xmax,] 

    counter = np.zeros((Y.shape))
    
    for i in range(xflen):
        xmin = nx*i
        xmax = nx*(i+1)+1
        for j in range(yflen):
            ymin = ny*j
            ymax = ny*(j + 1)+1
            #print((xmin,xmax), (ymin,ymax))
            Y[ymin:ymax,xmin:xmax,] = Y[ymin:ymax,xmin:xmax,] + np.mean(X[ymin:ymax,xmin:xmax,], axis=(0,1))
            counter[ymin:ymax,xmin:xmax,] = counter[ymin:ymax,xmin:xmax,] + 1

    Y = Y/counter #We take the average of the points that appear more than once

    return Xnew, Y, yvals, xvals

def filt2(X, yvals, xvals, ny, nx):
    """Different type of spatial filter which returns the average over the neighbours within a window of size nx*ny"""

    Y = dofilter2(X,nx,ny)
   
    Xnew = dofilter2(X,nx%2,ny%2)
    xvalsnew = dofilter2(xvals,ny%2,nx%2)
    yvalsnew = dofilter2(yvals,ny%2,nx%2)

    return Xnew, Y, yvalsnew, xvalsnew

def dofilter2(X,ny,nx):
    
    Y = np.zeros((X.shape))
    Y = Y[(ny%2):,(nx%2):,] #If ny or nx is odd then we need to shift the grid on which Y is defined

    print("X.shape=",X.shape)
    print("Y.shape=",Y.shape)
    
    for i in range(Y.shape[1]):
        if (i%(Y.shape[1]//10)==0): print(str(100*i/Y.shape[1]) + '%')
        if nx%2==0:
            xmin = max(0, i-nx//2)
            xmax = min(Y.shape[1], i+1+nx//2)
        else:
            xmin = max(0, i-(nx-1)//2)
            xmax = min(Y.shape[1], i+1+(nx+1)//2)
        for j in range(Y.shape[0]):
            if ny%2==0:
                ymin = max(0, j-ny//2)
                ymax = min(Y.shape[0], j+1+ny//2)
            else:
                ymin = j-(ny-1)//2
                ymax = j+1 + (ny+1)//2 
                ymin = max(0, ymin)
                ymax = min(Y.shape[0], ymax ) 
            Y[j,i,] = np.mean( X[ymin:ymax,xmin:xmax,] , axis = (0,1) )
    
    return Y

def project(X,ny,nx):

    if (ny%2==0): 
        Xnew = X[ny//2: -ny//2,]
        ny2 = ny//2
    else:
        ny2=(ny-1)//2
    
    if (nx%2==0):
        Xnew = X[ny//2: -ny//2,]
        nx2 = nx//2
    else:
        nx2 = (nx-1)//2

    Xnew = dofilter2(X,ny%2,nx%2)
    
    if ny!=1:
        Xnew = Xnew[ny2:-ny2,]
    if nx!=1:
        Xnew = Xnew[:,nx2:-nx2,]
    
    return Xnew

def detrend(X):

    nym1 = X.shape[0]-1
    nxm1 = X.shape[1]-1
    
    Xdetrend = np.zeros((X.shape))
    Xdetrend = Xdetrend[1:,1:,]

    for j in range(0,nxm1):
        for i in range(0,nym1):
            Xdetrend[i,j,] = X[i,j,] - (X[nym1,j,]-X[0,j,])*i/nym1 - (X[i,nxm1,]-X[i,0,])*j/nxm1 + \
                    (X[nym1,nxm1,] - X[nym1,0,] - X[0,nxm1,] + X[0,0,])*i*j/(nym1*nxm1)

    return(Xdetrend)


def errico(Y,dx,dy):

    #See Errico 1985, page 1555

    nxm1 = Y.shape[1]
    nym1 = Y.shape[0]
    if len(Y.shape)==3:
        nz = Y.shape[2]
    else:
        print("")

    DX = np.mean(dx, axis=(0,1))
    DY = np.mean(dy, axis=(0,1))

    ap = DX*nxm1
    aq = DY*nym1

    A = max(ap,aq)

    if abs(A-ap)<0.00000000001:
        kmax = (nym1+1)//2
    else:
        kmax = (nxm1+1)//2

    try:
        S = np.zeros((kmax,nz))
    except:
        S = np.zeros((kmax))

    Yrescaled = Y/(nxm1*nym1) #This scaling just keeps the numbers from getting too large
    
    for k in range(0,kmax):
        if k%(kmax//10)==0: print(str(100*k/kmax) + '%')
        qmax = math.ceil(min(nxm1,aq*(k+0.5)/A))
        p0 = ap*(k-0.5)/A
        q0 = aq*(k-0.5)/A
        for q in range(0,qmax):
            pmin = ((k-0.5)/A)**2- (q/aq)**2
            if pmin < 0:
                pmin = 0    
            else:
                pmin = ap*np.sqrt(pmin)
            pmax = ap*np.sqrt(((k+0.5)/A)**2- (q/aq)**2)
            pmin = max(math.floor(pmin),0)
            pmax = min(math.ceil(pmax),nym1)
            for p in range(pmin,pmax):
                if abs( A*np.sqrt( (p/ap)**2 + (q/aq)**2 ) - k ) < 0.5:
                    S[k,] = S[k,] + np.absolute(Yrescaled[p,q,])
    
    K = np.array(range(0,kmax))

    K = 2*np.pi*K/A
    
    return S, K

def delta(y):

    dy = np.zeros((y.shape))
    dy=dy[1:,] #the first dimension should be one less than that of y

    for n in range(0,dy.shape[0]):
        dy[n,] = y[n+1,] - y[n,]

    return dy

def deg2km(x,y):
    R = 6371
    pi = np.pi

    xnew = x*R*np.cos(y*pi/180)*pi/180
    ynew = y*R*pi/180

    return xnew, ynew


def deletefile(filename_full):
    try:
        f = open(filename_full, 'r')
        f.close()
    except:
        print("File " + filename_full + " does not exist")
    else:
        print('Removing file ' + filename_full + '...')
        try:
            os.remove(filename_full)
            print('File removed.')
        except:
            print("Could not remove file " + filename_full)


def gettimes(datapath,date):
    prefix = "wrfout_xhka_"
    times = []
    for f in os.listdir(datapath):
        if (f.startswith(prefix+date + ".t") and f.endswith(".nc")):
            t = f.replace(prefix+date+".t","")
            t = t.replace(".nc","")
            times = times + [t]

    return times

