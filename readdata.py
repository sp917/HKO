import numpy as np
from netCDF4 import Dataset
import math
import matplotlib.pyplot as plt
import sys
import os

def main():


    var = sys.argv[1]
    year = '%0.4d' % int(sys.argv[2])
    month = '%0.2d' % int(sys.argv[3])
    day = '%0.2d' % int(sys.argv[4])
    hour = '%0.2d' % int(sys.argv[5])
    t = '%0.2d' % int(sys.argv[6])

    date = year + '-' + day + '-' + month + '_' + hour + ':00:00.t' + t
    filename = 'UV_interp_' + date + '.nc'
    path = '/home/sp917/data/'
    
    print("Reading data for variable " + var + " from file " + filename)
    f = Dataset(path+filename,'r', format='NETCDF4_CLASSIC')

    vardata = f.variables[var]
    X = np.array(vardata)

    print(X.shape) # the shape of x should be time*height*lat*lon* or time*lat*lon

    X = X[0,]
    if len(X.shape) ==3:
        X = np.transpose(X,(1,2,0))
    elif len(X.shape)==2:
        print("")
    else:
        raise RuntimeError("Shape of X not understood")

    #the size of the time dimension of x and y is 1
    x = np.array(f.variables['XLONG'])[0,:,:]
    y = np.array(f.variables['XLAT'])[0,:,:]
    
    x,y = deg2km(x,y)

    dy = delta(y)
   
    x = np.transpose(x)
    dx = delta(x)
    dx = np.transpose(dx)
    x = np.transpose(x)

    print("Removing linear trend from data...")

    X = detrend(X)

    print("Done.")

    print("Performing fast Fourier transform on de-trended data...")

    Y = np.fft.fftn(X, axes=(0,1))
    
    print("Done.")

    print("Calculating spectrum...")

    S,K = errico(Y,dx,dy)

    if len(S.shape)==2:
        S = np.transpose(S)
    elif len(S.shape)==1:
        print("")
    else:
        raise RuntimeError("Something went wrong!")
    
    print("Done.")

    print("Saving spectrum to netCDF file...")
  
    newfile_name = "Spectrum_"+ var + "_" +  date + ".nc"
   
    deletefile(path+newfile_name)

    print("Creating file " + path+newfile_name)
    newfile = Dataset(path+newfile_name,'w',format='NETCDF4_CLASSIC')

    print(S.shape)

    if len(S.shape)==2:
        newfile.createDimension('wavenumber', K.shape[0])
        newfile.createDimension('bottom_top', S.shape[0])
        spectrum = newfile.createVariable('spectrum', np.float64, ('bottom_top', 'wavenumber'))
        spectrum[:,:] = S
    elif len(S.shape)==1:
        newfile.createDimension('wavenumber', K.shape[0])
        spectrum = newfile.createVariable('spectrum', np.float64, ('wavenumber'))
        spectrum[:] = S

    wavenumber = newfile.createVariable('k', np.float64, 'wavenumber')
    wavenumber[:] = K

    print("Done.")
    
    return

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

main()

  
