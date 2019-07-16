import sys
import os
import wrf
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick


def getdata(str_id):

    """ Reads in data and interoplates onto fixed horizontal levels """

    datapath = '/home/sp917/data/'
    files = [Dataset(datapath+f) for f in os.listdir(datapath) if f.startswith('wrfout_xhka')]    

    X = wrf.getvar(files, str_id, timeidx=wrf.ALL_TIMES, method='join', meta=True)
    z = wrf.getvar(files, 'z', timeidx=wrf.ALL_TIMES, method='join', meta=False)

    zmax = np.max(z)
    nz = z.shape[1]
    dz = zmax/(nz-1)
    levels=np.arange(0,zmax+dz,dz)

    try:
        X = wrf.vinterp(files, field=X, vert_coord='ght_msl', interp_levels=levels, field_type='z', timeidx=wrf.ALL_TIMES)
        print('Interpolating onto vertical levels')
    except:
        print('')

    return X

def plotsmooth(X, Xsmooth, yvals, xvals, plottitle, plotname, plotpath):

    """ X: numpy array of dimensions (nt,ny,nx).
        Xsmooth: result of applying smoothing to X
        yvals: numpy array of dimensions (nt,ny,nx).
        xvals: numpy array of same dimensions as xvals.
        plottitle: string, self-explanatory.
        plotname: string, name under which the plot is to be saved.
        plotpath: string, directory where plot is to be saved.
    """

    dX = X - Xsmooth

    Xplot = [X, Xsmooth, dX]

    nt = X.shape[0]
    npl = len(Xplot) #number of plots at each time

    t = np.array(getdata('xtimes')).flatten()/60

    fig, axs = plt.subplots(npl, nt, figsize=(nt*(5+2),5*npl), sharey='row', sharex='col')

    for i in range(nt):
        print(i)
        for j in range(npl):
            print(j)
            Y = Xplot[j][i,]
            print("a")
            v_min = np.min(Y)
            v_max = np.max(Y)
            
            im = axs[j,i].pcolormesh(xvals[i,],yvals[i,], Y, cmap='seismic',vmin=v_min,vmax=v_max)
            print("b")
            axs[j,i].set_title('t'+ '%0.2d' % t[i])
            axs[j,i].set_xlim([np.min(xvals),np.max(xvals)])
            axs[j,i].set_ylim([np.min(yvals),np.max(yvals)])
            axs[j,i].set_xticks([k for k in np.linspace(np.min(xvals),np.max(xvals),5)])
            axs[j,i].set_yticks([k for k in np.linspace(np.min(yvals),np.max(yvals),10)])
            axs[j,i].xaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f'))
            axs[j,i].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2f'))
            print("c")
            cbar = fig.colorbar(im,ax=axs[j,i])
            cbar.set_ticks([i for i in np.linspace(v_min,v_max,9)])
            print("d")

    fig.suptitle(plottitle)
    print("Saving plot as " + plotpath+plotname + ".png")
    fig.savefig(plotname,bbox_inches="tight")


def makesmooth(str_id, passes):

    X = getdata(str_id)
    Xsmooth = wrf.smooth2d(X,passes)

    return Xsmooth



def detrend(X):

    nym1 = X.shape[0]-1
    nxm1 = X.shape[1]-1
    
    Xdetrend = np.zeros((X.shape))
    Xdetrend = Xdetrend[1:,1:,]

    I = np.arange(0,nxm1)
    J = np.arange(0,nym1)

    JX = np.tensordot(J,X[nym1,:-1,] - X[0,:-1,],0)
    XI = np.tensordot(I,X[:-1,nxm1,] - X[:-1,0,],0)
    XI = np.moveaxis(XI, [0,1], [1,0])
    JI = np.tensordot(J,I,0)
    JIX = np.tensordot(JI,X[nym1,nxm1,] - X[nym1,0,] - X[0,nxm1,] + X[0,0,],0)

    Xdetrend[:,:,] = X[:-1,:-1,] - JX/nym1 - XI/nxm1 + JIX/(nym1*nxm1)
    
    return(Xdetrend)


def errico(Y,dx,dy):

    #See Errico 1985, page 1555. Y should be a de-trended field

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

    """Convert lat-lon coordinates to km. This assumes the Earth is spherical."""

    R = 6371
    pi = np.pi

    xnew = x*R*np.cos(y*pi/180)*pi/180
    ynew = y*R*pi/180

    return xnew, ynew


def deletefile(filename_full):
    
    """Input argument should be a string giving the full path to the file. 
       This function checks whether the corresponding file exists then if so
       deletes it.
    """
    
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
