import sys
import os
import wrf
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import math

plotpath = '/mnt/c/Users/hko/Desktop/plots/'
datapath = '/home/sp917/data/'


def getdata(str_id, interp=False, height=1):

    """ Reads in data and interoplates onto fixed horizontal levels """

    files = [Dataset(datapath+f) for f in os.listdir(datapath) \
            if f.startswith('wrfout_xhka')]    

    X = wrf.getvar(files, str_id, timeidx=wrf.ALL_TIMES, method='join', meta=True)
   
    if str_id.startswith('wspd_wdir'):
        X = X[0,]

    if interp:
        z = wrf.getvar(files, 'z', timeidx=wrf.ALL_TIMES, method='join', meta=False)
        z = z/1000 #convert to km
        nz = z.shape[1]
        zmin = np.min(np.max(z,axis=(2,3)))
        zmax = np.max(np.min(z, axis=(2,3)))
        if height < zmin or height > zmax:
            print(" Warning: Height not in appropriate range. \n Must be between " + str(zmin) + " and " + str(zmax))
        
        print('Interpolating vertically.')
        X = wrf.vinterp(files, field=X, vert_coord='ght_msl', interp_levels = [height], timeidx=wrf.ALL_TIMES)

        X = X[:,0,] #The height coordinate should have dimension 1, so we remove it
    nanwarning(X)
    print(X.shape)

    return X

def plotsmooth(X, Xsmooth, yvals, xvals, plottitle, plotname):


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

    print("xvals.shape =", xvals.shape)
    xvals = np.transpose(xvals, (0,2,1))
    yvals = np.transpose(yvals, (0,2,1))

    for i in range(nt):
        xplot = xvals[i,]
        yplot = yvals[i,]
        for j in range(npl):
            Y = Xplot[j][i,]
        
            v_min = np.min(Y)
            v_max = np.max(Y)
            
            im = axs[j,i].pcolormesh(xplot, yplot, np.transpose(Y), cmap='seismic',vmin=v_min,vmax=v_max)

            axs[j,i].set_title('t'+ '%0.2d' % t[i])
            axs[j,i].set_xlim([np.min(xvals),np.max(xvals)])
            axs[j,i].set_ylim([np.min(yvals),np.max(yvals)])
            axs[j,i].set_xticks([k for k in np.linspace(np.min(xvals),np.max(xvals),5)])
            axs[j,i].set_yticks([k for k in np.linspace(np.min(yvals),np.max(yvals),10)])
            axs[j,i].xaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f'))
            axs[j,i].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2f'))
        
            cbar = fig.colorbar(im,ax=axs[j,i])
            cbar.set_ticks([i for i in np.linspace(v_min,v_max,9)])
        

    fig.suptitle(plottitle)
    print("Saving plot as " + plotpath+plotname+'.png')
    fig.savefig(plotpath+plotname+'.png',bbox_inches="tight", format='png')
    plt.close()

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



def errico(Y,dy,dx):

    #See Errico 1985, page 1555. Y should be the Fourier transform of a de-trended two-dimensional field. dy and dx should have dimensions (ny-1)*nx and ny*(nx-1) respectively.

    nxm1 = Y.shape[1]
    nym1 = Y.shape[0]

    A, ap, aq, kmax = getAkmaxK(dy,dx)[0:4]

    S = np.zeros(kmax)
    
    Yrescaled = Y/(nxm1*nym1) #This scaling just keeps the numbers from getting too large
    
    for k in range(0,kmax):
        #if k%(kmax//10)==0: print(str(100*k/kmax) + '%')
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
    
    
#    K = np.array(range(1,kmax))
#    K = 2*np.pi*K/A

    return S[1:]


def delta(y):

    dy = np.zeros((y.shape))
    dy=dy[1:,] #the first dimension should be one less than that of y

    for n in range(0,dy.shape[0]):
        dy[n,] = y[n+1,] - y[n,]

    return dy


def deg2km(y,x):

    """Convert lat-lon coordinates to km. This assumes the Earth is spherical."""

    R = 6371
    pi = np.pi

    xnew = x*R*np.cos(y*pi/180)*pi/180
    ynew = y*R*pi/180

    return  ynew, xnew


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



def nanwarning(X):
    numnans = int(np.sum(np.isnan(X)))
    if numnans > 0:
        print( "Warning: Variable contains " + str(numnans) + " nans.")

def plotspectra(Splot, K, plottitle, plotname, A=np.ones(3)):

    """
    Splot is a list; each element of the list is also a list containing numpy arrays of various sizes.
    K is a list of numpy arrays of the same dimensions of those contained in the lists in Splot.
    """

    npl = len(Splot) #number of plots at each time

    t = np.array(getdata('xtimes')).flatten()/60

    nt = len(t)

    fig, axs = plt.subplots(1, npl, figsize=((5+2)*npl, 5+2), sharey='row', sharex='col')

    lims = []
    for j in range(npl):
        for i in range(nt):
            Y = Splot[j][i]
            Kplot = K[i]
       
            v_min = np.min(Y)
            v_max = np.max(Y)
            
            im = axs[j].plot(Kplot,Y,label=r'$t$'+str(t[i]))
            
            lims = lims + [np.min(Y)]+[np.max(Y)]

        Kline = np.arange(0.1*np.min(Kplot),10*np.max(Kplot))
        line = A[j]*np.power(Kline, -5/3)

        axs[j].plot(Kline, line, 'k--', label=r'$k^{-5/3}$')

        axs[j].set_xscale('log')
        axs[j].set_yscale('log')

        axs[j].set_ylim(min(lims),max(lims))
        axs[j].set_xlim(axs[j].get_xlim())

        axs[j].legend(loc='best')

        ax2 = axs[j].twiny()
        ax2.plot(2*np.pi/Kplot,Kplot,alpha=0.0)
        ax2.set_xscale('log')
        ax2.set_xlim(ax2.get_xlim()[::-1])


    fig.suptitle(plottitle)
    print("Saving plot as " + plotpath+plotname+'.png')
    fig.savefig(plotpath+plotname+'.png',bbox_inches="tight", format='png')
    plt.close()

def quickplots(Xplots, yvals, xvals, plottitle, plotname, titles=[]):

    """ Xplots should be a list of two-dimensional arrays, each of which has the same dimensions as yvals and xvals."""

    npl = len(Xplots)

    fig, axs = plt.subplots(1, npl, figsize=((5+2)*npl, 5))


    for j in range(npl):
        Y = Xplots[j]
    
        v_min = np.min(Y)
        v_max = np.max(Y)
        
        im = axs[j].pcolormesh(np.transpose(xvals),np.transpose(yvals), Y, cmap='seismic',vmin=v_min,vmax=v_max)

        axs[j].set_xlim([np.min(xvals),np.max(xvals)])
        axs[j].set_ylim([np.min(yvals),np.max(yvals)])
        axs[j].set_xticks([k for k in np.linspace(np.min(xvals),np.max(xvals),5)])
        axs[j].set_yticks([k for k in np.linspace(np.min(yvals),np.max(yvals),10)])
        axs[j].xaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f'))
        axs[j].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2f'))
    
        cbar = fig.colorbar(im,ax=axs[j])
        cbar.set_ticks([i for i in np.linspace(v_min,v_max,9)])

        if not len(titles)==0:
            axs[j].set_title(titles[j])
    

    fig.suptitle(plottitle + '\n \n')
    print("Saving plot as " + plotpath+plotname+'.png')
    fig.savefig(plotpath+plotname+'.png',bbox_inches="tight", format='png')
    plt.close() 


def quickplotspectra(Splot, K, plottitle, plotname, A=np.ones(5), labels=[]):

    """
    Splot is a list; each element of the list is a numpy array.
    K is a numpy array of the same dimensions .
    """

    npl = len(Splot) #number of plots at each time

    fig, axs = plt.subplots(1, npl, figsize=((5+2)*npl, 5+2), sharey='row', sharex='col')

    lims = []
    for j in range(npl):
        Y = Splot[j]
   
        v_min = np.min(Y)
        v_max = np.max(Y)
        
        im = axs[j].plot(K,Y)
        
        lims = lims + [np.min(Y)]+[np.max(Y)]

        Kline = np.arange(0.1*np.min(K),10*np.max(K))
        line = A[j]*np.power(Kline, -5/3)

        axs[j].plot(Kline, line, 'k--', label=r'$k^{-5/3}$')

        axs[j].set_xscale('log')
        axs[j].set_yscale('log')

        axs[j].set_ylim(min(lims),max(lims))
        axs[j].set_xlim(axs[j].get_xlim())

        axs[j].legend(loc='best')

        if len(labels)==npl:
            axs[j].set_title(labels[j])

        ax2 = axs[j].twiny()
        ax2.plot(2*np.pi/K,K,alpha=0.0)
        ax2.set_xscale('log')
        ax2.set_xlim(ax2.get_xlim()[::-1])


    fig.suptitle(plottitle + '\n \n')
    print("Saving plot as " + plotpath+plotname+'.png')
    fig.savefig(plotpath+plotname+'.png',bbox_inches="tight", format = 'png')
    plt.close()

def spect1d(X,dy,dx):
    """ 
    Inputs:

    X (numpy array), dimensions ny*nx
    dy (numpy array), dimensions (ny-1)*nx
    dx (numpy array), dimensions ny*(nx-1)

    Outputs:

    S (numpy array)

    """

    X = detrend(X)
    Y = np.fft.fftn(X,axes=(0,1))
    S = errico(Y,dy,dx)

    return S

def spectrum(X, dy, dx):

    """
    Inputs:

    X (numpy array), dimensions nt*ny*nx
    dy (numpy array), dimensions (ny-1)*nx*nt
    dx (numpy array), dimensions ny*(nx-1)*nt

    Outputs:

    S (list), contains nt numpy arrays of possibly different dimensions
    
    """

    S = []

    X = np.moveaxis(X, [-2,-1], [0,1])    
    X = detrend(X)
    
    nt = X.shape[2]

    Y = np.fft.fftn(X,axes=(0,1))

    for t in range(nt):
        dyt = dy[:,:,t]
        dxt = dx[:,:,t]
        S1 = errico(Y[:,:,t],dyt,dxt)
        S = S + [S1]

    return S

def deltayx(yvals,xvals):

    """
        yvals and xvals should be numpy arrays of shape nt*nx*ny
    
    """

    yvals = np.moveaxis(yvals, [-2,-1], [0,1])
    xvals = np.moveaxis(xvals, [-2,-1], [0,1])

    yvals, xvals = deg2km(yvals,xvals)

    dy = delta(yvals)

    xvals = np.moveaxis(xvals, 1,0)
    dx = delta(xvals)
    xvals = np.moveaxis(xvals, 1,0)
    dx = np.moveaxis(dx,1,0)

    return dy,dx


def getAkmaxK(dy,dx):
    
    nym1 = dy.shape[0]
    nxm1 = dx.shape[1]

    DX = np.mean(dx, axis=(0,1))
    DY = np.mean(dy, axis=(0,1))

    ap = DX*nxm1
    aq = DY*nym1

    A = max(ap,aq)

    if abs(A-ap)<0.00000000001:
        kmax = (nym1+1)//2
    else:
        kmax = (nxm1+1)//2

    K = np.array(range(1,kmax))
    K = 2*np.pi*K/A

    return A, ap, aq, kmax, K

def getK(dy,dx):

    """
        dy and dx should have dimensions (ny-1)*nx*nt and ny*(nx-1)*nt respectively

    """
    K = []
    for t in range(dy.shape[2]):
        K1 = getAkmaxK(dy[:,:,t],dx[:,:,t])[4]
        K = K + [K1]

    return K
    

def smooth_and_plot(var, passes, interp=False, height=False, doplot=True):

    X = getdata(var, interp=interp, height=height)

    yvals = getdata('lat')
    xvals = getdata('lon')

    print("Smoothing...")
    Xsmooth = wrf.smooth2d(X,passes)
    print("done")
    X = wrf.to_np(X)
    Xsmooth = wrf.to_np(Xsmooth)
    dX = X - Xsmooth

    yvals = wrf.to_np(yvals)
    xvals = wrf.to_np(xvals)

    nanwarning(X)
    nanwarning(Xsmooth)

    plottitle = str(passes) + "times smoothed " + var 
    plotname = "Smooth" + str(passes) + "_" + var
    
    if not (interp==False):
        plottitle = plottitle + " at height " + str(height) + "km"
        plotname = plotname + "_height_" + str(height)

    if doplot:
        plotsmooth(X,Xsmooth,yvals,xvals,plottitle,plotname)

    return X, Xsmooth, dX, yvals, xvals

def spectrum_and_plot(X, Xsmooth, dX, dy, dx, var, passes, height, interp, doplot=True):

    S = spectrum(X,dy,dx)
    Ssmooth = spectrum(Xsmooth,dy,dx)
    dS = spectrum(dX,dy,dx)
    K = getK(dy,dx)
 
    plottitle = "Spectrum for " + str(passes) + "-times smoothed " + var 
    plotname = "Spectrum_smooth" + str(passes) + "_" + var
    
    if not (interp==False):
        plottitle = plottitle + " at height " + str(height) + "km"
        plotname = plotname + "_height_" + str(height)

    Splot = [S,Ssmooth,dS]
    if doplot:
        print("doing plot")
        plotspectra( Splot, K, plottitle, plotname)

    return Splot, K


def add_turb_plot(X, Xsmooth, dX, yvals, xvals, var, passes, mu, alt, height, interp, doplot=True, a=8, b=0.01,c=0.07,d=0.4,e=0.003):

    isalt=str(alt)

    strmu = '%g' % mu

    if not height==False:
        scaling = 1-np.tanh(a*(height-b)/(c-d*(height-e)))
    else:       
        scaling = 1

    if alt==0:
        Xhy = X[0,] + mu*scaling*dX[2,]
    elif alt==1:
        Xhy = X[0,]*(1 + mu*dX[2,]/Xsmooth[2,])
    elif alt==2:
        Xhy = X[0,]*np.tanh(dX[2,]/Xsmooth[2,])
    elif alt==3:
        X0  = X[0,]
        R = np.abs(np.random.normal(size=(X0.shape),scale=mu))
        Xhy = (X0 + R*dX[2,])/(1 + R)
    elif alt==4:
        X0 = X[0,]
        #R = np.abs(np.random.normal(size=(X0.shape),scale=mu, loc =1))
        Xhy = X0*(1 + mu*scaling*(dX[2,]/np.max(np.abs(dX[2,]))))
    elif alt==5:
        X0 = X[0,]
        dX2 = dX[2,]
        Xsmooth2=Xsmooth[2,]
        Xhy = X0*(1 + mu*dX2/(Xsmooth2 + np.abs(dX2)))
    else:
        raise RuntimeError("Not a valid choice of method.")
    
    Xplots = [ X[0,], Xhy, dX[2,],  X[2,], Xsmooth[2,]]


    nanwarning(Xhy)

    titles = ["Field at t=0", "Hybrid field" , "Turbulence at t=12", "Field at t=12", "Smoothed field at t=12"]

    plotname = "Turbulence"+isalt+"_added_smooth" + str(passes) +  "_" + "mu" + strmu + "_" + var
    plottitle = "Turbulence adding from " + str(passes) + " smoothings.  mu =" + strmu

    
    if not (interp==False):
        plottitle = plottitle + " at height " + str(height) + "km"
        plotname = plotname + "_height_" + str(height)


    if doplot:
        quickplots(Xplots,yvals[0,],xvals[0,], plottitle,plotname, titles)
        
    return Xhy

def spect_turb_plot(Xhy, Splot, K, dy, dx, var, passes, mu, alt, height, interp, doplot=True, A=False):

    isalt=str(alt) 

    strmu = '%g' % mu
    
    Sturb = Splot[2][2]
    Sinit = Splot[0][0]
    Smooth= Splot[1][2]
    Sfin = Splot[0][2]

    Shy = spect1d(Xhy,dy[:,:,2],dx[:,:,2])

    Kplot = K[0]

    Splots = [Sinit, Shy, Sturb, Sfin, Smooth]

    plottitle = "Spectra after " + str(passes) + " smoothings. mu =" + strmu
    labels = ["Field at t=0", "Hybrid field" , "Turbulence at t=12", "Field at t=12", "Smoothed field at t=12"]
    plotname = "Turbulence"+isalt+"_added_spectra_smooth" + str(passes) + "_" + "mu" + strmu + "_" + var
    
    if not (interp==False):
        plottitle = plottitle + " at height " + str(height) + "km"
        plotname = plotname + "_height_" + str(height) 

    if A==False:
        A = np.ones(len(Splots))

    if doplot:
        quickplotspectra(Splots, Kplot, plottitle, plotname, labels=labels, A=A)

    return Splots
