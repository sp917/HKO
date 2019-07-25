import wrf
import numpy as np
import math
import os
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import pandas as pd

plotpath = '/mnt/c/Users/hko/Desktop/plots/'
datapath = '/home/sp917/data/'

######################################################################################
######################################################################################
######################################################################################

class DATA:
   
    def __init__(self, str_id, h=None, passes=None, pars=None, initialise='Full', \
            autoplot = True, scaling_type='tanh', load_coords=True):

        """
        str_id (str) :: String specifying name of variable to read.

        """
        
        self.str_id = str_id

        if str_id in ['ke', 'speed']:
            self.U = DATA('ua', h=h, passes=passes, pars=pars,\
                    autoplot=False, initialise=False, \
                    scaling_type=scaling_type, load_coords=False)
            self.V = DATA('va', h=h, passes=passes, pars=pars,\
                    autoplot=False, initialise=False, \
                    scaling_type=scaling_type, load_coords=False)
            U = np.array(wrf.to_np(self.U.X))
            V = np.array(wrf.to_np(self.V.X)) 
            self.derived=True
            self.X = switch[str_id](U,V)
        else:
            self.X = getdata(str_id)
            self.derived = False

        if load_coords:
            self.yvals = getdata('lat')
            self.xvals = getdata('lon')
            self.dx, self.dy = deltayx(self.yvals, self.xvals)
            self.K = getK(self.dy, self.dx)


            #Note: wspd_dir and sqrt(ua^2 +va^2) differ by order 10^-6 
                    # but we take them to be effectively the same

        self.z = getdata('z')
        self.isvertical=self.checkifvertical()

        self.Xsmooth = None
        self.passes = passes

        self.Xinterp = None

        self.h = h
        self.Xsmoothinterp = None
        self.dX = None
        
        self.S = None
        self.Ssmooth = None
        self.dS = None

        self.Xhy = None 
        self.pars = pars
        self.scaling_type = scaling_type

        self.Shy = None

        self.autoplot = autoplot
        self.A = 1
 
        if (initialise == 'Partial') or (initialise=='Full'):
            self.set_smoothness(self.passes)
            self.set_height(self.h)
            if pars==None:
                self.set_parameters()
            self.add_turbulence()
        if initialise == 'Full':
            self.calc_spectra()
            self.spectrum_turb_added()
  
  ###############################################################################
    
    def checkifvertical(self):

        if np.all(self.X==None):
            return None
        else:
            if self.X.shape[-4:]==self.z.shape:
                print("Variable " + self.str_id + " has shape " + str(self.X.shape) +\
                        ". Assuming defined on vertical layers.")
                isvertical=True
            else:
                print("Variable " + self.str_id + " has shape " + str(self.X.shape) + \
                        ". Assuming not defined on vertical layers.")
                isvertical=False
        
            return isvertical

    ###############################################################################
    
    def set_smoothness(self, passes, update=True):

        """
        passes (int) :: number of times to apply smoothing to field.

        """
        if passes==None:
            passes=100

        filename = datapath + self.str_id + "_" + str(passes)
        
        try:
            f = open(filename + '.npy')
            f.close()
        except:    
            print("Smoothing " + self.str_id + " " + str(passes) + " times...")
            self.Xsmooth = wrf.smooth2d(self.X, passes)
            print("Done.")
            np.save(filename, self.Xsmooth)
        else:
            print("Reading data from " + filename)
            Xsmooth = np.load(filename + '.npy')
            self.Xsmooth=wrf.smooth2d(self.X,1)
            self.Xsmooth[:,] = Xsmooth

        self.passes = passes

        if self.derived:
            self.U.set_smoothness(passes, update=False)
            self.V.set_smoothness(passes, update=False)
        
        if update:
            self.update(1)

        return
    
    ###############################################################################
    
    def set_height(self, h=1000, update=True):

        """
            h (float) :: height in metres.
        """

        if h==None:
            h=1000

        if np.any(self.Xsmooth==None):
            print("Field must first be smoothed. Run set_smoothness().") 
            return

        if self.isvertical==False:
            print("Variable not defined on vertical layers. Not interpolating.")
            self.Xinterp=self.X
            self.Xsmoothinterp=self.Xsmooth
        else:
            zmin = np.array(wrf.to_np(np.min(np.max(self.z,axis=(2,3)))))
            zmax = np.array(wrf.to_np(np.max(np.min(self.z,axis=(2,3)))))

            if (h<zmin) or (h>zmax):
                print("Height = " + str(h) + \
                        "\nMay produce missing values if not between " + str(zmin) + " and " + str(zmax)) 

            print("Interpolating " +self.str_id+ " onto height " + str(h) + "...")

            self.Xinterp = wrf.interplevel(self.X, self.z, h)            
            self.Xsmoothinterp = wrf.interplevel(self.Xsmooth, self.z, h)
            
            print("Done.")

            nans = np.sum(np.isnan(np.array(wrf.to_np(self.Xinterp))))
            if nans > 0:
                print("Warning: Field contains " + str(nans) + " NaNs.")
            elif nans==np.size(self.Xinterp):
                print("Warning: Field only contains NaNs.")
            else:
                print("Field contains no NaNs.")
            self.h = h


        self.dX = self.Xinterp - self.Xsmoothinterp
        
        if self.derived:
            self.U.set_height(h, update=False)
            self.V.set_height(h, update=False)

        if update:
            self.update(2)
            self.update(3)

        return


    ###############################################################################

    def calc_spectra(self):
 
        if np.any(self.Xsmoothinterp==None):
            print("Field must first be interpolated. Run set_height().")
            return
        
        print("Calculating spectrum 1 of 3...")
        self.S = spectrum(self.Xinterp, self.dy, self.dx)
        print("Calculating spectrum 2 of 3...")
        self.Ssmooth = spectrum(self.Xsmoothinterp, self.dy, self.dx)
        print("Calculating spectrum 3 of 3...")
        self.dS = spectrum(self.dX, self.dy, self.dx)
        print("Done.")

        return

    ###############################################################################
    
    def set_parameters(self, a=1, b=2000, c=0, d=-1, e=1, mu=0.5, alt=0, update=True):

        """
        a, b, c, d, e, mu (float) :: parameters for scaling the addition of 
                                     turbulence.
        """

        self.pars = params(a,b,c,d,e,mu,alt)
        
        if self.derived:
            self.U.set_parameters(a,b,c,d,e,mu,alt, update=False)
            self.V.set_parameters(a,b,c,d,e,mu,alt, update=False)
        
        if update:
            self.update(3)
        
        return
    
    ###############################################################################
    
    def change_parameters(self, a=None, b=None, c=None, d=None, e=None, mu=None,\
            alt=None, update=True):

        """
        a, b, c, d, e, mu (float) :: parameters for scaling the addition of 
                                     turbulence. This function is the same 
                                     as set_parameters but the default values 
                                     are set to those already in use.
        """
        if a==None:
            a = self.pars.a
        if b==None:
            b = self.pars.b
        if c==None:
            c = self.pars.c
        if d==None:
            d = self.pars.d
        if e==None:
            e = self.pars.e
        if mu==None:
            mu = self.pars.mu
        if alt==None:
            alt = self.pars.alt

        self.pars = params(a,b,c,d,e,mu,alt)

        if self.derived:
            self.U.set_parameters(a,b,c,d,e,mu,alt, update=False)
            self.V.set_parameters(a,b,c,d,e,mu,alt, update=False)
        
        if update:
            self.update(3)

        return

    ###############################################################################

    def calc_scaling(self,z):

        mu=self.pars.mu
        a = self.pars.a
        b = self.pars.b
        c = self.pars.c
        d = self.pars.d
        e = self.pars.e

        if self.scaling_type in ["constant", "const", "c", "Constant", "Const", "C"]:
            scaling = 2*mu
            self.scaling_type = "constant"
        elif self.scaling_type in ["tanh", "Tanh"]:
            scaling = mu*(1-np.tanh(a*(z-b)/(c*z-d))/e)
            self.scaling_type = "tanh"
        elif self.scaling_type in ["linear","Linear", "lin", "l", "Lin", "L"]:
            scaling = mu*(1 - a*z)
            self.scaling_type = "linear"
        else:
            scaling = 2*mu
            self.scaling_type = "constant"

        if self.derived:
            self.U.calc_scaling(z)
            self.V.calc_scaling(z)

        return scaling
    
    ###############################################################################
    
    def add_turbulence(self, update=True):

        if np.any(self.Xsmoothinterp==None) or (self.pars==None):
            print("Field must first be interpolated and parameters" \
                   + " must be defined. Run set_height(), set_parameters().")
            return
            
        if self.str_id=='ke':
            self.U.add_turbulence(update=False)
            self.V.add_turbulence(update=False)
            
            Xhy = 0.5*(self.U.Xhy*self.U.Xhy + self.V.Xhy*self.V.Xhy)
        
        else:

            alt=self.pars.alt

            print("Adding turbulence field to intial field...")

            X0 = self.Xinterp[0,]
            dXfin = self.dX[-1,]
            Xsmoothfin = self.Xsmoothinterp[-1,]

            if self.str_id.endswith('10'): 
                scaling = self.calc_scaling(10)
            else:
                scaling = self.calc_scaling(self.h)

            print(scaling)

            if alt==0:
                Xhy = X0 + scaling*dXfin
            elif alt==1:
                dXfin = dXfin - np.mean(dXfin)
                Xhy = X0*(1 + scaling*dXfin/np.max(np.abs(dXfin)))
            elif alt==2:
                Xhy = X0*(1 + scaling*dXfin/np.max(np.abs(dXfin)))
            elif alt==3:
                dXfin = dXfin - np.mean(dXfin)
                Xhy = X0 + scaling*dXfin
            elif alt==4:
                Xhy = X0*(1 + scaling*dXfin/np.max(np.abs(Xsmoothfin)))
            else:
                alt = 0
                Xhy = X0 + scaling*dXfin
        
            print("Done.")

        self.Xhy = Xhy
        
        if update:
            self.update(4)

        if self.autoplot:
            plot_field_turb(self)

        return

    ###############################################################################
    
    def spectrum_turb_added(self):


        if np.any(self.Xhy)==None:
            print("Hybrid field must first be calculated. " \
                    +" Run add_turbulence().")
            return

        print("Calculating spectrum for turbulence-added field...")

        self.Shy = spect1d(self.Xhy, self.dy[-1,], self.dx[-1])

        print("Done.")

        print("h = ", self.h)
        print(np.max(self.Shy-self.S[-1]))

        if self.autoplot:
            plot_spectra_turb(self)

        return

    ###############################################################################

    def update(self,where=0):

        """
        where (int) :: specify when to begin updating.

        """

        if where==0:
            if not (self.passes==None or self.Xsmooth==None):
                self.set_smoothness(self.passes)
        if where==1:
            if not (self.h==None or np.any(self.Xsmooth==None) or self.Xinterp==None):
                self.set_height(self.h)
        if where==2:
            if not (np.all(self.Xsmoothinterp==None) or self.S==None):
                self.calc_spectra()
        if where==3:
            if not (np.all(self.Xsmoothinterp==None) or self.pars==None \
                    or np.all(self.Xhy==None)):
                self.add_turbulence()
        if where==4:
            if not (np.all(self.Xhy==None) or np.all(self.Shy==None)):
                self.spectrum_turb_added()

        return

    ###############################################################################

    def change_options(self, autoplot=None, scaling_type=None):

        if not autoplot==None:
            self.autoplot = autoplot
        if not scaling_type==None:
            self.scaling_type = scaling_type
            self.update(3)

        return

    ###############################################################################

    def print_info(self):

        print("str_id = ", self.str_id)
        print("isvertical = ", self.isvertical)
        print("h = ", self.h)
        print("X.shape = ", self.X.shape)

        if np.any(self.Xinterp==None):
            print("Xinterp = ", self.Xinterp)
            print("Xsmoothinterp = ", self.Xsmoothinterp)
        else:
            print("Xinterp.shape = ", self.Xinterp.shape)
            print("Xsmoothinterp.shape = ", self.Xsmoothinterp.shape)
        
        if np.any(self.Xsmooth==None):
            print("Xsmooth = ", self.Xsmooth)
            print("dX = ", self.dX)
        else:
            print("Xsmooth.shape = ", self.Xsmooth.shape)
            print("dX.shape = ", self.dX.shape)

        print("passes = ", self.passes)

        if self.S==None:
            print("S = ", self.S)
            print("Ssmooth = ", self.Ssmooth)
            print("dS = ", self.dS)
        else:
            print("len(S) = ", len(self.S))
            print("len(Ssmooth)= ", len(self.Ssmooth))
            print("len(dS) = ", len(self.dS))

        if np.any(self.Xhy==None):
            print("Xhy = ", self.Xhy)
        else:
            print("Xhy.shape = ", self.Xhy.shape)

        print("pars = ", self.pars)
        
        if not self.pars==None:
            print("a = ", self.pars.a)
            print("b = ", self.pars.b)
            print("c = ", self.pars.c)
            print("d = ", self.pars.d)
            print("e = ", self.pars.e)
            print("mu = ", self.pars.mu)
            print("alt = ", self.pars.alt)

        if np.any(self.Shy)==None:
            print("Shy = ", self.Shy)
        else:
            print("Shy.shape = ", self.Shy.shape)


    ###############################################################################


######################################################################################
######################################################################################

def plot_field_turb(data):

    fig, axs = plt.subplots(1,5, figsize = (35, 5))

    Xplots = [ data.Xinterp[0,], data.Xinterp[-1,], data.Xsmoothinterp[-1,], \
            data.dX[-1,], data.Xhy]

    titles = ["initial", "final", "smooth", "turbulence", "hybrid"]

    xvals = data.xvals[0,]
    yvals = data.yvals[0,]

    for i in range(5):

        Y = Xplots[i]

        Y = np.array(wrf.to_np(Y))

        wherenans = np.isnan(Y)
        notnans = (1-wherenans)==1

        v_min = np.min(Y[notnans])
        v_max = np.max(Y[notnans])

        Y[wherenans] = -10e20

        im = axs[i].pcolormesh(xvals,yvals, Y, cmap='seismic', vmin=v_min, vmax=v_max)
        axs[i].set_xlim([np.min(xvals),np.max(xvals)])
        axs[i].set_ylim([np.min(yvals),np.max(yvals)])
        axs[i].set_xticks([k for k in np.linspace(np.min(xvals),np.max(xvals),5)])
        axs[i].set_yticks([k for k in np.linspace(np.min(yvals),np.max(yvals),10)])
        axs[i].xaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f'))
        axs[i].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2f'))
        
        cbar = fig.colorbar(im,ax=axs[i])
        cbar.set_ticks([i for i in np.linspace(v_min,v_max,9)])

        axs[i].set_title(titles[i])
    

    if data.isvertical:
        plottitle = data.str_id +  " at height " + str(data.h) + " with " \
                + str(data.passes) + " smoothings."
        plotname = data.str_id + "_height" + str(data.h) + "_smooth" \
                + str(data.passes) + "_turb_added_field_method" + str(data.pars.alt) 
    else:
        plottitle = data.str_id +  " with " + str(data.passes) + " smoothings."
        plotname = data.str_id + "_smooth" + str(data.passes)  \
                + "_turb_added_field_method" + str(data.pars.alt)
    
    if data.scaling_type=='tanh':
        plotname=plotname + "_params" + str(data.pars.a) + "_" \
                + str(data.pars.b) + "_" + str(data.pars.c) + "_" \
                + str(data.pars.d) + "_" + str(data.pars.e) \
                + "_" + str(data.pars.mu) 
    elif data.scaling_type=='constant':
       plotname = plotname + "_mu" + str(data.pars.mu)  
    elif data.scaling_type=='linear': 
       plotname = plotname + "_mu" + str(data.pars.mu)+ "_a" + str(data.pars.a)
    else:
       plotname = plotname + "_unknown_scaling_type" 

    print("Saving plot as " + plotpath + plotname)
    
    fig.suptitle( plottitle )
    fig.savefig(plotpath+plotname+'.png', bbox_inches="tight",format='png')
    plt.close()

    return

######################################################################################

def plot_spectra_turb(data):
    
    Xplots = [ data.S[-1], data.S[0], data.dS[-1], data.Shy ] 
    labels = ["Final field", "Initial field", "Turbulence", "Hybrid field"]
    

    fig = plt.figure()
    ax1 = plt.subplot(111)

    lims = []
    K = data.K[0]
    for i in range(len(Xplots)):

        Y = Xplots[i]
        Y = np.array(wrf.to_np(Y))

        v_min = np.min(Y)
        v_max = np.max(Y)
    
        im = ax1.plot(K,Y, label=labels[i])

        lims = lims + [v_min,v_max]

    Kline = np.arange(0.1*np.min(K), 10*np.max(K))
    line = data.A*np.power(Kline, -5/3)

    ax1.plot(Kline, line, 'k--', label=r'$k^{-5/3}$')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_ylim(min(lims),max(lims))
    ax1.set_xlim(ax1.get_xlim())
    ax1.legend(loc='best')
    ax1.set_xlabel("Wavenumber")
    ax1.set_ylabel("S(k)")
    
    ax2 = ax1.twiny()
    ax2.plot(2*np.pi/K,K,alpha=0.0)
    ax2.set_xscale('log')
    ax2.set_xlim(ax2.get_xlim()[::-1])
    ax2.set_xlabel("Wavelength (km)")

    
    if data.isvertical:
        plottitle = "Spectra of " + data.str_id +  " at height " + str(data.h) \
                + " with " + str(data.passes) + " smoothings."
        plotname = data.str_id + "_spectra_height" + str(data.h) + "_smooth" \
                + str(data.passes) + "_turb_added_field_method" + str(data.pars.alt) 
    else:
        plottitle = "Spectra of " + data.str_id + " with " + str(data.passes) \
                + " smoothings."
        plotname = data.str_id + "_spectra_smooth" + str(data.passes) \
                + "_turb_added_field_method" + str(data.pars.alt) 
    
    if data.scaling_type=='tanh':
        plotname=plotname + "_params" + str(data.pars.a) + "_" \
                + str(data.pars.b) + "_" + str(data.pars.c) + "_" \
                + str(data.pars.d) + "_" + str(data.pars.e) \
                + "_" + str(data.pars.mu) 
    elif data.scaling_type=='constant':
       plotname = plotname + "_mu" + str(data.pars.mu)  
    elif data.scaling_type=='linear': 
       plotname = plotname + "_mu" + str(data.pars.mu)+ "_a" + str(data.pars.a)
    else:
       plotname = plotname + "_unknown_scaling_type" 
    
    plt.title(plottitle, y = 1.2)

    print("Saving plot as " + plotpath + plotname)
    fig.savefig(plotpath+plotname+'.png', bbox_inches="tight",format='png')
    plt.close()

    return

######################################################################################
######################################################################################

class params:
    def __init__(self,a=0.001,b=1000,c=0,d=-1,e=1,mu=0.5,alt=0):
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.e = e
        self.mu = mu
        self.alt = alt
      
    def __str__(self):
        return "(%g, %g, %g, %g, %g, %g, %d)" %  (self.a,self.b,self.c,self.d,self.e,self.mu,self.alt)

    def __eq__(self,other):
        if not (type(other)==type(self)):
            return False
        else:
            alleq = (self.a==other.a) and (self.b==other.b) and (self.c==other.c) and \
                    (self.d==other.d) and (self.e==other.e) and (self.mu==other.mu) and\
                    (self.alt==other.alt)
            return alleq

######################################################################################
######################################################################################


def getdata(str_id):

    files = [Dataset(datapath+f) for f in os.listdir(datapath) \
        if f.startswith('wrfout_xhka')]    

    X = wrf.getvar(files, str_id, timeidx=wrf.ALL_TIMES, method='join', meta=True)
   
    if str_id.startswith('wspd_wdir'):
        X = X[0,] #in this case we choose wind speed

    return X

######################################################################################

def deltayx(yvals,xvals):

    """
        yvals and xvals should be numpy arrays of shape nt*ny*nx
    
    """

    yvals, xvals = deg2km(yvals,xvals)

    dy = np.zeros((yvals.shape))
    dx = np.zeros((xvals.shape))

    dy=dy[:,1:,:] #the y dimension should be one less than that of yvals
    dx=dx[:,:,1:]

    for n in range(dy.shape[1]):
        dy[:,n,:] = yvals[:,n+1,:] - yvals[:,n,:]
    for m in range(dx.shape[2]):
        dx[:,:,m] = xvals[:,:,m+1] - xvals[:,:,m]

    return dy,dx

######################################################################################

def deg2km(y,x):

    """Convert lat-lon coordinates to km. This assumes the Earth is spherical."""

    R = 6371
    pi = np.pi

    xnew = x*R*np.cos(y*pi/180)*pi/180
    ynew = y*R*pi/180

    return  ynew, xnew

######################################################################################

def getK(dy,dx):

    """
        dy and dx should have dimensions nt*(ny-1)*nx and nt*ny*(nx-1) respectively

    """
    K = []
    for t in range(dy.shape[0]):
        K1 = getAkmaxK(dy[t,],dx[t,])[4]
        K = K + [K1]

    return K

#######################################################################################

def getAkmaxK(dyt,dxt):
    
    nym1 = dyt.shape[0]
    nxm1 = dxt.shape[1]

    DY = np.mean(dyt, axis=(0,1))
    DX = np.mean(dxt, axis=(0,1))

    aq = DY*nym1
    ap = DX*nxm1

    A = max(aq,ap)

    if abs(A-ap)<0.00000000001:
        kmax = (nym1+1)//2
    else:
        kmax = (nxm1+1)//2

    K = np.array(range(1,kmax))
    K = 2*np.pi*K/A

    return A, ap, aq, kmax, K

######################################################################################

def spectrum(X, dy, dx):

    """
    Inputs:

    X (numpy array), dimensions nt*ny*nx
    dy (numpy array), dimensions nt*(ny-1)*nx
    dx (numpy array), dimensions nt*ny*(nx-1)

    Outputs:

    S (list), contains nt numpy arrays of possibly different dimensions
    
    """

    S = []

    X = detrend(X)

    Y = np.fft.fftn(X,axes=(1,2))

    for t in range(Y.shape[0]):
        S1 = errico(Y[t,],dy[t,],dx[t,])
        S = S + [S1]

    return S

######################################################################################

def fill_nans(X):

    """X should be a numpy array of dims (ny,nx)"""
    
    Xnew = pd.DataFrame(X, index = np.arange(X.shape[0]))
#    Xnew = Xnew.fillna(method="ffill", axis=0)
#    Xnew = Xnew.fillna(method="bfill", axis=0)
#    Xnew = Xnew.fillna(method="ffill", axis=1)
#    Xnew = Xnew.fillna(method="bfill", axis=1)

    Xnew = Xnew.interpolate(limit_direction='both')

    Xnew = np.array(Xnew)

    return Xnew

def detrend(X):

    """ X should have shape (nt,ny,nx) or (ny,nx) """
    
    X = np.array(wrf.to_np(X))

    X = np.moveaxis(X, [-2,-1], [0,1])

    #wherenans = np.isnan(X)
    #X[wherenans] = 0
    
    numnans = np.sum(np.isnan(X))

    if numnans>0:
        print(str(numnans) + " NaNs in field. Interpolating so that spectra" \
                + " may be calculated.")
        if len(X.shape)==3:
            for t in range(X.shape[0]):
                X[t,] = fill_nans(X[t,])
        else:
            X = fill_nans(X) 
        numnans = np.sum(np.isnan(X))
        if numnans==0:
            print("Interpolation successful: " + str(numnans) + " NaNs in field.")
        else:
            print("Warning: field still contains " + str(numnans) + " NaNs.")

    ny = X.shape[0]
    nx = X.shape[1]

    I = np.arange(0,nx)
    J = np.arange(0,ny)

    JX = np.tensordot(J,X[-1,] - X[0,],0)
    XI = np.tensordot(I,X[:,-1,] - X[:,0,],0)
    XI = np.moveaxis(XI, [0,1], [1,0])
    JI = np.tensordot(J,I,0)
    JIX = np.tensordot(JI,X[-1,-1,] - X[-1,0,] - X[0,-1,] + X[0,0,],0)

    Xdetrend = X - JX/(ny-1) - XI/(nx-1) + JIX/((ny-1)*(nx-1))
    
    Xdetrend = Xdetrend[:-1,:-1,]

    Xdetrend = np.moveaxis(Xdetrend, [0, 1], [-2,-1] )

    return(Xdetrend)

######################################################################################

def errico(Y,dyt,dxt):

    #See Errico 1985, page 1555. Y should be the Fourier transform of a de-trended two-dimensional field. dy and dx should have dimensions (ny-1)*nx and ny*(nx-1) respectively.

    nym1 = Y.shape[0]
    nxm1 = Y.shape[1]

    A, ap, aq, kmax = getAkmaxK(dyt,dxt)[0:4]

    S = np.zeros(kmax)
    
    Yrescaled = Y/(nxm1*nym1)
    
    for k in range(0,kmax):
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

    return S[1:]

######################################################################################

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

#####################################################################################

def KE(U,V):
    return 0.5*(U*U + V*V)

def speed(U,V):
    return np.sqrt(U*U +V*V)

switch = { 'ke' : KE,
           'speed' : speed
        }

#####################################################################################
#####################################################################################
