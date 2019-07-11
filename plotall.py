import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import sys
import os

def main():
    
    var = sys.argv[1]
    year = '%0.4d' % int(sys.argv[2])
    month = '%0.2d' % int(sys.argv[3])
    day = '%0.2d' % int(sys.argv[4])
    hour = '%0.2d' % int(sys.argv[5])
    layer = sys.argv[6]
    a11 = float(sys.argv[7])
    a12 = float(sys.argv[8])
    a21 = float(sys.argv[9])
    a22 = float(sys.argv[10])
    alt = 1

    path = "/home/sp917/"
    datapath = path + "data/"
    plotpath = path + "plots/"

    date =  year + '-' + day + '-' + month + '_' + hour + ':00:00'

    prefix = 'Spectrum_' + var + '_' + date + '.t'

    files = [ f for f in os.listdir(datapath) \
            if (f.startswith(prefix) and f.endswith('.nc') ) ]

    fig = plt.figure(1)
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twiny()
   
    lims = []

    for f in files:
        data = Dataset(datapath + f, "r", format="NETCDF4_CLASSIC")

        t = f.replace(prefix,'')
        t = t.replace('.nc','')
        print(t)
        
        layer = int(layer)
        S = np.array(data.variables['spectrum'])
        
        if len(S.shape)==2:
            layers=True
            S = S[layer,1:]
        elif len(S.shape)==1:
            layers=False
            S = S[1:]
        else:
            raise RuntimeError("Something went wrong!")
            
        lims = lims + [np.min(S),np.max(S)]
        K = np.array(data.variables['k'])[1:]
        ax1.plot(K,S, label='Spectrum for t' + t)
   
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    
    lims = np.array(lims)

    ax1.set_ylim(np.min(lims),np.max(lims))

    ax1.set_xlim(ax1.get_xlim())

    Kline = np.arange(0.1*np.min(K),10*np.max(K))
    line1 = a11*np.power(Kline,-5/3)
    line2 = a12*np.power(Kline,-5/3)
    ax1.plot(Kline,line1, 'k--', label=r'$k^{-5/3}$')
    ax1.plot(Kline,line2, 'k--')

    ax1.set_xlabel(r'$k$ $(km^{-1})$')
    ax1.set_ylabel('$S(k)$')

    ax1.legend(loc='best')
    
    ax2.plot(2*np.pi/K, K, alpha=0.0)
    ax2.set_xscale('log')
    ax2.set_xlim(ax2.get_xlim()[::-1])
    ax2.set_xlabel("Wavelength (km)")
  
    if layers:
        Title = "Spectrum of " + var + " in layer " + str(layer) + "\n"
    else: 
        Title = "Spectrum of " + var + "\n"
    
    plt.title(Title)

    fig = plt.gcf()
    plotname='spectra_' + var + '_layer_' + str(layer) + '_' + date + '.png'
    fig.savefig(plotpath+plotname, bbox_inches="tight")



    if alt==1:

        files = [ f for f in os.listdir(datapath) \
                if (f.startswith("UV_interp_"+date) and f.endswith('.nc') ) ]

        fignew = plt.figure(2)
        ax1 = fignew.add_subplot(111)
        ax2 = ax1.twiny()
   
        lims = []

        for f in files:
            data = Dataset(datapath + f, "r", format="NETCDF4_CLASSIC")

            t = f.replace("UV_interp_"+date+".t",'')
            t = t.replace('.nc','')
            print(t)
            
            layer = int(layer)
            S = np.array(data.variables[var+'_spectrum'])
            if len(S.shape)==3:
                S = S[0,layer,]
            elif len(S.shape)==2:
                S = S[0,]
            else:
                raise RuntimeError("Something went wrong!")
                
            lims = lims + [np.min(S),np.max(S)]
            K = np.array(data.variables['wave_number'])
            ax1.plot(K,S, label='Spectrum for t' + t)
       
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        
        lims = np.array(lims)

        ax1.set_ylim(np.min(lims),np.max(lims))

        ax1.set_xlim(ax1.get_xlim())

        Kline = np.arange(0.1*np.min(K),10*np.max(K))
        line1 = a21*np.power(Kline,-5/3)
        line2 = a22*np.power(Kline,-5/3)
        ax1.plot(Kline,line1, 'k--', label=r'$k^{-5/3}$')
        ax1.plot(Kline,line2, 'k--')

        ax1.set_xlabel(r'$k$ $(km^{-1})$')
        ax1.set_ylabel('$S(k)$')

        ax1.legend(loc='best')
        
        ax2.plot(2*np.pi/K, K, alpha=0.0)
        ax2.set_xscale('log')
        ax2.set_xlim(ax2.get_xlim()[::-1])
        ax2.set_xlabel("Wavelength (km)")
       
        if layers:
            Title = "Spectrum of " + var + " in layer " + str(layer) + "\n"
        else: 
            Title = "Spectrum of " + var + "\n"
        
        plt.title(Title)

        fignew = plt.gcf()
        plotname='alt_spectra_' + var + '_layer_' + str(layer) + '_' + date + '.png'
        fignew.savefig(plotpath+plotname, bbox_inches="tight")


main()
