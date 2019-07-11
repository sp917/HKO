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
    a1 = float(sys.argv[7])
    a2 = float(sys.argv[8])

    path = "/home/sp917/"
    datapath = path + "data/"
    plotpath = path + "plots/"

    date =  year + '-' + day + '-' + month + '_' + hour + ':00:00'

    if layer=='average2':
        prefix = 'Spectrum_vertav_' + var + '_' + year + '-' + day + '-' + month + '_' + hour + ':00:00.t'
    else:
        prefix = 'Spectrum_' + var + '_' + year + '-' + day + '-' + month + '_' + hour + ':00:00.t'

    files = [ f for f in os.listdir(datapath) if (f.startswith(prefix) and f.endswith('.nc') ) ]

    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    ax2 = ax1.twiny()
   
    lims = []
    for f in files:
        data = Dataset(datapath + f, "r", format="NETCDF4_CLASSIC")

        t = f.replace(prefix,'')
        t = t.replace('.nc','')
        print(t)
        
        if layer=='average1':
            S=np.array(data.variables['spectrum'])[0,:,1:]
            S=np.mean(S,axis=0)
        elif layer=='average2':
            S=np.array(data.variables['spectrum'][0,1:])
        else:
            layer = int(layer)
            S = np.array(data.variables['spectrum'])[0,layer,1:]

        lims = lims + [np.min(S),np.max(S)]
        K = np.array(data.variables['k'])[1:]
        ax1.plot(K,S, label='Spectrum for t' + t)
   
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    
    lims = np.array(lims)

    ax1.set_ylim(np.min(lims),np.max(lims))

    ax1.set_xlim(ax1.get_xlim())

    Kline = np.arange(0.1*np.min(K),10*np.max(K))
    line1 = a1*np.power(Kline,-5/3)
    line2 = a2*np.power(Kline,-5/3)
    ax1.plot(Kline,line1, 'k--', label=r'$k^{-5/3}$')
    ax1.plot(Kline,line2, 'k--')

    ax1.set_xlabel(r'$k$ $(km^{-1})$')
    ax1.set_ylabel('$S(k)$')

    ax1.legend(loc='best')
    
    ax2.plot(2*np.pi/K, K, alpha=0.0)
    ax2.set_xscale('log')
    ax2.set_xlim(ax2.get_xlim()[::-1])
    ax2.set_xlabel("Wavelength (km)")
   
    if layer=='average1':
        Title = "Vertically-averaged spectrum of " + var
    elif layer=='average2':
        Title = "Spectrum of vertically-averaged " + var
    else:
        Title = "Spectrum of " + var + " in layer " + str(layer)

    Title = Title + "\n"

    plt.title(Title)

    fig = plt.gcf()
    plotname='spectra_' + var + '_layer_' + str(layer) + '_' + date + '.png'
    fig.savefig(plotpath+plotname, bbox_inches="tight")

main()
