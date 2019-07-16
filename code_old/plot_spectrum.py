import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import sys

def main():
    
    var = sys.argv[1]
    year = '%0.4d' % int(sys.argv[2])
    month = '%0.2d' % int(sys.argv[3])
    day = '%0.2d' % int(sys.argv[4])
    hour = '%0.2d' % int(sys.argv[5])
    t = '%0.2d' % int(sys.argv[6])
    layer = sys.argv[7]
    a = float(sys.argv[8])

    path = "/home/sp917/"
    datapath = path + "data/"
    plotpath = path + "plots/"

    date = year + '-' + day + '-' + month + '_' + hour + ':00:00.t' + t

    if layer=='average2': 
        datafile = datapath + "Spectrum_vertav_" + var + "_" + date + ".nc"
        f = Dataset(datafile, "r", format="NETCDF4_CLASSIC")
        K = np.array(f.variables['k'])[1:]
        S = np.array(f.variables['spectrum'])[0,1:]
        print(S.shape)
    else:
        layer=int(layer)
        datafile = datapath + "Spectrum_" + var + "_" + date + ".nc"
        f = Dataset(datafile, "r", format="NETCDF4_CLASSIC")
        K = np.array(f.variables['k'])[1:]
        S = np.array(f.variables['spectrum'])[0,layer,1:]

    #a = np.max(S)*np.power(K[0],5/3)
    line = a*np.power(K,-5/3)


    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    ax2 = ax1.twiny()

    p1 = ax1.plot(K,S,'k', label='Spectrum')
    p2 = ax1.plot(K,line, 'k--', label=r'$k^{-5/3}$')

    ax1.set_xscale('log')
    ax1.set_yscale('log')

    ax1.set_ylim([np.min(S),np.max(S)])

    ax1.set_xlabel(r'$k$ $(km^{-1})$')
    ax1.set_ylabel('$S(k)$')

    ax1.legend(loc='best')
    
    ax2.plot(2*np.pi/K, K, alpha=0.0)
    ax2.set_xscale('log')
    ax2.set_xlim(ax2.get_xlim()[::-1])
    ax2.set_xlabel("Wavelength (km)")
    
    fig = plt.gcf()
    plotname='spectrum_' + var + '_layer_' + str(layer) + '_' + date + '.png'
    fig.savefig(plotpath+plotname, bbox_inches="tight")

main()
