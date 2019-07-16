import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import sys
import os
import routines 
def main():

    var = sys.argv[1]
    year = '%0.4d' % int(sys.argv[2])
    month = '%0.2d' % int(sys.argv[3])
    day = '%0.2d' % int(sys.argv[4])
    hour = '%0.2d' % int(sys.argv[5])
    layer = int(sys.argv[6])
    ny = int(sys.argv[7])
    nx = int(sys.argv[8])
    filt = int(sys.argv[9])

    path = "/home/sp917/"
    datapath = path + "data/"
    plotpath = path + "plots/"
    date =  year + '-' + day + '-' + month + '_' + hour + ':00:00'
    prefix="Filtered" + str(filt) + "_" + str(ny) + "_" + str(nx) + "_" 
    t = routines.gettimes(datapath,date)

    filename = prefix + var + '_' + date + '.nc'
    data = Dataset(datapath+filename,"r",format="NETCDF4_CLASSIC")
    
    X = np.array(data.variables[var])
    Y = np.array(data.variables[var+"_filtered"])
    xvals = np.array(data.variables["XLONG"])
    yvals = np.array(data.variables["XLAT"])

    if len(X.shape)==4:
        X = X[:,layer,:,:]
        Y = Y[:,layer,:,:]
        layers = True
    else:
        layers = False
    
    nt = X.shape[0]
    fig, axs = plt.subplots(3, nt, figsize=(nt*(5+2),5*3), sharey='row', sharex='col')
    Xall = [X, Y, X-Y] 
    print(t)   
    for i in range(nt):
        for j in range(len(Xall)):
            Xplot = Xall[j][i,]
            v_min = np.min(Xplot)
            v_max = np.max(Xplot)
            print(Xplot.shape)
            im = axs[j,i].pcolormesh(xvals[i,],yvals[i,],Xplot, cmap='seismic',vmin=v_min,vmax=v_max)
            axs[j,i].set_title('t'+t[i])
            axs[j,i].set_xlim([np.min(xvals),np.max(xvals)])
            axs[j,i].set_ylim([np.min(yvals),np.max(yvals)])
            axs[j,i].set_xticks([k for k in np.linspace(np.min(xvals),np.max(xvals),5)])
            axs[j,i].set_yticks([k for k in np.linspace(np.min(yvals),np.max(yvals),10)])
            axs[j,i].xaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f'))
            axs[j,i].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2f'))
            cbar = fig.colorbar(im,ax=axs[j,i])
            cbar.set_ticks([i for i in np.linspace(v_min,v_max,9)])
    
    if layers:
        plotname = prefix + var + '_layer_' + str(layer) +'_' + date
        title=var + " in layer " + str(layer)
    else: 
        plotname = prefix + var + '_' + date 
        title=var

    fig.suptitle(title)
    print("Saving plot as " + plotname + ".png")
    fig.savefig(plotpath+plotname,bbox_inches="tight")



main()




"""
        im1 = axs[0,i].pcolormesh(xvals[i,],yvals[i,],X[i,],cmap='seismic')
        im2 = axs[1,i].pcolormesh(xvals[i,],yvals[i,],Y[i,],cmap='seismic')
        im3 = axs[2,i].pcolormesh(xvals[i,],yvals[i,],X[i,]-Y[i,],cmap='seismic')
        cbar1 = fig.colorbar(im1,ax=axs[0,i])
        cbar2 = fig.colorbar(im2,ax=axs[1,i])
        cbar3 = fig.colorbar(im3,ax=axs[2,i])
"""
