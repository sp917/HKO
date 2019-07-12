import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import sys
import os

def main():

    print("running python script")
    
    var = sys.argv[1]
    year = '%0.4d' % int(sys.argv[2])
    month = '%0.2d' % int(sys.argv[3])
    day = '%0.2d' % int(sys.argv[4])
    hour = '%0.2d' % int(sys.argv[5])
    layer = int(sys.argv[6])

    onecb = False

    path = "/home/sp917/"
    datapath = path + "data/"
    plotpath = path + "plots/"
    date =  year + '-' + day + '-' + month + '_' + hour + ':00:00'
    prefix="UV_interp_" 

    files = [ f for f in os.listdir(datapath) \
            if (f.startswith(prefix + date + ".t") and f.endswith(".nc")) ]
    
    print(prefix + date + ".t")
    
    nf = len(files)
    
    if onecb:
        lims = []
        for i in range(nf):
            f = files[i]
            data = Dataset(datapath+f,'r',format="NETCDF4_CLASSIC")
            X = np.array(data.variables[var])[0,]
            lims = lims + [np.max(X), np.min(X)]

        v_max = max(lims)
        v_min = min(lims)
    
    fig, axs = plt.subplots(1,nf, figsize=(nf*(5+2),5), sharey='row')
    for i in range(nf):
        f = files[i]
        data = Dataset(datapath+f,'r',format="NETCDF4_CLASSIC")
        t=f.replace(prefix+date+".t","")
        t=t.replace(".nc","")
        print(t)
        X = np.array(data.variables[var])[0,]
        xvals = np.array(data.variables['XLONG'])[0,]
        yvals = np.array(data.variables['XLAT'])[0,]

        if len(X.shape)==3:
            layers=True
            X = X[layer,]
        else:
            layers=False

        if not onecb:
            v_min = np.min(X)
            v_max = np.max(X)
        
        im = axs[i].pcolormesh(xvals,yvals,X, cmap='seismic', vmin=v_min, vmax=v_max)
        axs[i].set_title('t'+t) 
        axs[i].set_xlim([np.min(xvals),np.max(xvals)])
        axs[i].set_xticks([i for i in np.linspace(np.min(xvals),np.max(xvals),5)])
        axs[i].set_ylim([np.min(yvals),np.max(yvals)])
        axs[i].set_yticks([i for i in np.linspace(np.min(yvals),np.max(yvals),10)])
        axs[i].xaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f'))
        axs[i].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2f'))

        if not onecb:
            cbar = fig.colorbar(im, ax = axs[i])
            cont = [v_min, v_max, 0.001]
            cbar.set_ticks( [i for i in np.linspace(cont[0],cont[1],9)])

    
    fig.subplots_adjust(right=0.8)
    
    if onecb: 
        cbar_ax = fig.add_axes([0.85, 0.125, 0.01, 0.725])
        cont = [v_min, v_max, 0.001]
        cbar = fig.colorbar(im, cax=cbar_ax)
        cbar.set_ticks( [i for i in np.linspace(cont[0],cont[1],9)])
        
    if layers:
        plotname = 'Contour_' + var + '_layer_' + str(layer) +'_' + date
        title=var + " in layer " + str(layer)
    else: 
        plotname = 'Contour_' + var + '_' + date 
        title=var
        
    fig.suptitle(title)
    print("Saving plot as " + plotname + '.png') 
    fig.savefig(plotpath+plotname,bbox_inches="tight")

main()
