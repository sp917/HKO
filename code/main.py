import sys
import routines
import wrf
import numpy as np


def main():
    
    var = sys.argv[1]      #str, specifying variable
    passes = int(sys.argv[2]) #number of times to apply smoothing filter
    try:
        height = float(sys.argv[3])    # height in km
    except:
        height = False
        interp = False
    else:
        interp=True

    X = routines.getdata(var, interp=interp, height=height)

    yvals = routines.getdata('lat')
    xvals = routines.getdata('lon')
    #yvals = routines.getdata('XLAT')
    #xvals = routines.getdata('XLONG')

    Xsmooth = wrf.smooth2d(X,passes)

    X = wrf.to_np(X)
    Xsmooth = wrf.to_np(Xsmooth)
    yvals = wrf.to_np(yvals)
    xvals = wrf.to_np(xvals)
    
    routines.nanwarning(X)
    routines.nanwarning(Xsmooth)

    plottitle = str(passes) + "times smoothed " + var 
    plotname = "Smooth" + str(passes) + "_" + var
    
    if not (height==False):
        plottitle = plottitle + " \n at height " + str(height) + "km"
        plotname = plotname + "_height_" + str(height)

    routines.plotsmooth(X,Xsmooth,yvals,xvals,plottitle,plotname)


    dX = X - Xsmooth

    dy,dx = routines.deltayx(yvals, xvals)
   
    print(dy.shape,dx.shape)

    S = routines.spectrum(X,dy,dx)
    Ssmooth = routines.spectrum(Xsmooth,dy,dx)
    dS = routines.spectrum(dX,dy,dx)

    K = routines.getK(dy,dx)
    
    plottitle = "Spectrum for " + str(passes) + "-times smoothed " + var 
    plotname = "Spectrum_smooth" + str(passes) + "_" + var
    
    if not (height==False):
        plottitle = plottitle + " \n at height " + str(height) + "km"
        plotname = plotname + "_height_" + str(height)

    Splot = [S,Ssmooth,dS]
    routines.plotspectra( Splot, K, plottitle, plotname)

    return 0

main()
