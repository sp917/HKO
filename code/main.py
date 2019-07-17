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

    print(yvals.shape,xvals.shape)

    Xsmooth = wrf.smooth2d(X,passes)

    print(int(np.sum(np.isnan(Xsmooth))))

    routines.nanwarning(X)
    routines.nanwarning(Xsmooth)

    plottitle = str(passes) + "times smoothed " + var 
    plotname = "Smooth" + str(passes) + "_" + var
    
    if not (height==False):
        plottitle = plottitle + " \n at height " + str(height) + "km"
        plotname = plotname + "_height_" + str(height)

    routines.plotsmooth(X,Xsmooth,yvals,xvals,plottitle,plotname)

    X = wrf.to_np(X)
    Xsmooth = wrf.to_np(Xsmooth)
    yvals = wrf.to_np(yvals)
    xvals = wrf.to_np(xvals)
    
    X = np.moveaxis(X, [-2,-1], [0,1])
    Xsmooth = np.moveaxis(Xsmooth, [-2,-1], [0,1])
    dX = X-Xsmooth

    X = routines.detrend(X)
    Xsmooth = routines.detrend(Xsmooth)
    dX = routines.detrend(dX)

    yvals = np.moveaxis(yvals, [-2,-1], [0,1])
    xvals = np.moveaxis(xvals, [-2,-1], [0,1])

    yvals, xvals = routines.deg2km(yvals,xvals)

    dy = routines.delta(yvals)

    xvals = np.moveaxis(xvals, 1,0)
    dx = routines.delta(xvals)
    xvals = np.moveaxis(xvals, 1,0)
    dx = np.moveaxis(dx,1,0)
   
    print(dy.shape,dx.shape)

    nt = X.shape[2]

    S = []
    Ssmooth = []
    Sd = []
    K = []

    Y = np.fft.fftn(X,axes=(0,1))
    Ysmooth = np.fft.fftn(Xsmooth,axes=(0,1))
    dY = np.fft.fftn(dX,axes=(0,1))
    

    for t in range(nt):
        dyt = dy[:,:,t]
        dxt = dx[:,:,t]
        print(X.shape)
        S1,K1 = routines.errico(Y[:,:,t],dyt,dxt)
        print(S1.shape)
        S2, K2 = routines.errico(Ysmooth[:,:,t],dyt,dxt)
        S3, K3 = routines.errico(dY[:,:,t],dyt,dxt)
        S = S + [S1]
        Ssmooth = Ssmooth + [S2]
        Sd = Sd + [S3]
        K = K + [K1] #All three of K1, K2, K3 should be the same
 
    plottitle = "Spectrum for " + str(passes) + "-times smoothed " + var 
    plotname = "Spectrum_smooth" + str(passes) + "_" + var
    
    if not (height==False):
        plottitle = plottitle + " \n at height " + str(height) + "km"
        plotname = plotname + "_height_" + str(height)

    Splot = [S,Ssmooth,Sd]
    routines.plotspectra( Splot, K, plottitle, plotname)

    return 0

main()
