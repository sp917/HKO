import numpy as np
from netCDF4 import Dataset
import math
import matplotlib.pyplot as plt
import sys
import os
import routines

def main():


    var = sys.argv[1]
    year = '%0.4d' % int(sys.argv[2])
    month = '%0.2d' % int(sys.argv[3])
    day = '%0.2d' % int(sys.argv[4])
    hour = '%0.2d' % int(sys.argv[5])
    t = '%0.2d' % int(sys.argv[6])

    date = year + '-' + day + '-' + month + '_' + hour + ':00:00.t' + t
    filename = 'UV_interp_' + date + '.nc'
    path = '/home/sp917/data/'
    
    print("Reading data for variable " + var + " from file " + filename)
    f = Dataset(path+filename,'r', format='NETCDF4_CLASSIC')

    vardata = f.variables[var]
    X = np.array(vardata)

    print(X.shape) # the shape of x should be time*height*lat*lon* or time*lat*lon

    X = X[0,]
    if len(X.shape) ==3:
        X = np.transpose(X,(1,2,0))
    elif len(X.shape)==2:
        print("")
    else:
        raise RuntimeError("Shape of X not understood")

    #the size of the time dimension of x and y is 1
    x = np.array(f.variables['XLONG'])[0,:,:]
    y = np.array(f.variables['XLAT'])[0,:,:]
    
    x,y = routines.deg2km(x,y)

    dy = routines.delta(y)
   
    x = np.transpose(x)
    dx = routines.delta(x)
    dx = np.transpose(dx)
    x = np.transpose(x)

    print("Removing linear trend from data...")

    X = routines.detrend(X)

    print("Done.")

    print("Performing fast Fourier transform on de-trended data...")

    Y = np.fft.fftn(X, axes=(0,1))
    
    print("Done.")

    print("Calculating spectrum...")

    S,K = routines.errico(Y,dx,dy)

    if len(S.shape)==2:
        S = np.transpose(S)
    elif len(S.shape)==1:
        print("")
    else:
        raise RuntimeError("Something went wrong!")
    
    print("Done.")

    print("Saving spectrum to netCDF file...")
  
    newfile_name = "Spectrum_"+ var + "_" +  date + ".nc"
   
    routines.deletefile(path+newfile_name)

    print("Creating file " + path+newfile_name)
    newfile = Dataset(path+newfile_name,'w',format='NETCDF4_CLASSIC')

    print(S.shape)

    if len(S.shape)==2:
        newfile.createDimension('wavenumber', K.shape[0])
        newfile.createDimension('bottom_top', S.shape[0])
        spectrum = newfile.createVariable('spectrum', np.float64, ('bottom_top', 'wavenumber'))
        spectrum[:,:] = S
    elif len(S.shape)==1:
        newfile.createDimension('wavenumber', K.shape[0])
        spectrum = newfile.createVariable('spectrum', np.float64, ('wavenumber'))
        spectrum[:] = S

    wavenumber = newfile.createVariable('k', np.float64, 'wavenumber')
    wavenumber[:] = K

    print("Done.")
    
    return

main()

  
