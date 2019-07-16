import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import sys
import os
import routines

def main():

    try:
        var = sys.argv[1]
        year = '%0.4d' % int(sys.argv[2])
        month = '%0.2d' % int(sys.argv[3])
        day = '%0.2d' % int(sys.argv[4])
        hour = '%0.2d' % int(sys.argv[5])
        ny = int(sys.argv[6])
        nx = int(sys.argv[7])
        filt = int(sys.argv[8])
    except:
        var = 'speed10'
        year = '%0.4d' %  2019
        month = '%0.2d' % 9
        day = '%0.2d' % 7
        hour = '%0.2d' % 8
        nx = 4
        ny = 4
        filt = 1

    path = "/home/sp917/"
    datapath = path + "data/"
    plotpath = path + "plots/"
    date =  year + '-' + day + '-' + month + '_' + hour + ':00:00'
    prefix="UV_interp_" 

    files = [ f for f in os.listdir(datapath) \
            if (f.startswith(prefix + date + ".t") and f.endswith(".nc")) ]
    
    for i in range(len(files)):
        f = datapath + files[i]
        data = Dataset(f, 'r', format="NETCDF4_CLASSIC")
        if i==0:
            X = np.array(data.variables[var])
            yvals = np.array(data.variables['XLAT'])
            xvals = np.array(data.variables['XLONG'])
        else:
            X = np.concatenate( (X,np.array(data.variables[var])), axis=0 ) 
            yvals = np.concatenate((yvals,np.array(data.variables['XLAT'])),axis=0)
            xvals = np.concatenate((xvals,np.array(data.variables['XLONG'])),axis=0)



    #Shape of X should be time*z*y*x or time*y*x. We need to reshape so that the y and x dimensions are first

    X = np.moveaxis(X, [-2, -1], [0, 1])
    yvals = np.moveaxis(yvals, [-2,-1], [0,1])
    xvals = np.moveaxis(xvals, [-2,-1], [0,1])

    if filt == 2:
        X, Y, yvals, xvals = routines.filt2(X, yvals, xvals, ny, nx) 
    else:
        filt=1
        X, Y, yvals, xvals = routines.filt1(X, yvals, xvals, ny, nx) 
 
    X = np.moveaxis(X, [0, 1], [-2,-1])
    Y = np.moveaxis(Y,  [0, 1], [-2,-1]) 
    yvals = np.moveaxis(yvals,  [0, 1], [-2,-1])
    xvals = np.moveaxis(xvals, [0,1], [-2,-1])

    newfile_name = "Filtered" + str(filt) + "_" + str(ny) + "_" + str(nx) + "_" + var + "_" + date + ".nc"

    routines.deletefile(datapath+newfile_name)

    print("Creating file " + datapath+newfile_name + "...")

    newfile = Dataset(datapath+newfile_name, 'w', format="NETCDF4_CLASSIC")

    newfile.createDimension('south_north', X.shape[-2])
    newfile.createDimension('west_east', X.shape[-1])
    newfile.createDimension('time')

    if len(X.shape)==4:
        newfile.createDimension('interp_levels', X.shape[3])
        oldvar = newfile.createVariable(var,np.float64,('Time','interp_levels','south_north','west_east'))
        newvar = newfile.createVariable(var+'_filtered',np.float64,('time','interp_levels','south_north','west_east')) 
        oldvar[:,:,:,:] = X
        newvar[:,:,:,:] = Y
    else:
        oldvar = newfile.createVariable(var,np.float64,('time','south_north','west_east'))
        newvar = newfile.createVariable(var+'_filtered',np.float64,('time','south_north','west_east'))    
        oldvar[:,:,:] = X
        print(Y.shape, xvals.shape,yvals.shape)
        newvar[:,:,:] = Y

    lons = newfile.createVariable('XLONG',np.float64,('time','south_north','west_east'))
    lats = newfile.createVariable('XLAT',np.float64,('time','south_north','west_east'))
    
    lons[:,:,:] = yvals
    lats[:,:,:] = xvals

    print("Done.")

    return




main()
