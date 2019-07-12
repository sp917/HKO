import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import sys
import os

def main():

    print("running python script")
    
    try:
        var = sys.argv[1]
        year = '%0.4d' % int(sys.argv[2])
        month = '%0.2d' % int(sys.argv[3])
        day = '%0.2d' % int(sys.argv[4])
        hour = '%0.2d' % int(sys.argv[5])
    except:
        var = 'speed10'
        year = '%0.4d' %  2019
        month = '%0.2d' % 9
        day = '%0.2d' % 7
        hour = '%0.2d' % 8

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
        else:
            X = np.concatenate( (X,np.array(data.variables[var])), axis=0 )
    print(X.shape)

    #Shape of X should be time*z*y*x or time*y*x. We need to reshape so that the y and x dimensions are first

    X = np.moveaxis(X, [-2, -1], [0, 1])
    print(X.shape)


            
        
        

main()
