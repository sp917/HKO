from netCDF4 import Dataset
import wrf
import numpy as np
import os

datapath = "../data/"

def modfields(passes):

    fnew = Dataset("wrfinput_d01_new.nc", 'r+')
    
    filenames = [datapath+f for f in os.listdir(datapath) \
                if f.startswith('wrfout_xhka')]

    f0 = Dataset(filenames[0])
    f2 = Dataset(filenames[-1])

    scalingdict = {
                    'U'      : 0.75,
                    'V'      : 0.75,
                    'U10'    : 0.75,
                    'V10'    : 0.75,
                    'T'      : 0.5,
                    'TH2'    : 1.0,
                    'QVAPOR' : 0.5,
    }

    
    for str_id in ['U', 'V', 'T', 'QVAPOR']:

        Xnew = fnew.variables[str_id][0]

        X0 = wrf.getvar(f0, str_id, meta=False)
        X2 = wrf.getvar(f2, str_id, meta=False)

        filename = "../data/" + str_id + "_final_time_" + str(passes)
        try:
            f = open(filename + '.npy')
            f.close()
        except:    
            print("Smoothing " + str_id + " " + str(passes) + " times...")
            Xsmooth = wrf.smooth2d(X2, passes)
            print("Done.")
            np.save(filename, Xsmooth)
        else:
            print("Reading data from " + filename)
            Xsmooth = np.load(filename + '.npy')
            print('done')
            
        dX = X2 - Xsmooth     

        mu = scalingdict[str_id]

        print('calculating new field...')

        turb = dX/np.max(np.abs(Xsmooth))

        newfield = X0*(1 + 2*mu*turb)
 
        del turb

        Xnew[:,] = 0

        Xnew[:,] = newfield[:,]

        del newfield
  
        print('done')

    fnew.close()

    return

