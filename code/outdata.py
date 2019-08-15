import wrf
import numpy as np
import os
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import pandas as pd
import xarray as xa

plotpath = '/mnt/c/Users/hko/Desktop/plots/'
datapath = '/home/sp917/data/'

######################################################################################
######################################################################################
######################################################################################

class DERIVE:

    def __init__(self):

        derived_list = ['ke', 'speed', 'ke10', 'speed10']
        
        function_switch = { 'ke' : [self.getke, ['ua', 'va'] ]
                            'speed' : [self.getspeed, ['ua', 'va'] ]
                            }

       

    def getke(u,v):
        return 0.5*(u*u + v*v)
 
    def getspeed(u,v):
        return np.sqrt(u*u + v*v)


class DATA:
   
    def __init__(self, str_id):

        self.dependencies
        self.D = DERIVE()

    def get_dependencies(self):
        
        if self.str_id in self.D.derived_list:
            
        
    
    
