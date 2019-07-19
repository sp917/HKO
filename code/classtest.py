import routines
import wrf

class DATA:
    def __init__(self, str_id):

        self.X = routines.getdata(str_id)

        self.yvals = routines.getdata('lat')
        self.xvals = routines.getdata('lon')

        self.Xsmooth = None

    def getsmooth(self, passes):

        self.Xsmooth = wrf.smooth2d(self.X, passes)

        return
