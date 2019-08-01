import routines
import sys

def main():
   
    var, passes, mu, alt, height, interp = getsysargs()

    X, Xsmooth, dX, yvals, xvals = routines.smooth_and_plot(var,passes,interp,height)
    
    dy,dx = routines.deltayx(yvals,xvals)

    Splot, K = routines.spectrum_and_plot(X, Xsmooth, dX, dy, dx, var, passes, height, interp)
    
    Xhy = routines.add_turb_plot(X,Xsmooth,dX, yvals, xvals, var, passes, mu, alt, height, interp)

    Splots = routines.spect_turb_plot(Xhy, Splot, K, dy, dx, var, passes, mu, alt, height, interp)

    return 0


def getsysargs():
    
    var = sys.argv[1]      #str, specifying variable
    passes = int(sys.argv[2]) #number of times to apply smoothing filter
    mu = float(sys.argv[3])
    alt = int(sys.argv[4])
    try:
        height = float(sys.argv[5])    # height in km
    except:
        if var.endswith('10'):
            height=0.01
        else:
            height=False
        interp = False
    else:
        interp=True


    return var, passes, mu, alt, height, interp

main()
