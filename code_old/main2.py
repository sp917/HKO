import classtest
import sys

def main():
    
    print("\nReading arguments.")
    str_id, h, passes, a, b, c, d, e, mu, alt = getargs()
    
    print("\nGetting data.")
    data = classtest.DATA(str_id)

    data.set_height(h)

    print("\nSmoothing " + str(passes) + " times.")
    data.set_smoothness(passes)

    print("\nCalculating spectra.")
    data.calc_spectra()

    print("\nSetting scaling parameters.")
    data.set_parameters(a,b,c,d,e,mu)

    print("\nAdding turbulence from final state to initial state.")
    data.add_turbulence(alt)

    print("\nCalculating the spectrum of the new field.")
    data.spectrum_turb_added()

    return

def getargs():

    str_id = sys.argv[1]
    h = float(sys.argv[2])
    passes = int(sys.argv[3])
    a = float(sys.argv[4])
    b = float(sys.argv[5])
    c = float(sys.argv[6])
    d = float(sys.argv[7])
    e = float(sys.argv[8])
    mu = float(sys.argv[9])
    alt = int(sys.argv[10])

    return str_id, h, passes, a, b, c, d, e, mu, alt

main()
