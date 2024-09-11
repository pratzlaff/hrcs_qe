import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
rc('text', usetex=True)

def main():
    wav, neg, pos = np.loadtxt('NewOldRatio.out', usecols=(0,1,2), skiprows=2, unpack=True)
    ratio = 0.5 * (neg + pos)

    plt.plot(wav, ratio)
    plt.xlabel(r'$\textrm{Wavelength (\AA)}$')
    plt.ylabel(r'$\textrm{N0004 / N0003 EEFRAC, default HRC-S/LETG region}$')
    plt.tight_layout()
#    plt.show()
    plt.savefig('hrcs_eefrac_ratio_N0004_N0003.png')

if __name__ == '__main__':
    main()
