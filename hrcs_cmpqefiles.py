import astropy.io.fits
import numpy as np
import argparse
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib
matplotlib.rc('text', usetex=True)
from scipy import interpolate

def read_qe(qeufile, hdunum):
    hdulist = astropy.io.fits.open(qeufile)
    data = hdulist[hdunum].data

    regionid = data.field('regionid')
    energy = data.field('energy')
    qe = data.field('qe')

    hdulist.close()

    return regionid, energy, qe

def main():
    parser = argparse.ArgumentParser(
        description='Compare HRC-S QE files.'
    )
    parser.add_argument('file1')
    parser.add_argument('file2')
    args = parser.parse_args()

    for hdunum in range(1,4):

        region1, energy1, qe1 = read_qe(args.file1, hdunum)
        region2, energy2, qe2 = read_qe(args.file2, hdunum)

        assert np.abs(region1-region2).sum()==0
        assert np.abs(energy1-energy2).sum()==0

        plot_dims = (energy1.shape[0], 1)

        for j in range(energy1.shape[0]):

            plt.subplot2grid(plot_dims, (j, 0))

            ratio = qe2[j] / qe1[j]
            wav = 12.39854 / energy2[j]
            plt.plot(wav, ratio)
            #plt.ylim(0.5, 1)

        plt.tight_layout()
        plt.show()
        plt.close()



if __name__ == '__main__':
    main()
