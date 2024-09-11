import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits

def main():

    hdulist = astropy.io.fits.open('hrcsD2012-03-29qeN0014.fits')

#    plot_dims = (4, 2)
    plot_dims = (1, 1)

    k=0
#    for i in xrange(3):
    for i in xrange(1):

        data = hdulist[i+1].data

#        for j in xrange(data.field('energy').shape[0]):
        for j in xrange(1):
            row = int(k/plot_dims[1])
            col = k % plot_dims[1]
            plt.subplot2grid(plot_dims, (row, col))

            plt.plot(data.field('energy')[j], data.field('qe')[j])
            plt.xlim(0.28, 0.31)
            plt.xlabel('Energy (keV)')
            plt.ylabel('QE')
#            plt.title('k = {}'.format(k))
            plt.title('HRC-S')

            k += 1

    plt.show()
    hdulist.close()

if __name__ == '__main__':
    main()
