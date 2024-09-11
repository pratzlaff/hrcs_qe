import numpy as np
import astropy.io.fits
import argparse
import glob

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from matplotlib import rc, rcParams
rc('text', usetex=True)
rcParams.update({'font.size': 18})

def read_rmf(rmf):
    hdulist = astropy.io.fits.open(rmf)
    data = hdulist['matrix'].data

    f_chan = data.field('f_chan')
    n_chan = data.field('n_chan')
    matrix = data.field('matrix')

    n = f_chan.shape[0]

    resp = np.zeros(n)

    for i in xrange(n):
        resp[f_chan[i]-1:f_chan[i]-1+n_chan[i]] += matrix[i]

    energy = 0.5 * (data.field('energ_lo') + data.field('energ_hi'))

    hdulist.close()

    return resp, energy

def read_garf(garf):
    hdulist = astropy.io.fits.open(garf)
    data = hdulist['specresp'].data

    specresp = data.field('specresp')
    energy = 0.5 * (data.field('energ_lo') + data.field('energ_hi'))

    hdulist.close()

    return specresp, energy

def garf_file(obsid, cal, order):
    if order == 1: ostr = '1'
    elif (order == -1): ostr = '-1'
    else: print "whoops" ; exit(1)
    return glob.glob('/data/legs/rpete/flight/xcal_hrcsi/data/hz43/'+obsid+'/tg_reprocess_qe'+cal+'/hrcf*'+obsid+'N*LEG_'+ostr+'_garf.fits')[0]

def rmf_file(obsid, cal, order):
    if order == 1: ostr = '1'
    elif (order == -1): ostr = '-1'
    else: print "whoops" ; exit(1)
    return '/data/legs/rpete/flight/xcal_hrcsi/data/hz43/'+obsid+'/tg_reprocess_qe'+cal+'/LEG_'+ostr+'.rmf'

def main():

    parser = argparse.ArgumentParser(
        description='Plot N0014 / N0013 response'
    )
    parser.add_argument('-p', '--pdf', help='Save plot to named file.')
    parser.add_argument('obsid')
    args = parser.parse_args()

    data = { 'N0013' : { }, 'N0014' : { } }
    orders = { 'pos' : 1, 'neg' : -1 }
    labels = { 'pos' : 'TG\_M = +1', 'neg' : 'TG\_M = -1'}
    colors = { 'pos' : '-k', 'neg' : '-r' }

    for cal in data:
        for order in orders:

            data[cal][order] = { }

            rmf, e_rmf = read_rmf(rmf_file(args.obsid, cal, orders[order]))
            data[cal][order]['rmf'] = rmf
            data[cal][order]['e_rmf'] = e_rmf

            garf, e_garf = read_garf(garf_file(args.obsid, cal, orders[order]))
            data[cal][order]['garf'] = garf
            data[cal][order]['e_garf'] = e_garf

            data[cal][order]['response'] = garf * rmf
            data[cal][order]['lambda'] = 12.39854 / e_garf

    for order in orders:
        outfile = args.obsid + '_garf_rmf_' + order + '.txt'
        np.savetxt(outfile,
                   np.transpose([
                       data['N0013'][order]['e_garf'],
                       data['N0013'][order]['garf'],
                       data['N0014'][order]['garf'],
                       data['N0013'][order]['rmf'],
                       data['N0014'][order]['rmf']]),
                   fmt=['%6g', '%4g', '%4g', '%4g', '%4g'],
                   header='energy\tN0013_garf\tN0014_garf\tN0013_rmf\tN0014_rmf',
                   delimiter='\t')

    if args.pdf:
        pdf = PdfPages(args.pdf)
        fig = plt.figure(figsize=(11, 8.5))

    for order in orders:
        wav = data['N0013'][order]['lambda']
        ratio = data['N0014'][order]['response'] / data['N0013'][order]['response']

        plt.plot(wav, ratio, colors[order], label=labels[order])
        
    plt.legend(loc='lower right', fontsize=10)
    plt.ylabel(r'$\textrm{N0014 / N0013 Response Ratio}$')
    plt.xlabel(r'$\textrm{Wavelength (\AA)}$')
    plt.ylim(.994,1.006)

    plt.tight_layout()

    if args.pdf:
        pdf.savefig(fig)
        pdf.close()

    else:
        plt.show()

if __name__ == '__main__':
    main()

