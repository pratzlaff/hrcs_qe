import numpy as np
import astropy.io.fits
import argparse
import glob
import os
import response

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from matplotlib import rc, rcParams
rc('text', usetex=True)
rcParams.update({'font.size': 18})

def main():

    parser = argparse.ArgumentParser(
        description='Plot ratio of responses'
    )
    parser.add_argument('-p', '--pdf', help='Save plot to named file.')
    parser.add_argument('obsid')
    args = parser.parse_args()

    rmfpaths = { 'N0013' : '/data/legs/rpete/flight/rmfs/HRC_lsfparm_N0003',
                 'N0014' : '/data/legs/rpete/flight/rmfs',
                 }

    data = { 'N0013' : { }, 'N0014' : { } }
    arms = { 'pos' : 1, 'neg' : -1 }
    labels = { 'pos' : 'TG\_M = +1', 'neg' : 'TG\_M = -1'}
    colors = { 'pos' : '-k', 'neg' : '-r' }

    for cal in data:

        os.environ['RMFPATH'] = rmfpaths[cal]
        os.environ['ARFPATH'] = './garfs/S_qe_' + cal

        bin_lo, bin_hi, resp = response.get_response(args.obsid, 'HRC-S', 'LEG', maxorder=1)
        wav = 0.5 * (bin_lo + bin_hi)

        for arm in arms:
            data[cal][arm] = { }
            data[cal][arm]['response'] = resp[arm][0]
            data[cal][arm]['lambda'] = wav[0]

    #         rmf, e_rmf = read_rmf(rmf_file(args.obsid, cal, orders[order]))
    #         data[cal][order]['rmf'] = rmf
    #         data[cal][order]['e_rmf'] = e_rmf

    #         garf, e_garf = read_garf(garf_file(args.obsid, cal, orders[order]))
    #         data[cal][order]['garf'] = garf
    #         data[cal][order]['e_garf'] = e_garf

    #         data[cal][order]['response'] = garf * rmf
    #         data[cal][order]['lambda'] = 12.39854 / e_garf

    # for arm in arms:
    #     outfile = args.obsid + '_garf_rmf_' + order + '.txt'
    #     np.savetxt(outfile,
    #                np.transpose([
    #                    data['N0013'][arm]['e_garf'],
    #                    data['N0013'][arm]['garf'],
    #                    data['N0014'][arm]['garf'],
    #                    data['N0013'][arm]['rmf'],
    #                    data['N0014'][arm]['rmf']]),
    #                fmt=['%6g', '%4g', '%4g', '%4g', '%4g'],
    #                header='energy\tN0013_garf\tN0014_garf\tN0013_rmf\tN0014_rmf',
    #                delimiter='\t')

    if args.pdf:
        pdf = PdfPages(args.pdf)
        fig = plt.figure(figsize=(11, 8.5))

    for arm in arms:
        wav = data['N0013'][arm]['lambda']
        ratio = data['N0014'][arm]['response'] / data['N0013'][arm]['response']

        plt.plot(wav, ratio, colors[arm], label=labels[arm])
        
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

