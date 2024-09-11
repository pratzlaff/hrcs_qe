import argparse
import astropy.io.fits
import glob as gl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import re

def read_obsid_det_date(f):
    hdulist = astropy.io.fits.open(f)
    header = hdulist[1].header
    obsid = header['obs_id']
    detnam = header['detnam']
    mjd = header['mjd-obs']

    hdulist.close()

    date = 1998 + (mjd-50814)/365.2422

    return obsid, detnam, date

def get_pifiles_info(files):
    obsids, dets, dates, srcs = [], [], [], []
    for f in files:
        obsid, det, date = read_obsid_det_date(f)
        obsids.append(obsid)
        dets.append(det)
        dates.append(date)
        srcs.append(re.findall('/\d+_(\w+)_0001', f)[0])
    return obsids, dets, dates, srcs

def pifile_from_fluxfile(fluxfile):
    return fluxfile[:-10]+'_0001.pi'

def get_files(globs):
    fluxfiles = []
    pifiles = []
    for glob in globs:
        fluxfiles_ = gl.glob(glob)
        pifiles_ = [pifile_from_fluxfile(f) for f in fluxfiles_]
        fluxfiles.extend(fluxfiles_)
        pifiles.extend(pifiles_)
    return fluxfiles, pifiles

def mk_data_struct(fluxfiles, obsids, dets, dates, srcs):
    data = { }
    for det in dets:
        data[det] = { }
        for src in srcs:
            data[det][src] = {
                'obsid' : [],
                'flux' : [],
                'errlo' : [],
                'errhi' : [],
                'date' : []
            }

    for i in range(len(fluxfiles)):
        fluxfile = fluxfiles[i]
        obsid = obsids[i]
        det = dets[i]
        date = dates[i]
        src = srcs[i]
        try:
            flux, errlo, errhi = read_fluxes(fluxfile)
            data[det][src]['obsid'].append(obsid)
            data[det][src]['flux'].append(flux)
            data[det][src]['errlo'].append(errlo)
            data[det][src]['errhi'].append(errhi)
            data[det][src]['date'].append(date)
        except:
            pass
    for det in data:
        for src in data[det]:
            for key in data[det][src]:
                data[det][src][key] = np.array(data[det][src][key])

    return data

def read_fluxes(fluxfile):
    hdulist = astropy.io.fits.open(fluxfile)
    data = hdulist[1].data
    flux = data.field('net_mflux_aper')[0]
    errlo = flux - data.field('net_mflux_aper_lo')[0]
    errhi = data.field('net_mflux_aper_hi')[0] - flux
    hdulist.close()
    return flux, errlo, errhi

def main():
    parser = argparse.ArgumentParser(
        description='Plot fluxes from HRC-S G21.5-0.9 observatinos, for two CalDB versions.'
    )
    parser.add_argument('-o', '--outfile', help='Output file name.')
    args = parser.parse_args()

    plt.ylabel('Flux $(\mathregular{erg\;cm^{-2}\;s^{-1}})$')
    plt.xlabel('Date')

    qes = ['N0014', 'N0015']
    labels = { 'N0014' : 'QE N0014', 'N0015' : 'QE N0015' }
    fmts = { 'N0014' : 'o', 'N0015' : 's' }

    globs = [
        '/data/legs/rpete/flight/g21.5-0.9/srcflux/qe_N0014_qeu_N0013/*.flux',
        '/data/legs/rpete/flight/g21.5-0.9/srcflux/qe_N0015_qeu_N0013/*.flux',
    ]

    for qe in qes:
        glob = '/data/legs/rpete/flight/g21.5-0.9/srcflux/qe_{}_qeu_N0013/*.flux'.format(qe)
        print(glob)

        fluxfiles, pifiles = get_files((glob,))
        obsids, dets, dates, srcs = get_pifiles_info(pifiles)
        sorti = sorted(range(len(dates)), key=lambda ix: dates[ix])

        fluxfiles = [fluxfiles[i] for i in sorti]
        pifiles = [pifiles[i] for i in sorti]
        obsids = [obsids[i] for i in sorti]
        dets = [dets[i] for i in sorti]
        dates = [dates[i] for i in sorti]
        srcs = [srcs[i] for i in sorti]

        data = mk_data_struct(fluxfiles, obsids, dets, dates, srcs)

        det, src = 'HRC-S', 'plerion'
        x = data[det][src]['date']
        y = data[det][src]['flux']
        errlo = data[det][src]['errlo']
        errhi = data[det][src]['errhi']

        mask = x > 1998
        if qe == 'N0015':
            mask = x > 2010

        plt.errorbar(x[mask], y[mask],
                     yerr=[errlo[mask], errhi[mask]],
                     fmt=fmts[qe], label=labels[qe])

        if qe == 'N0014':
            mask = x < 2010
            mean = np.mean(y[mask])
            std = np.std(y[mask])
            plt.gca().axhline(mean, color='k', ls='--')
            eb = plt.errorbar((0.5*(x[0]+x[-1]),), (mean,), (std,), fmt='', ecolor='k', capsize=12)
            eb[-1][0].set_linestyle('--')

    plt.legend(loc='upper left')
    plt.tight_layout()
    if (args.outfile):
        plt.savefig(args.outfile)
    else:
        plt.show()

    plt.close()

if __name__ == '__main__':
    main()
