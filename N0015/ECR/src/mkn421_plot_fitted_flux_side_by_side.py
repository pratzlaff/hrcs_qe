import argparse
import pprint
import matplotlib.pyplot as plt
import numpy as np
import re

import xcal

from matplotlib import rc, rcParams
#rcParams.update({'font.size': 12})

ANGSTROM, LAMBDA = "Åλ"

def plot(args):

    rdbs = { 'n0014' : (args.acis, args.n0014),
             'n0015' : (args.acis, args.n0015),
    }

    data = {}
    for s in rdbs:
        hdr, data[s] = xcal.read_fits_rdbs(rdbs[s])

    # the assumption is that n0014 and n0015 rdb files contain the same groups
    s = 'n0014'

    groups = xcal.group_days(data[s]['days'])

    fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(8.5, 11), sharey='row') #, sharex=True)

    try:
        wlo, whi = xcal.band_wav_limits(int(hdr['band']))
    except:
        wlo, whi = hdr['wlo'], hdr['whi']
        
    cols = { 'n0014' : 0, 'n0015' : 1 }

    #for i in 6, 7, 8:

    #    row = i - 6
    for i in -3, -2, -1:

        row = 6 + i - 3

        date_obs = xcal.date_obs_obsid(data[s]['obsid'][groups[i]][0])

        if np.unique(data[s]['inst'][groups[i]]).size == 1:
            continue

        for s in data:
            col = cols[s]
            ax = axs[row, col]

            plot_group(data[s], groups[i], args, ax)

            if col == 0:
                ax.set_title(date_obs[:10])

            if row == 0:
                ax.set_title('QE '+s.upper(), loc='right')

            if row == 1 and col == 0:
                global ANGSTROM
                ylabel = r'{:.2g}-{:.2g} {} Flux $(\mathregular{{erg\;cm^{{-2}}\;s^{{-1}}}})$'.format(float(wlo), float(whi), ANGSTROM)
                ax.set_ylabel(ylabel)

            if row == 2:
                ax.set_xlabel('Time Offset (ksec)')

            if row==0 and col==0:
                ax.legend(loc='upper left', fontsize=10)

    plt.tight_layout()
    if args.outfile:
        plt.savefig(args.outfile)
    else:
        plt.show()

def plot_group(data, group, args, ax):

    doff = data['days'][group][0]
    days_ = data['days'][group] - doff
    dfudge = .003 * days_[-1]

    colors = {
        'l' : 'k',
        'm' : 'r',
        'h' : 'm',
        's' : 'b'
    }

    markers = { 'neg' : 'v' , 'pos' : '^' }
    linestyles = { 'neg' : '--' , 'pos' : '-' } # for the error bars

    param = 'flux'

    days = { 'neg' : days_ - dfudge, 'pos' : days_ + dfudge }
    for inst in colors:

        if inst == 'h' and args.noheg:
            continue

        if inst == 'm' and args.nomeg:
            continue

        for orders in 'neg', 'pos':
            ii = np.where(
                (data['inst'][group]==inst) & (data['order'][group]==orders)
            )[0]

            if not ii.size:
                continue

            best = data[param][group][ii]
            low = data[param+'_min'][group][ii]
            high = data[param+'_max'][group][ii]
                                    
            label=None
            if orders=='pos':
                det = list(xcal.configs[inst].keys())[0]
                part = xcal.configs[inst][det]
                label=det+'/'+part
                        
            yerr=None
            if not args.noerrb:
                yerr = [best-low, high-best]
                yerravg = 0.5*(yerr[1]+yerr[0])
                mask = yerravg < 2*np.median(yerravg) 
                eb = ax.errorbar(days[orders][ii][mask]*86.4,
                                  best[mask],
                                  yerr=(yerr[i][mask] for i in range(len(yerr))),
                                  ecolor=colors[inst],
                                  fmt=colors[inst]+markers[orders],
                                  capthick=0,
                                  label=label
                )

                # eb[-1][0] is the LineCollection objects of the errorbar lines
                eb[-1][0].set_linestyle(linestyles[orders])

            else:
                ax.plot(days[orders][ii]*86.4, best, colors[inst]+markers[orders], label=label)

                
def main():

    parser = argparse.ArgumentParser(
        description='Plot interleaved calibration observation fitted fluxes.'
    )
    parser.add_argument('-o', '--outfile', help='Output file name.')
    parser.add_argument('--noerrb', action='store_true', help='Do not plot error bars.')
    parser.add_argument('--noheg', action='store_true', help='Do not plot HEG fluxes.')
    parser.add_argument('--nomeg', action='store_true', help='Do not plot MEG fluxes.')
    parser.add_argument('acis', help='ACIS RDB file containing fit results')
    parser.add_argument('n0014', help='HRC-S RDB file containing fit results for QE N0014')
    parser.add_argument('n0015', help='HRC-S RDB file containing fit results for QE N0015')
    args = parser.parse_args()

    plot(args)

if __name__ == '__main__':
    main()
