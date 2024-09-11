import argparse
import pprint
import matplotlib.pyplot as plt
import numpy as np
import re

import xcal

from matplotlib import rc, rcParams
#rcParams.update({'font.size': 12})

def plot(args):

    hdr, data = xcal.read_fits_rdbs(args.rdbs)

    groups = xcal.group_days(data['days'])

    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(11, 8.5))

    try:
        wlo, whi = xcal.band_wav_limits(int(hdr['band']))
    except:
        wlo, whi = hdr['wlo'], hdr['whi']
        
    coords = { 5 : (0, 0),
               6 : (0, 1),
               7 : (1, 0),
               8 : (1, 1),
    }

    for i in 5, 6, 7, 8:

        date_obs = xcal.date_obs_obsid(data['obsid'][groups[i]][0])

        if np.unique(data['inst'][groups[i]]).size == 1:
            continue

        row, col = coords[i]
        ax = axs[row, col]

        plot_group(hdr, data, groups[i], args, ax)
        if row == 1:
            ax.set_xlabel('Time Offset (ks)')
        if col == 0:
            ylabel = r'{:.2g}-{:.2g} Å Flux'.format(float(wlo), float(whi))
            ax.set_ylabel(ylabel)

        ax.set_title(date_obs[:10])
        if row==0 and col==0:
            ax.legend(loc='upper left', fontsize=10)
            if args.caldb:
                ax.set_title('QE {}'.format(args.caldb), loc='right')

    plt.tight_layout()
    if args.outfile:
        plt.savefig(args.outfile)
    else:
        plt.show()

def plot_group(hdr, data, group, args, ax):

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
    parser.add_argument('-c', '--caldb', help='QE CalDB string.')
    parser.add_argument('--noerrb', action='store_true', help='Do not plot error bars.')
    parser.add_argument('--noheg', action='store_true', help='Do not plot HEG fluxes.')
    parser.add_argument('--nomeg', action='store_true', help='Do not plot MEG fluxes.')
    parser.add_argument('rdbs', nargs='+', help='RDB files containing fit results')
    args = parser.parse_args()

    plot(args)

if __name__ == '__main__':
    main()
