import argparse
import astropy.io.fits
import glob
import matplotlib.pyplot as plt
import numpy as np
import sys

import response
import util
import hz43
import flux

def wav_ranges(detnam):

    # wavelength intervals will be (e.g.):
    # TG_M = -1 : 60-160, 60-65, 65-70, ..., 155-160
    # TG_M = +1 : 70-170, 70-75, 75-80, ..., 165-170

    # FIXME: starting with N0011, obsid 1011 negative-order residuals
    # are not being calculated, and wavelength range doesn't seem to
    # be the culprit

    increment = {
        'HRC-S' : { 'neg' : 2, 'pos' : 2 },
        'HRC-I' : { 'neg' : 10, 'pos' : 10 },
    }

    steps = {
        'HRC-S' : { 'neg' : 50, 'pos' : 52 },
        'HRC-I' : { 'neg' : 1, 'pos' : 1 },
    }

    minimum = {
        'HRC-S' : { 'neg' : 57, 'pos' : 69 },
        'HRC-I' : { 'neg' : 60, 'pos' : 60 },
    }

    w1 = { order : (np.arange(steps[detnam][order]+1)-1)*increment[detnam][order]+minimum[detnam][order]
           for order in increment['HRC-S'] }

    w2 = {}
    for order in w1:
        w2[order] = w1[order] + increment[detnam][order]
        w1[order][0] = w1[order][1]
        w2[order][0] = w2[order][-1]

    return w1, w2


def dispersed_flux(src, bg, resp, bin_lo, bin_hi, wav_lo, wav_hi, hdr, predicted, factor, args):
    ind = np.where((bin_lo>=wav_lo) & (bin_hi<wav_hi))
    rate, rate_err, flux_, flux_err = flux.flux_summed(src[ind], bg[ind], resp[ind], bin_lo[ind], bin_hi[ind], hdr, factor[ind])

    ratio = rate / predicted[ind].sum()
    ratio_err = rate_err / predicted[ind].sum()

    if args.flux:
        ratio = flux_ / predicted[ind].sum()
        ratio_err = flux_err / predicted[ind].sum()

    return rate, rate_err, flux_, flux_err, ratio, ratio_err

def dispersed_lc_obsid(obsid, w1, w2, tg_reprocess, args):

    bin_lo, bin_hi, resp = response.get_response(obsid, 'LEG', maxorder=args.maxorder)
    wav = 0.5*(bin_lo+bin_hi)
    energy = util.w2e(wav)

    if 'model_flux' not in dispersed_lc_obsid.__dict__:
        model_wav, model_flux = hz43.model()
        dispersed_lc_obsid.model_flux = np.interp(wav, model_wav, model_flux) * (bin_hi - bin_lo)

    # FIXME: figure out how python can assign this variable without this assignment!
    model_flux = dispersed_lc_obsid.model_flux

    # ratio of 1st-order counts / all orders
    predicted_rate = { order : model_flux / 1.602e-9 / energy * resp[order] for order in resp }
    maxi = predicted_rate['pos'].shape[0]
    factor = { order : predicted_rate[order][0] / predicted_rate[order].sum(axis=0) for order in resp }

    rates, rate_errs, fluxes, flux_errs, ratios, ratio_errs = ({} for  i in range(6))
    for order in w1:
        rates[order] = np.zeros((w1[order].size,))
        rate_errs[order] = rates[order].copy()
        fluxes[order] = rates[order].copy()
        flux_errs[order] = rates[order].copy()
        ratios[order] = rates[order].copy()
        ratio_errs[order] = rates[order].copy()

    # read PHA2
    d, h = util.read_pha2(util.pha2_file(obsid, tg_reprocess=tg_reprocess))

    rows = {}
    orders = { 'neg':-1, 'pos':+1 }
    for order in w1:

        ind = np.where(d['tg_m']==orders[order])[0][0]
        src = d['counts'][ind]
        bg = d['background_up'][ind] + d['background_down'][ind]
            
        predicted = predicted_rate[order][0]
        if args.flux:
            predicted = model_flux[0]
            
        for j in range(w1[order].size):
            r, rerr, f, ferr, ra, raerr = dispersed_flux(src, bg, resp[order][0], bin_lo[0], bin_hi[0], w1[order][j], w2[order][j], h, predicted, factor[order], args)

            rates[order][j] = r
            rate_errs[order][j] = rerr
            fluxes[order][j] = f
            flux_errs[order][j] = ferr
            ratios[order][j] = ra
            ratio_errs[order][j] = raerr

    return h['date-obs'][:10], rates, rate_errs, fluxes, flux_errs, ratios, ratio_errs

# get HRC-S/LETG counts light curves for dispersed orders
def dispersed_lc(args, detnam, tg_reprocess='tg_reprocess'):
    obsids, years = hz43.obsids_years(detnam)
    w1, w2 = wav_ranges(detnam)
    date_str = []

    rates, rate_errs, fluxes, flux_errs, ratios, ratio_errs = ({} for  i in range(6))
    for order in w1:
        rates[order] = np.zeros((w1[order].size, years.size))
        rate_errs[order] = rates[order].copy()
        fluxes[order] = rates[order].copy()
        flux_errs[order] = rates[order].copy()
        ratios[order] = rates[order].copy()
        ratio_errs[order] = rates[order].copy()

    for i in range(obsids.size):
        obsid = obsids[i]
        dstr, r, rerr, f, ferr, ra, raerr = dispersed_lc_obsid(obsid, w1, w2, tg_reprocess, args);
        date_str.append(dstr)

        for order in w1:
            rates[order][:,i] = r[order]
            rate_errs[order][:,i] = rerr[order]
            fluxes[order][:,i] = f[order]
            flux_errs[order][:,i] = ferr[order]
            ratios[order][:,i] = ra[order]
            ratio_errs[order][:,i] = raerr[order]

    return obsids, years, date_str, w1, w2, rates, rate_errs, fluxes, flux_errs, ratios, ratio_errs

# get HRC-S/LETG counts light curves for dispersed orders
def dispersed_lc_old(args, detnam, tg_reprocess='tg_reprocess'):

    orders = { 'neg':-1, 'pos':+1 }

    obsids, years = hz43.obsids_years(detnam)
    w1, w2 = wav_ranges(detnam)

    model_flux = None

    date_str = []

    rates, rate_errs, fluxes, flux_errs, ratios, ratio_errs = ({} for  i in range(6))
    for order in w1:
        rates[order] = np.zeros((w1[order].size, years.size))
        rate_errs[order] = rates[order].copy()
        fluxes[order] = rates[order].copy()
        flux_errs[order] = rates[order].copy()
        ratios[order] = rates[order].copy()
        ratio_errs[order] = rates[order].copy()

    for i in range(obsids.size):
        obsid = obsids[i]

        bin_lo, bin_hi, resp = response.get_response(obsid, 'LEG', maxorder=args.maxorder)

        if model_flux is None:
            model_wav, model_flux = hz43.model()
            wav = 0.5*(bin_lo+bin_hi)
            energy = util.w2e(wav)
            model_flux = np.interp(wav, model_wav, model_flux) * (bin_hi - bin_lo)

        # ratio of 1st-order counts / all orders
        predicted_rate = { order : model_flux / 1.602e-9 / energy * resp[order] for order in resp }
        maxi = predicted_rate['pos'].shape[0]
        factor = { order : predicted_rate[order][0] / predicted_rate[order].sum(axis=0) for order in resp }

        # read PHA2
        d, h = util.read_pha2(util.pha2_file(obsid, tg_reprocess=tg_reprocess))

        date_str.append(h['date-obs'][:10])

        rows = {}
        for order in w1:

            ind = np.where(d['tg_m']==orders[order])[0][0]
            src = d['counts'][ind]
            bg = d['background_up'][ind] + d['background_down'][ind]
            
            predicted = predicted_rate[order][0]
            if args.flux:
                predicted = model_flux[0]
            
            for j in range(w1[order].size):
                r, rerr, f, ferr, ra, raerr = dispersed_flux(src, bg, resp[order][0], bin_lo[0], bin_hi[0], w1[order][j], w2[order][j], h, predicted, factor[order], args)
                
                rates[order][j][i] = r
                rate_errs[order][j][i] = rerr
                fluxes[order][j][i] = f
                flux_errs[order][j][i] = ferr
                ratios[order][j][i] = ra
                ratio_errs[order][j][i] = raerr

    return obsids, years, date_str, w1, w2, rates, rate_errs, fluxes, flux_errs, ratios, ratio_errs

# get HRC-S/LETG counts light curves for zeroth order
def zeroth_lc(detector, tg_reprocess='tg_reprocess'):
    if (detector == 'HRC-S'):
        obsids, years = hz43.obsids_years('HRC-S')
    elif (detector == 'HRC-I'):
        obsids, years = hz43.obsids_years('HRC-I')
    else:
        raise ValueError(det)

    rates, rate_errs = util.zeroth_rates(obsids, tg_reprocess=tg_reprocess)
    model_rates = hz43.predicted_rates(obsids)
    return years, rates, rate_errs, model_rates, rates/model_rates, rate_errs/model_rates

def plot_zero(d, args, label=None):
    x = d['year']
    y = d['rate']
    yerr = d['rate_err']

    y_2012_5 = np.interp(2012.5, x, y)
    y = y/y_2012_5
    yerr = yerr/y_2012_5

    if args.ratio:
        y = d['ratio']
        yerr = d['ratio_err']

    mask = x > 2012
    mask = x > 0

    plt.errorbar(x[mask], y[mask], yerr[mask], label=label)

def plot_dispersed(d, detnam, order, index, args, color=None):
    label_prefix = { 'HRC-S' : 'S - ', 'HRC-I' : 'I - ' }
    labels = { 'pos' : '+1st Order', 'neg' : '-1st Order' }

    linestyles = { 'pos' : '-', 'neg' : '--' }

    x = d['year']

    y = d['rate'][order][index]
    yerr = d['rate_err'][order][index]

    y_2012_5 = np.interp(2012.5, x, y)
    y = y/y_2012_5
    yerr = yerr/y_2012_5

    if args.ratio:
        y = d['ratio'][order][index]
        yerr = d['ratio_err'][order][index]

    mask = x > 2012
    mask = x > 0

    label = label_prefix[detnam]+labels[order] + ": {:.0f}-{:.0f} Å".format(d['bin_lo'][order][index], d['bin_hi'][order][index])
    plt.errorbar(x[mask], y[mask], yerr[mask], label=label, color=color, linestyle=linestyles[order])

def linest(x, y, yerr):

    y_2008_5 = np.interp(2008.5, x, y)
    y = y/y_2008_5
    yerr = yerr/y_2008_5

    mask = (x > 2012) & (x < 2020)
    mask = x > 0
    x = x[mask]
    y = y[mask]
    yerr = yerr[mask]

    p, cov = np.polyfit(x, y, 1, w=yerr, cov=True)

    print('b = {} +/- {}\nm = {} +/- {}'.format(
        p[1], np.sqrt(cov[1,1]),
        p[0], np.sqrt(cov[0,0]),
        ))

def main():
    parser = argparse.ArgumentParser(
        description='HRC/LETG HZ 43 observations',
    )
    parser.add_argument('--tg_reprocess_hrcs', default='tg_reprocess', help='tg_reprocess output directory for HRC-S.')
    parser.add_argument('--tg_reprocess_hrci', default='tg_reprocess', help='tg_reprocess output directory for HRC-I.')
    parser.add_argument('-o', '--outfile', help='Output file name.')
    parser.add_argument('-f', '--flux', help='Compare with model fluxes, rather than rates.', action='store_true')
    parser.add_argument('--noi', help='Do not plot HRC-I zeroeth order.', action='store_true')
    parser.add_argument('-r', '--ratio', help='Plot observed/predicted ratios rather than rates.', action='store_true')
    parser.add_argument('-m', '--maxorder', help='Maximum ARF/RMF order to read.', default=3, type=int)

    args = parser.parse_args()

    s_lc_0 = {}
    s_lc_0.update(zip(('year', 'rate', 'rate_err', 'model_rate', 'ratio', 'ratio_err'), zeroth_lc('HRC-S', tg_reprocess=args.tg_reprocess_hrcs)))

    linest(s_lc_0['year'],
           s_lc_0['rate'],
           s_lc_0['rate_err'],
           )

    s_lc_1 = {}
    s_lc_1.update(zip(('obsid', 'year', 'date', 'bin_lo', 'bin_hi', 'rate', 'rate_err', 'flux', 'flux_err', 'ratio', 'ratio_err'), dispersed_lc(args, 'HRC-S', tg_reprocess=args.tg_reprocess_hrcs)))

    if not args.noi:
        i_lc_0 = {}
        i_lc_0.update(zip(('year', 'rate', 'rate_err', 'model_rate', 'ratio', 'ratio_err'), zeroth_lc('HRC-I', tg_reprocess=args.tg_reprocess_hrci)))
        plot_zero(i_lc_0, args, label=r'I - 0th Order')

        i_lc_1 = {}
        i_lc_1.update(zip(('obsid', 'year', 'date', 'bin_lo', 'bin_hi', 'rate', 'rate_err', 'flux', 'flux_err', 'ratio', 'ratio_err'), dispersed_lc(args, 'HRC-I', tg_reprocess=args.tg_reprocess_hrci)))
        for order in i_lc_1['rate']:
            plot_dispersed(i_lc_1, 'HRC-I', order, 0, args)

        obsids, years = hz43.obsids_years('HRC-I', offaxis=True)
        ii = np.where((obsids==1516)|(obsids==1004)|(obsids==1005)|(obsids==2601)|(obsids==2603)|(obsids==18415))
        print('WHEE', obsids)
        obsids, years = obsids[ii], years[ii]
        w1 = {
            1516  : { 'pos' : np.array((65.,)) },
            1004  : { 'pos' : np.array((65.,)) },
            1005  : { 'neg' : np.array((65.,)) },
            2601  : { 'pos' : np.array((65.,)) },
            2603  : { 'pos' : np.array((65.,)) },
            18415 : { 'pos' : np.array((65.,)) },
        }
        w2 = {
            1516  : { 'pos' : np.array((100.,)) },
            1004  : { 'pos' : np.array((100.,)) },
            1005  : { 'neg' : np.array((100.,)) },
            2601  : { 'pos' : np.array((100.,)) },
            2603  : { 'pos' : np.array((100.,)) },
            18415 : { 'pos' : np.array((95.,)) },
        }
        for i in range(obsids.size):
            jnk, r, rerr, f, ferr, ra, raerr = dispersed_lc_obsid(obsids[i], w1[obsids[i]], w2[obsids[i]], args.tg_reprocess_hrci, args);
            order = list(w1[obsids[i]].keys())[0]
            plt.errorbar(years[i], ra[order][0], raerr[order][0], fmt='ko')

    plot_zero(s_lc_0, args, label=r'S - 0th Order')
    for order in s_lc_1['rate']:
        plot_dispersed(s_lc_1, 'HRC-S', order, 0, args)
    plt.legend()
    if args.ratio:
        plt.ylabel('Count Rate Ratio: Observed / Predicted')
    else:
        plt.ylabel('Count Rate Relative to 2012.5')
    plt.xlabel('Date')

    if args.outfile:
        plt.savefig(args.outfile)
    else:
        plt.show()

if __name__ == '__main__':
    main()
