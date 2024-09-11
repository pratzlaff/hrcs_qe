import argparse
import numpy as np
from scipy import interpolate
from matplotlib import pyplot as plt

# lambda is also chr(955)
ANGSTROM, LAMBDA = "Åλ"

# for years later than 2012.5, the maximum correction is 2.1%/yr at 0 Ang, and everything from
# 20 AA down to 0 AA linearly transisitions from 1.0 to that maximum
def corr_old(y, wav):
    y0 = 2012.5
    w0, w1 = 0, 20
    rate = 0.021

    max = 1 - (y - y0) * rate

    corr = np.ones_like(wav)
    if y < y0: return corr
    corr[wav<w0] = max
    ii = np.where((wav>=w0) & (wav<w1))
    corr[ii] = 1 - (1-max) * (w1-wav[ii]) / (w1-w0)
    corr = 1/corr
    return corr

def corr(y, wav):
    y0 = 2012.5
    w0, w1 = 0, 20
    rate = 0.021

    max = (y - y0) * rate

    corr = np.ones_like(wav)
    if y < y0: return corr
    corr[wav<w0] = max
    ii = np.where((wav>=w0) & (wav<w1))
    corr[ii] = 1 + max * (w1-wav[ii]) / (w1-w0)
    return corr

def plot(args):
    global ANSTROM, LAMBDA
    
    y = np.arange(5)+2018.
    wav = np.arange(191)/10.+1
    plt.plot((wav[0],wav[-1]), (1, 1), 'k-', label='prior to 2013')
    for i in range(y.size):
        c = corr(y[i], wav)
        plt.plot(wav, c, '--', label='{:g}'.format(y[i]))
    plt.xlabel(r'{} ({})'.format(LAMBDA, ANGSTROM))
    plt.ylabel('QE Multiplicative Factor')
    plt.legend()

    plt.tight_layout()

    if args.outfile:
        plt.savefig(args.outfile)
    else:
        plt.show()

    plt.close()
    exit()


    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(8.5,11))

    fmts = ['--', '--', '--']

    ax = axes[0]

    for i in range(corrs.size):
        ax.plot(years, [1.0, corrs[i]], fmts[i], label=r'{:.3g} Å'.format(wavmid[i]))
    ax.set_xlabel('Date')
    ax.set_ylabel('QE Divided By')
    ax.set_title('Inverse Corrective Factor vs Date')
    ax.legend()

    ax = axes[1]

    y = np.arange(5)+2018.
    wav = np.arange(191)/10.+1
    for i in range(y.size):
        c = corr(y[i], wav)
        ax.plot(wav, c, '--', label='{}'.format(y[i]))
    ax.set_xlabel(r'$\lambda$')
    ax.set_ylabel('Corr')
    ax.set_title('Corrective Factor vs Wavelength')
    ax.legend()

    plt.tight_layout()

    if args.outfile:
        plt.savefig(args.outfile)
    else:
        plt.show()

    plt.close()

def main():
    parser = argparse.ArgumentParser(
        description='Plot corrective curves for HRC-S short wavelengths.'
    )
    parser.add_argument('-o', '--outfile', help='Output file name.')
    args = parser.parse_args()
    plot(args)

if __name__ == '__main__':
    main()
