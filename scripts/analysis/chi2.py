'''LEE Significance Calculation.

Several functions to compute chi2 significance of a LEE model, and
produce various plots.

The input files and the set of functions that are run are set in the
"if __name__ == '__main__'" code block.

A. Mastbaum <mastbaum@uchicago.edu>, 2017/09/11
'''

import sys
import numpy as np
from numpy.core.umath_tests import inner1d
import scipy.stats
from matplotlib import pyplot as plt

def mc_chi2(data, exp, cov, n=250000, n2=0, plot=False):
    '''Calculate a p-value with a Monte Carlo using a chi2-like statistic.

    :param data: The observed data histogram
    :param exp: The expected (background) histogram
    :param cov: The covariance matrix
    :param n: Number of Monte Carlo draws to evaluate significance
    :param plot: If True, plot the distribution of the test statistic
    :returns: A p-value, for P(data|exp,cov)
    '''
    # Build the MC chi2 distribution, working in chunks to keep the memory
    # usage within reason (to be improved).
    c2ss = []
    for ii in range(100):
        nn = int(1.0 * n / 100)

        # The mean values for n universes, sampled according to the covariance
        means = np.random.multivariate_normal(exp, cov-np.diag(exp), size=nn)
        means[means<0] = 0  # We have to cap these at zero to Poisson sample

        # Fake-data observations for universes, Poisson-distributed
        obs = np.random.poisson(means).astype(np.float32)

        # Inverse covariance matrix
        icov = np.linalg.inv(cov)

        # Calculate the correlated chi2s
        v = obs - exp
        c2s = inner1d(np.dot(v, icov), v)

        c2ss.append(c2s)

    c2ss = np.sort(np.hstack(c2ss))  # Sorted list of all sampled chi2s

    # Optionally plot the MC distribution with a standard chi2
    if plot:
        # Histogram the statistics for the sampled universes
        hist, edges = np.histogram(c2ss, bins=250)
        hist = hist.astype(np.float32) / np.sum(hist) / np.diff(edges)[0]

        # Try to fit this thing to a chi2, floating the dof
        def fcn(dof):
            return np.abs(hist - scipy.stats.chi2.pdf(edges[:-1], df=dof))
        r = scipy.optimize.leastsq(fcn, len(data))
        bfdof = r[0]

        plt.bar(edges[:-1], hist, label='Monte Carlo')
        plt.plot(edges[:-1], scipy.stats.chi2.pdf(edges[:-1], df=len(data)),
                 color='black', label='$\chi^2(dof=N)$')
        plt.plot(edges[:-1], scipy.stats.chi2.pdf(edges[:-1], df=bfdof), '--',
                 color='black', label='$\chi^2(dof=%1.2f)$' % bfdof)
        plt.xlabel('$\chi^2$ Statistic')
        plt.ylabel('Probability/bin')
        plt.xlim(0, 100)
        plt.legend(loc='upper right')
        plt.show()

    # Calculate chi2 for the observed data (or many fake datasets if n2>0)
    if n2 > 0:
        sigs = []
        abv = []

        for ii, d in enumerate(np.random.poisson(np.tile(data, (n2, 1)))):
            v = d - exp
            c2d = inner1d(np.dot(v, icov), v)

            i = np.argmin(np.abs(c2ss - c2d))
            above = len(c2ss[i:])
            p = 1.0 * above / len(c2ss)
            sigma = scipy.stats.norm.ppf(1.0 - p)

            sigs.append(sigma)
            abv.append(above)

        sigma = np.median(sigs)
        above = np.mean(abv)

        print('mc: sigma=%1.5f, above=%i' % (sigma, above))

    else:
        v = data - exp
        c2d = inner1d(np.dot(v, icov), v)

        i = np.argmin(np.abs(c2ss - c2d))
        above = len(c2ss[i:])
        p = 1.0 * above / len(c2ss)
        sigma = scipy.stats.norm.ppf(1.0 - p)

        print('mc: sigma=%1.5f, above=%i' % (sigma, above))

    return sigma


def load(em_file, signal_file):
    '''Load histograms and covariance matrices from disk.

    Nominally these are the output from TSCovariance.

    :param em_file: File with e and mu samples
    :param signal_file: File the the signal sample
    :returns: Tuple of histograms, cov
    '''
    import ROOT

    fi = ROOT.TFile(em_file)
    fs = ROOT.TFile(signal_file)

    henu_e = fi.Get('enu_nue')
    henu_m = fi.Get('enu_numu')
    henu_s = fs.Get('hg')
    hcov = fi.Get('cov')

    # Load raw histograms and covariance matrix
    ebins = range(1, henu_e.GetNbinsX()+1)
    mbins = range(1, henu_m.GetNbinsX()+1)
    em = np.array([henu_m.GetBinContent(i) for i in ebins])
    ee = np.array([henu_e.GetBinContent(i) for i in mbins])

    eb = np.maximum(1e-40, np.hstack((em, ee)))
    cov = np.empty(shape=(len(eb), len(eb)), dtype=eb.dtype)
    for i in range(len(eb)):
        for j in range(len(eb)):
            cov[i][j] = hcov.GetBinContent(i, j)

    # Select desired bins and recompute the covariances
    gmbins = np.arange(3, 17)
    gebins = np.arange(1, len(ee)-1)

    fcov = cov / np.einsum('i,j', eb, eb)
    good = np.hstack((gmbins, len(mbins)+gebins))
    em = em[gmbins]
    ee = ee[gebins]
    eb = eb[good]
    fcov = fcov[good[:,np.newaxis],good]
    cov = fcov * np.einsum('i,j', eb, eb)

    emx = np.array([henu_m.GetBinCenter(int(i)) for i in gmbins])
    eex = np.array([henu_e.GetBinCenter(int(i)) for i in gebins])

    sscale = 1.0 / 10  # Arbitrary scaling to get the signal ish-right
    es = np.array([henu_s.GetBinContent(int(i)) for i in good]) * sscale
    ess = es[len(emx):]

    return emx, em, eex, ee, es, ess, cov


def significance(em, ee, es, fcov,
                 pot=5e19, syst=True, mc=False, eff_m=0.3, eff_e=0.3):
    '''Calculate the significance for the signal, 

    There are two options for evaluating the significance: a standard chi2
    distribution (i.e. assuming that the requirement for Gaussian bin errors
    is satisfied) or by performing a Monte Carlo to build the distribution
    of the test statistic.

    In either case, the resulting p-value is converted to a "sigma" by
    evaluating the corresponding quantile on the normal distribution.

    :param pot: Exposure (POT)
    :param syst: Include systematic uncertainties (or just statistical)
    :param mc: Calculate the p-value with a Monte Carlo
    :param eff_m: 1m1p efficiency
    :param eff_e: 1e1p efficiency
    :returns: Significance of the excess, as (p-value, sigma)
    '''
    sf = pot / 5e19

    # Rescaled efficienies. Magic numbers are the actual selection efficiency
    # (reported by the truth selection code).
    em = em * sf * eff_m / 0.382086
    ee = ee * sf * eff_e / 0.854848
    es = es * sf * eff_e / 0.854848
    eb = np.hstack((em, ee))

    if syst:
        covr = fcov * np.einsum('i,j', eb, eb) + np.diag(eb)
    else:
        covr = np.diag(eb)

    if mc:
        return mc_chi2(eb+es, eb, covr, n=100000000, n2=1000)

    else:
        means = np.random.multivariate_normal(eb, covr-np.diag(eb), size=100000)
        means[means<0] = 0
        v = np.random.poisson(means + es).astype(np.float64) - eb
        c2s = inner1d(np.dot(v, np.linalg.inv(covr)), v)
        ps = scipy.stats.chi2.sf(c2s, len(es))
        sigmas = scipy.stats.norm.ppf(1.0 - ps)

        return np.median(sigmas)


def plot_sigma_v_pot(em, ee, es, fcov,
                     pots, syst=True, mc=False, eff_m=0.3, eff_e=0.3):
    '''Create a plot of significance (sigma) vs. exposure (POT).

    :param pots: An iterable of POTs to evaluate.
    :returns: The pyplot Lines (for further decoration)
    '''
    pots = np.array(pots)
    ys = np.empty_like(pots)
    for i, pot in enumerate(pots):
        args = (em, ee, es, fcov, pot, syst, mc, eff_m, eff_e)
        ys[i] = significance(*args)
        print('%1.1e POT: %1.3f sigma %s' % (pot, ys[i], 'syst' if syst else ''))

    plt.xlabel('Exposure (POT)')
    plt.ylabel('Significance ($\sigma$)')
    plot = plt.plot(pots, ys)
    return plot


def plot_cor(cov, outfile='corr.pdf'):
    '''Plot the correlation matrix.

    :param cov: Covariance matrix
    '''
    cor = cov / np.sqrt(np.einsum('ii,jj->ij', cov, cov))
    plt.matshow(cor, interpolation='nearest', origin='lower', cmap='jet')
    plt.gca().xaxis.tick_bottom()
    plt.colorbar()
    if outfile:
        plt.savefig(outfile)
    plt.show()


def plot_spectra(emx, em, eex, ee, ess, fcov, pot=5e19, eff_m=0.3, eff_e=0.3):
    '''Plot the signal and background spectra.

    :param pot: Exposure (POT)
    :param eff_m: 1m1p efficiency
    :param eff_e: 1e1p efficiency
    '''
    sf = pot / 5e19

    em =  em * sf * eff_m / 0.382086
    ee =  ee * sf * eff_e / 0.854848
    es = ess * sf * eff_e / 0.854848
    eb = np.hstack((em, ee))

    covr = fcov * np.einsum('i,j', eb, eb) #+ np.diag(eb)

    fig, (ax1, ax2) = plt.subplots(1, 2)

    # 1mu1p
    poisi = np.array(scipy.stats.poisson.interval(mu=em, alpha=0.682))
    stat = [em - poisi[0], poisi[1] - em]
    syst = np.sqrt(np.diag(covr[:len(emx),:len(emx)]))
    tot  = [np.sqrt(np.square(stat[0]) + np.square(syst)),
            np.sqrt(np.square(stat[1]) + np.square(syst))]
    ax1.bar(emx, em, width=np.diff(emx)[0], color='gray', label='$1\mu1p$')
    ax1.errorbar(emx, em, stat, fmt=',',
                 color='black', capsize=3, label='Error (stat)')
    ax1.errorbar(emx, em, tot,
                 fmt=',', color='black', label='Error (total)')
    #ax1.errorbar(emx, em,
    #             np.sqrt(np.diag((covr+np.diag(eb))[:len(emx),:len(emx)])),
    #             fmt=',', color='red', label='Error (total)')
    ax1.set_xlabel('$E_\\nu^{CCQE}$ (MeV)')
    ax1.set_ylabel('Events/%1.1f MeV bin/%1.1E POT' % (np.diff(emx)[0], pot))
    ax1.set_xlim(emx[0]-np.diff(emx)[0]/2, emx[-1]+np.diff(emx)[0]/2)
    ax1.set_ylim(0, 1.35 * np.max(em))
    ax1.legend(loc='upper right')
    
    # 1e1p + LEE signal
    poisi = np.array(scipy.stats.poisson.interval(mu=ee, alpha=0.682))
    stat = [ee - poisi[0], poisi[1] - ee]
    syst = np.sqrt(np.diag(covr[len(emx):,len(emx):]))
    tot  = [np.sqrt(np.square(stat[0]) + np.square(syst)),
            np.sqrt(np.square(stat[1]) + np.square(syst))]
    print(stat)
    print(syst)
    print(tot)

    ax2.bar(eex, ee+es, width=np.diff(eex)[0], color='darkred', label='Signal')
    ax2.bar(eex, ee, width=np.diff(eex)[0], color='gray', label='$1e1p$')
    ax2.errorbar(eex, ee, stat, fmt=',',
                 color='black', capsize=3, label='Error (stat)')
    ax2.errorbar(eex, ee, tot,
                 fmt=',', color='black', label='Error (total)')
    #ax2.errorbar(eex, ee,
    #             np.sqrt(np.diag((covr+np.diag(eb))[len(emx):,len(emx):])),
    #             fmt=',', color='red', label='Error (total)', capsize=3)

    import matplotlib.ticker as ticker
    ax2.set_xlabel('$E_\\nu^{CCQE}$ (MeV)')
    #ax2.set_ylabel('Events/%1.1f MeV bin (%1.1E POT)' % (np.diff(eex)[0], pot)
    ax2.set_xlim(eex[0]-np.diff(eex)[0]/2, eex[-1]+np.diff(eex)[0]/2)
    ax2.xaxis.set_major_locator(ticker.MultipleLocator(500))
    ax2.set_ylim(0, 1.1 * np.max(ee+es))
    ax2.legend(loc='upper right')
   
    plt.show()


def plot_sigma_v_eff(em, ee, es, fcov, pot=5e19, mc=False):
    '''Plot significance as a function of mu and e selection efficiencies.

    :param pot: Exposure (POT)
    :param mc: Calculate the p-value with a Monte Carlo
    '''
    eeffs = np.arange(0.1, 0.31, 0.05)
    meffs = np.linspace(0.0001, 0.01, 10)
    see = np.empty(shape=(len(eeffs), len(meffs)), dtype=np.float32)

    for i, effm in enumerate(meffs):
        for j, effe in enumerate(eeffs):
            sigma = significance(em, ee, es, fcov,
                                 pot, mc=mc, eff_m=effm, eff_e=effe)
            print('m %1.2f / e %1.2f: %f' % (effm, effe, sigma))
            see[j][i] = sigma

    xlo1, xhi1 = meffs[0], meffs[-1] + np.diff(meffs)[0]
    xlo2, xhi2 = eeffs[0], eeffs[-1] + np.diff(eeffs)[0]
    plt.matshow(see, interpolation='nearest', origin='lower',
                cmap='jet', extent=(xlo1, xhi1, xlo2, xhi2), aspect='auto')
    plt.gca().xaxis.tick_bottom()
    plt.xlabel('$1\\mu1p$ efficiency')
    plt.ylabel('$1e1p$ efficiency')
    cb = plt.colorbar()
    cb.set_label('Significance ($\sigma$)')
    plt.show()


def constrain_nue(emx, em, eex, ee, ess, fcov, pot=5e19, eff_m=0.3, eff_e=0.3):
    '''Plot the constrained 1e1p systematics for different 1mu1p
    efficiencies.

    This doesn't seem to be working right, statistically. The nue constraint
    for small numu statistics seems too good.
    '''
    sf = pot / 5e19
    ee =  ee * sf * eff_e / 0.854848

    # Unconstrained 1e1p errors
    emr =  em * sf * 1e-10 / 0.382086  # mu efficiency here doesn't matter
    eb = np.hstack((emr, ee))
    covr = fcov * np.einsum('i,j', eb, eb)
    e1 = np.sqrt(np.diag(covr))[len(emr):] / ee
    plt.plot(eex, e1, label='Pre-constraint', ls='-', color='red')

    # Constrained errors, for various 1mu1p efficiencies
    ls = ['--', '-', ':']
    for i, ef in enumerate([1e-3, 1.0]):
        emr =  em * sf * ef / 0.382086
        eb = np.hstack((emr, ee))
        covr = fcov * np.einsum('i,j', eb, eb)
        icov = np.linalg.inv(covr)

        # Total number of mu events in some range
        print('N 1mu1p 200<Eqe<1600 MeV: %f',
              np.sum(emr[(emx>200)&(emx<1600)]))

        # The constraint, a la MiniBooNE (Alexis's thesis)
        ib = np.copy(icov)
        ib[:len(emr),:len(emr)] += 1.0 / emr  # Add 1/N to the muon block
        b = np.linalg.inv(ib)

        errs = np.sqrt(np.diag(b))[len(emr):] / ee
        plt.plot(eex, errs,
                 label='Post-constraint, $\\eta_\\mu=%1.0E$' % ef,
                 ls=ls[i], color='black')

    plt.xlabel('$E_\\nu^\mathrm{CCQE}$ (MeV)')
    plt.ylabel('$1e1p$ fractional systematic uncertainty, %1.1e POT' % pot)
    plt.xlim(0, 3000)
    plt.ylim(0, 0.3)
    plt.axes().yaxis.grid(which='major')
    plt.legend(loc='best')
    plt.show()


def plot_bin_errors(emx, em, eex, ee, ess, fcov,
                    pot=5e19, eff_m=0.3, eff_e=0.3):
    '''Plot the (unconstrained) bin errors.'''
    sf = pot / 5e19
    ee =  ee * sf * eff_e / 0.854848
    em =  em * sf * eff_m / 0.382086
    eb = np.hstack((em, ee))
    covr = fcov * np.einsum('i,j', eb, eb)

    f, (ax2, ax1) = plt.subplots(1, 2, sharey='row', sharex=True)

    ax1.errorbar(eex, np.zeros_like(eex),
                 np.sqrt(np.diag(covr))[len(em):] / ee,
                 fmt=',', ls='-', color='black', label='$1e1p$')

    ax2.errorbar(emx, np.zeros_like(emx),
                 np.sqrt(np.diag(covr))[:len(em)] / em,
                 fmt=',', ls='-', color='black', label='$1\mu1p$')
    ax2.set_ylabel('Fractional systematic uncertainty, %1.1e POT' % pot)

    for ax in [ax1, ax2]:
        ax.set_xlabel('$E_\\nu^\mathrm{CCQE}$ (MeV)')
        #plt.xlim(0, 3000)
        ax.axes.yaxis.grid(which='major')
        ax.legend(loc='upper left')

    plt.show()


if __name__ == '__main__':

    # argv1: merged background output from cov.py
    # argv2: merged signal output from cov.py 
    # Load data from files
    enu_mx, enu_m, enu_ex, enu_e, enu_s, enu_ss, cov = \
        load(sys.argv[1] + '/cov_all.root', sys.argv[2] + '/cov_all.root')

    eb = np.hstack((enu_m, enu_e))
    fcov = cov / np.einsum('i,j', eb, eb)

    # Plot significance vs. POT
    do_sigma_pot_plot = True
    if do_sigma_pot_plot:
        #pots = np.linspace(5e19, 6.6e20, 3)
        pots = np.linspace(5e19, 6.8e20, 10)
        p1 = plot_sigma_v_pot(enu_m, enu_e, enu_s, fcov, pots, mc=False)[0]
        p1.set_label('Stat+Syst')
        p2 = plot_sigma_v_pot(enu_m, enu_e, enu_s, fcov, pots, mc=False, syst=False)[0]
        p2.set_label('Stat Only')
        plt.legend(loc='best')
        plt.savefig('sigma_v_pot.pdf')
        plt.show()

    # Plot significance as a function of mu and e selection efficiencies
    do_plot_sigma_v_eff = True
    if do_plot_sigma_v_eff:
        plot_sigma_v_eff(enu_m, enu_e, enu_s, fcov, pot=6.6e20)

    # Print interesting significances
    if True:
        em = enu_m
        ee = enu_e
        es = enu_s
        #print('6.6e20 syst chi2: %f' %
        #       significance(em, ee, es, fcov,
        #       pot=66e19, syst=True, mc=False, eff_m=0.3, eff_e=0.3))

        #print('6.6e20 stat chi2: %f' %
        #       significance(em, ee, es, fcov,
        #       pot=66e19, syst=False, mc=False, eff_m=0.3, eff_e=0.3))

        #print('5.0e19 syst chi2: %f' %
        #       significance(em, ee, es, fcov,
        #       pot=5e19, syst=True, mc=False, eff_m=0.3, eff_e=0.3))

        #print('5.0e19 stat chi2: %f' %
        #       significance(em, ee, es, fcov,
        #       pot=5e19, syst=False, mc=False, eff_m=0.3, eff_e=0.3))

        print('6.6e20 stat mc: %f' %
               significance(em, ee, es, fcov,
               pot=66e19, syst=False, mc=False, eff_m=0.3, eff_e=0.3))

        print('6.6e20 syst mc: %f' %
               significance(em, ee, es, fcov,
               pot=66e19, syst=True, mc=False, eff_m=0.3, eff_e=0.3))

        #print('5.0e19 syst mc: %f' %
        #       significance(em, ee, es, fcov,
        #       pot=5e19, syst=True, mc=True, eff_m=0.3, eff_e=0.3))

        #print('5.0e19 stat mc: %f' %
        #       significance(em, ee, es, fcov,
        #       pot=5e19, syst=False, mc=True, eff_m=0.3, eff_e=0.3))

    # Draw the total correlation matrix
    draw_cor = True
    if draw_cor:
        plot_cor(cov)

    do_constrain_nue = True
    if do_constrain_nue:
        constrain_nue(enu_mx, enu_m, enu_ex, enu_e, enu_ss, fcov,
                      pot=6.6e20, eff_e=0.3)

    plot_bin_errors(enu_mx, enu_m, enu_ex, enu_e, enu_ss, fcov,
                  pot=6.6e20, eff_e=0.3)

    # Draw the signal and background spectra
    draw_spectra = True
    if draw_spectra:
        #pot=1e21
        #e_eff = 0.3
        ##mu_effs = np.linspace(0.001, 0.2, 10) # [0.005, 0.01, 0.05, 0.1, 0.5]
        ##mu_effs = np.linspace(0.1, 1.0, 10)
        #mu_effs = np.logspace(-3, -0.5, 10) # [0.005, 0.01, 0.05, 0.1, 0.5]
        #argh = []
        #for mu_eff in mu_effs:
        #    #p, s = significance(enu_m, enu_e, enu_s, fcov, syst=False, pot=pot, eff_m=mu_eff, eff_e=e_eff)
        #    #print('chi2 stat: p=%1.5f, sigma=%1.5f' % (p, s))
        #    #p, s = significance(enu_m, enu_e, enu_s, fcov, pot=pot, eff_m=mu_eff, eff_e=e_eff)
        #    #print('chi2 syst: p=%1.5f, sigma=%1.5f' % (p, s))

        #    p, s = significance(enu_m, enu_e, enu_s, fcov, pot=pot, eff_m=mu_eff, eff_e=e_eff, mc=False)
        #    print('mc syst: %f  p=%1.5f, sigma=%1.5f' % (mu_eff, p, s))
        #    #p, s = significance(enu_m, enu_e, enu_s, fcov, syst=False, pot=6.6e20, mc=True, eff_m=mu_eff, eff_e=e_eff)
        #    #print('mc stat:   p=%1.5f, sigma=%1.5f' % (p, s))

        #    argh.append(s)

        #    #plot_spectra(enu_mx, enu_m, enu_ex, enu_e, enu_ss, fcov, pot=pot, eff_m=mu_eff, eff_e=e_eff)
        plot_spectra(enu_mx, enu_m, enu_ex, enu_e, enu_ss, fcov, pot=6.6e20)

        #plt.semilogx(mu_effs, argh)
        #plt.xlabel('$1\mu1p$ efficiency')
        #plt.ylabel('Significance ($\sigma$), $6.6\\times10^{20}$ POT')
        #plt.show()

    # Solo systematics
    solo_systs = True
    if solo_systs:
        wnames = [
            'expskin_FluxUnisim',
            'genie_AGKYpT_Genie',
            'genie_AGKYxF_Genie',
            'genie_DISAth_Genie',
            'genie_DISBth_Genie',
            'genie_DISCv1u_Genie',
            'genie_DISCv2u_Genie',
            'genie_FermiGasModelKf_Genie',
            'genie_FermiGasModelSf_Genie',
            'genie_FormZone_Genie',
            'genie_IntraNukeNabs_Genie',
            'genie_IntraNukeNcex_Genie',
            'genie_IntraNukeNel_Genie',
            'genie_IntraNukeNinel_Genie',
            'genie_IntraNukeNmfp_Genie',
            'genie_IntraNukeNpi_Genie',
            'genie_IntraNukePIabs_Genie',
            'genie_IntraNukePIcex_Genie',
            'genie_IntraNukePIel_Genie',
            'genie_IntraNukePIinel_Genie',
            'genie_IntraNukePImfp_Genie',
            'genie_IntraNukePIpi_Genie',
            'genie_NC_Genie',
            'genie_NonResRvbarp1pi_Genie',
            'genie_NonResRvbarp2pi_Genie',
            'genie_NonResRvp1pi_Genie',
            'genie_NonResRvp2pi_Genie',
            'genie_ResDecayEta_Genie',
            'genie_ResDecayGamma_Genie',
            'genie_ResDecayTheta_Genie',
            'genie_ccresAxial_Genie',
            'genie_ccresVector_Genie',
            'genie_cohMA_Genie',
            'genie_cohR0_Genie',
            'genie_ncelAxial_Genie',
            'genie_ncelEta_Genie',
            'genie_ncresAxial_Genie',
            'genie_ncresVector_Genie',
            'genie_qema_Genie',
            'genie_qevec_Genie',
            'horncurrent_FluxUnisim',
            'kminus_PrimaryHadronNormalization',
            'kplus_PrimaryHadronFeynmanScaling',
            'kzero_PrimaryHadronSanfordWang',
            'nucleoninexsec_FluxUnisim',
            'nucleonqexsec_FluxUnisim',
            'nucleontotxsec_FluxUnisim',
            'piminus_PrimaryHadronSWCentralSplineVariation',
            'pioninexsec_FluxUnisim',
            'pionqexsec_FluxUnisim',
            'piontotxsec_FluxUnisim',
            'piplus_PrimaryHadronSWCentralSplineVariation',
            'xsr_scc_Fv3_XSecRatio',
            'xsr_scc_Fa3_XSecRatio',
            'genie',
            'flux',
            'xsr',
            'all',
        ]

        s = significance(enu_m, enu_e, enu_s, fcov, pot=6.6e20, syst=False)
        print('stat %1.5f' % (s))

        for w in reversed(wnames):
            enu_mx, enu_m, enu_ex, enu_e, enu_s, enu_ss, cov = \
                load(sys.argv[1] + ('/cov_%s.root' % w), sys.argv[2] + ('/cov_%s.root' % w))  

            eb = np.hstack((enu_m, enu_e))
            fcov = cov / np.einsum('i,j', eb, eb)

            s = significance(enu_m, enu_e, enu_s, fcov, pot=6.6e20)
            print('%s %1.5f' % (w, s))

