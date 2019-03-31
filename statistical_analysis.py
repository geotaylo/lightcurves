"""

Run a statistical anaylsis of fitted vs true parameters from model lightcurves of type 1a SN, fitted with sncosmo.
Plots distributions of residuals and errors.
G. Taylor, 2018

- Some versions in here might be a bit odd, we'll work on fixing it.

"""

import csv
from prettytable import PrettyTable
import os
import matplotlib.pyplot as plt
import pylab as py
import numpy as np
from scipy.optimize import curve_fit
import copy
from matplotlib.ticker import MaxNLocator
import heapq
from scipy import stats


def analyse(folder, set, fails=[], wipe_fails=False):
    """
        Calculates residuals (true - fitted value) of SN parameters for one filter set.

        :param folder: path to specific filter set data (e.g. '~\kst\')
        :param set: name of filter set, for labelling purposes (e.g. 'kst')
        :param fails: associated failed SN from another filter set.  For example, if kst has failed to fit a lightcurve to
                        SN10, the combined sm+kst fit will not be able to use kst's fitted t0 value as an initial guess,
                        and the fit will be affected.  SN10 should be counted as a 'fail' for sm+kst, as the fit has
                        underperformed.
        :param wipe_fails: if True, the 'failed' SN will not be included in the output of residuals, showing the best
                            possible performance of the code.  If False, the entire SN set is included, regardless of fit
                            success.
        :return: error set (list of 'failed' SN) and diffs, a table of absolute residuals.
        """

    # Naming convention - 'no_fails' indicates fails have been removed from statistical set.
    if wipe_fails:
        set = set + '_nofails'

    # Check 'failed' error SN and flag (these will all be handled better at some point):
    error_file = folder + 'error_sn.txt'
    snr_error_file = folder + 'snr_error_sn.txt'
    error_set = []
    # get error sn from error file
    f = open(error_file, "r")
    lines = f.readlines()
    for x in lines:
        y = x.translate(None, "SN :")
        error_set.append(int(y))
    f.close()
    # get error sn from snr_error file
    f = open(snr_error_file, "r")
    lines = f.readlines()
    for x in lines:
        y = x.translate(None, "SN :")
        error_set.append(int(y))
    f.close()

    # Import and strip true and fitted params
    fp = folder + 'fitted_parameters.txt'
    tp = folder + 'true_parameters.txt'

    fitted_c = []
    fitted_mwebv = []
    fitted_t0 = []
    fitted_x0 = []
    fitted_x1 = []
    fitted_z = []

    # Handling fitted SN
    with open(fp, 'rb') as f:
        reader = csv.reader(f, delimiter=' ')
        for row in reader:
            fitted_c.append(float(row[1].replace('c:', '')))
            fitted_mwebv.append(float(row[2].replace('mwebv:', '')))
            fitted_t0.append(float(row[3].replace('t0:', '')))
            fitted_x0.append(float(row[4].replace('x0:', '')))
            fitted_x1.append(float(row[5].replace('x1:', '')))
            fitted_z.append(float(row[6].replace('z:', '')))

    # Redundancy check for missed fails (shouldn't occur, but just in case
    for i in range(len(fitted_c)):
        if fitted_c[i] == 0 and fitted_t0[i] == 0 and fitted_x0[i] == 0 and fitted_x1[i] == 0:
            if i + 1 not in error_set:
                error_set.append(i + 1)
                print "Found a skipped error!"
    error_set.sort()

    # Handling true SN
    true_c = []
    true_mwebv = []
    true_t0 = []
    true_x0 = []
    true_x1 = []
    true_z = []

    with open(tp, 'rb') as file:
        reader = csv.reader(file, delimiter=' ')
        for row in reader:
            true_c.append(float(row[1].replace('c:', '')))
            true_mwebv.append(float(row[2].replace('mwebv:', '')))
            true_t0.append(float(row[3].replace('t0:', '')))
            true_x0.append(float(row[4].replace('x0:', '')))
            true_x1.append(float(row[5].replace('x1:', '')))
            true_z.append(float(row[6].replace('z:', '')))

    # Handling uncertainties
    # Import and strip true and fitted params
    up = folder + 'fitted_errors.txt'

    uncert_c = []
    uncert_t0 = []
    uncert_x0 = []
    uncert_x1 = []

    with open(up, 'rb') as f:
        reader = csv.reader(f, delimiter=' ')
        for row in reader:
            # these are actually errors, I should be less lazy
            uncert_c.append(float(row[1].replace('c:', '')))
            uncert_t0.append(float(row[2].replace('t0:', '')))
            uncert_x0.append(float(row[3].replace('x0:', '')))
            uncert_x1.append(float(row[4].replace('x1:', '')))

    # Calculating residuals
    diff_c = []
    diff_t0 = []
    diff_x0 = []
    diff_x1 = []
    diff_z = []

    # For percentage difference (not using right now)
    # for i in range(len(true_c)):
    #     diff_c.append((true_c[i] - fitted_c[i])/true_c[i]*100)
    #     diff_t0.append((true_t0[i] - fitted_t0[i]))
    #     diff_x0.append((true_x0[i] - fitted_x0[i])/true_x0[i]*100)
    #     diff_x1.append((true_x1[i] - fitted_x1[i])/true_x1[i]*100)
    #     diff_z.append((true_z[i] - fitted_z[i])/true_z[i]*100)

    for i in range(len(true_c)):
        diff_c.append((true_c[i] - fitted_c[i]))
        diff_t0.append((true_t0[i] - fitted_t0[i]))
        diff_x0.append((true_x0[i] - fitted_x0[i]))
        diff_x1.append((true_x1[i] - fitted_x1[i]))
        diff_z.append((true_z[i] - fitted_z[i]))

    # Create sn_num array (an array of all SN in set, with 'fails' flagged.
    sn_num = range(1, len(diff_c) + 1)
    for i in error_set:
        sn_num[i - 1] = 'fit_error' + str(i)

    # Flag kepler errors (only applicable for combined seeing where kst t0 couldn't be passed)
    for i in fails:
        sn_num[i - 1] = 'kst_error' + str(i)

    total_fails = filter(lambda x: x in error_set, fails)

    # Handles failures in both kst and current fit
    for i in total_fails:
        sn_num[i - 1] = 'fit_and_kst_error' + str(i)

    # remove fails from data
    if wipe_fails:
        for i in sorted(error_set + fails, reverse=True):
            del sn_num[i - 1]
            del diff_c[i - 1]
            del diff_t0[i - 1]
            del diff_x0[i - 1]
            del diff_x1[i - 1]
            del diff_z[i - 1]
            del uncert_c[i - 1]
            del uncert_t0[i - 1]
            del uncert_x0[i - 1]
            del uncert_x1[i - 1]

        # for checking how many SN are removed with further cuts
        uncut_len = len(uncert_c)
        #apply quality cuts - fitted c outside [-0.3, 0.3], fitted x1 outside [-3, 3]
        #                   - c uncertainty > 0.3, x1 uncertainty > 3, t0 uncertainty > 5
        for i in sorted(range(len(sn_num)), reverse=True):
            if (not -0.3 <= fitted_c[i] <= 0.3) or (not -3.0 <= fitted_x1[i] <= 3.0) or \
                    (uncert_c[i] > 0.3) or (uncert_x1[i] > 3.0) or (uncert_t0[i] > 5):
                del sn_num[i]
                del diff_c[i]
                del diff_t0[i]
                del diff_x0[i]
                del diff_x1[i]
                del diff_z[i]
                del uncert_c[i]
                del uncert_t0[i]
                del uncert_x0[i]
                del uncert_x1[i]

    # Creates residuals table
    t = PrettyTable()
    t.title = set
    t.add_column('SN', sn_num)
    t.add_column('c-diff', diff_c)
    t.add_column('t0-diff', diff_t0)
    t.add_column('x0-diff', diff_x0)
    t.add_column('x1-diff', diff_x1)
    t.add_column('z-diff', diff_z)

    table_txt = t.get_string()

    # Create uncertainties table
    u = PrettyTable()
    u.title = set
    u.add_column('SN', sn_num)
    u.add_column('c-uncertainty', uncert_c)
    u.add_column('t0-uncertainty', uncert_t0)
    u.add_column('x0-uncertainty', uncert_x0)
    u.add_column('x1-uncertainty', uncert_x1)

    table_txt_2 = u.get_string()

    # Writes results
    writefolder = folder + "stats/"
    if not os.path.isdir(writefolder):
        os.makedirs(writefolder)

    if wipe_fails:
        with open(writefolder + 'uncertainty_output_nofails.txt', 'w') as file:
            file.write(table_txt_2)
    else:
        with open(writefolder + 'uncertainty_output.txt', 'w') as file:
            file.write(table_txt_2)

    if wipe_fails:
        with open(writefolder + 'residual_output_nofails.txt', 'w') as file:
            file.write(table_txt)
    else:
        with open(writefolder + 'residual_output.txt', 'w') as file:
            file.write(table_txt)

    diffs = [diff_c, diff_t0, diff_x0, diff_x1, diff_z]
    uncerts = [uncert_c, uncert_t0, uncert_x0, uncert_x1]

    return error_set, diffs, uncerts, uncut_len


    # # remove outliers (top and bottom 1% or 10 entries)
    # final_diff_c = sorted(diff_c)[10:-10]
    # final_diff_t0 = sorted(diff_t0)[10:-10]
    # final_diff_x0 = sorted(diff_x0)[10:-10]
    # final_diff_x1 = sorted(diff_x1)[10:-10]
    #
    # final_uncert_c = sorted(uncert_c)[10:-10]
    # final_uncert_t0 = sorted(uncert_t0)[10:-10]
    # final_uncert_x0 = sorted(uncert_x0)[10:-10]
    # final_uncert_x1 = sorted(uncert_x1)[10:-10]
    #
    #
    # # inspect largest elements in set
    # print 'c res: %s' %heapq.nlargest(5, final_diff_c)
    # print 't0 res: %s' %heapq.nlargest(5, final_diff_t0)
    # print 'x0 res: %s' %heapq.nlargest(5, final_diff_x0)
    # print 'x1 res: %s' %heapq.nlargest(5, final_diff_x1)
    #
    # print 'c uncty: %s' %heapq.nlargest(5, final_uncert_c)
    # print 't0 uncty: %s' % heapq.nlargest(5, final_uncert_t0)
    # print 'x0 uncty: %s' % heapq.nlargest(5, final_uncert_x0)
    # print 'x1 uncty: %s' % heapq.nlargest(5, final_uncert_x1)


def f(x, a, b, c):
    """
    Equation for a Gaussian
    """
    return a * py.exp((-(x - b) ** 2.0) / (2 * c ** 2))

#
# def plot_diffs(scopes, labels, colour, folder):
# #     """
# #     Plots histograms and fitted gaussians of residuals for a set of filters, for each parameter.
# #     """
# #
# #     # C - set to [-0.2,0.2] with bins of 0.01
# #     crange = 0.2
# #     trange = 0.2
# #     x1range = 0.2
# #     try:
# #         figa, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True)
# #
# #         for i in range(0, len(scopes)):
# #             abs_c = scopes[i][0]
# #             trimmed_c = [x for x in abs_c if x >= -crange and x <= crange]
# #             bins_c = np.arange(min(trimmed_c), max(trimmed_c) + 0.01, 0.01)
# #             data = ax1.hist(trimmed_c, bins_c, histtype='step', color=colour[i], label=labels[i]+'(%s fits)'%str(len(trimmed_c)))
# #
# #             # Generate data from bins as a set of points
# #             x = [0.5 * (data[1][t] + data[1][t + 1]) for t in xrange(len(data[1]) - 1)]
# #             y = data[0]
# #
# #             popt, pcov = curve_fit(f, x, y)
# #
# #             x_fit = py.linspace(-crange, crange, 200)
# #             y_fit = f(x_fit, *popt)
# #
# #             ax2.plot(x_fit, y_fit, lw=2, color=colour[i])
# #
# #             ax3.hist(trimmed_c, bins_c, histtype='step', color=colour[i], label=labels[i]+'(%s fits)'%str(len(trimmed_c)/2))
# #             ax3.plot(x_fit, y_fit, lw=2, color=colour[i])
# #
# #         if not os.path.isdir(folder):
# #             os.makedirs(folder)
# #
# #         figa.subplots_adjust(hspace=0)
# #         plt.setp([a.get_xticklabels() for a in figa.axes[:-1]], visible=False)
# #         ax1.set_title('Colour (c) Residuals')
# #         plt.xlabel('Residual (true - fitted value of c)')
# #         plt.ylabel('Frequency')
# #         plt.xlim(-crange, crange)
# #         ax1.legend(fontsize = 'x-small')
# #         figa.savefig(folder + 'colour.png', dpi=1000)
# #         plt.close()
# #     except RuntimeError:
# #         print('An error occured.')
# #         plt.close()
# #
# #     # t_0 - set to [-0.25,0.25] with bins of 0.02
# #     try:
# #         figa2, (ax4, ax5, ax6) = plt.subplots(3, sharex=True, sharey=True)
# #
# #         for i in range(0, len(scopes)):
# #             abs_t0 = scopes[i][1]
# #             trimmed_t0 = [x for x in abs_t0 if x >= -trange and x <= trange]
# #             bins_t0 = np.arange(min(trimmed_t0), max(trimmed_t0) + 0.02, 0.02)
# #             data = ax4.hist(trimmed_t0, bins_t0, histtype='step', color=colour[i], label=labels[i]+'(%s fits)'%str(len(trimmed_t0)))
# #
# #             # Generate data from bins as a set of points
# #             x = [0.5 * (data[1][t] + data[1][t + 1]) for t in xrange(len(data[1]) - 1)]
# #             y = data[0]
# #
# #             popt, pcov = curve_fit(f, x, y)
# #
# #             x_fit = py.linspace(-trange, trange, 200)
# #             y_fit = f(x_fit, *popt)
# #
# #             ax5.plot(x_fit, y_fit, lw=2, color=colour[i])
# #
# #             ax6.hist(trimmed_t0, bins_t0, histtype='step', color=colour[i], label=labels[i])
# #             ax6.plot(x_fit, y_fit, lw=2, color=colour[i])
# #
# #         if not os.path.isdir(folder):
# #             os.makedirs(folder)
# #
# #         figa2.subplots_adjust(hspace=0)
# #         plt.figure(dpi=1200)
# #         plt.setp([a.get_xticklabels() for a in figa2.axes[:-1]], visible=False)
# #         ax4.set_title(r'Explosion Time ($t_0$) Residuals')
# #         plt.xlabel(r'Residual (true - fitted value of $t_0$)')
# #         plt.ylabel('Frequency')
# #         plt.xlim(-trange, trange)
# #         ax4.legend(fontsize = 'x-small')
# #         figa2.savefig(folder + 't0.png')
# #         plt.close()
# #     except RuntimeError:
# #         print('An error occured.')
# #         plt.close()
# #
# #     # x_0 - no range because it's a bit garbage.
# #     # try:
# #     #     figa3, (ax7, ax8, ax9) = plt.subplots(3, sharex=True, sharey=True)
# #     #
# #     #     for i in range(0, len(scopes)):
# #     #         abs_x0 = scopes[i][2]
# #     #         bins_x0 = 100
# #     #         data = ax7.hist(abs_x0, bins_x0, histtype='step', color=colour[i], label=labels[i]+'(%s fits)'%str(len(abs_x0)))
# #     #
# #     #         # Generate data from bins as a set of points
# #     #         x = [0.5 * (data[1][t] + data[1][t + 1]) for t in xrange(len(data[1]) - 1)]
# #     #         y = data[0]
# #     #
# #     #         popt, pcov = curve_fit(f, x, y)
# #     #
# #     #         x_fit = py.linspace(min(abs_x0), max(abs_x0), 200)
# #     #         y_fit = f(x_fit, *popt)
# #     #
# #     #         ax8.plot(x_fit, y_fit, lw=2, color=colour[i])
# #     #
# #     #         ax9.hist(abs_x0, bins_x0, histtype='step', color=colour[i], label=labels[i])
# #     #         ax9.plot(x_fit, y_fit, lw=2, color=colour[i])
# #     #
# #     #     if not os.path.isdir(folder):
# #     #         os.makedirs(folder)
# #     #
# #     #     figa3.subplots_adjust(hspace=0)
# #     #     plt.setp([a.get_xticklabels() for a in figa3.axes[:-1]], visible=False)
# #     #     ax7.set_title(r'$x_0$ Residuals')
# #     #     plt.xlabel(r'Residual (true - fitted value of $x_0$)')
# #     #     plt.ylabel('Frequency')
# #     #     ax7.legend(fontsize = 'x-small')
# #     #     plt.xlim(-1, 1)
# #     #     figa3.savefig(folder + 'x0.png')
# #     #     plt.close()
# #     #     # plt.show()
# #     # except RuntimeError:
# #     #     print('An error occured.')
# #     #     plt.close()
# #
# #     # x_1 - set to [-0.2, 0.2] with bins of 0.01
# #     try:
# #         figa4, (ax10, ax11, ax12) = plt.subplots(3, sharex=True, sharey=True)
# #
# #         for i in range(0, len(scopes)):
# #             abs_x1 = scopes[i][3]
# #             trimmed_x1 = [x for x in abs_x1 if x >= -x1range and x <= x1range]
# #             bins_x1 = np.arange(min(trimmed_x1), max(trimmed_x1) + 0.01, 0.01)
# #             data = ax10.hist(trimmed_x1, bins_x1, histtype='step', color=colour[i], label=labels[i]+'(%s fits)'%str(len(trimmed_x1)))
# #
# #             # Generate data from bins as a set of points
# #             x = [0.5 * (data[1][t] + data[1][t + 1]) for t in xrange(len(data[1]) - 1)]
# #             y = data[0]
# #
# #             popt, pcov = curve_fit(f, x, y)
# #
# #             x_fit = py.linspace(-x1range, x1range, 200)
# #             y_fit = f(x_fit, *popt)
# #
# #             ax11.plot(x_fit, y_fit, lw=2, color=colour[i])
# #
# #             ax12.hist(trimmed_x1, bins_x1, histtype='step', color=colour[i], label=labels[i])
# #             ax12.plot(x_fit, y_fit, lw=2, color=colour[i])
# #
# #         if not os.path.isdir(folder):
# #             os.makedirs(folder)
# #
# #         figa4.subplots_adjust(hspace=0)
# #         plt.figure(dpi=1200)
# #         plt.setp([a.get_xticklabels() for a in figa4.axes[:-1]], visible=False)
# #         ax10.set_title(r' Stretch ($x_1$) Residuals')
# #         plt.xlabel(r'Residual (true - fitted value of $x_1$)')
# #         plt.ylabel('Frequency')
# #         ax10.legend(fontsize = 'x-small')
# #         plt.xlim(-x1range, x1range)
# #         figa4.savefig(folder + 'x1.png')
# #         figa4.clf()
# #         # plt.show()
# #         plt.close()
# #     except RuntimeError:
# #         print('An error occured.')
# #         plt.close()
# #
#      return

def plot_boxplot(scopes, label, colour, folder, mode):
    """
    Plots histograms and fitted gaussians of residuals for a set of filters, for each parameter.
    """

    if not os.path.isdir(folder):
        os.makedirs(folder)

    plt.style.use('seaborn-paper')

    try:
        fig, ax = plt.subplots()
        index = 0
        for i in range(0, len(scopes)):
            abs_c = scopes[i][0]
            boxprops = dict(linestyle='-', linewidth=0.7, color=colour[i], facecolor='white')
            flierprops = dict(marker='s', markerfacecolor=colour[i], markersize=3,
                              linestyle='none', alpha=0.3)
            genprops=dict(linestyle='-', linewidth=0.7, color=colour[i])

            plt.boxplot(abs_c, positions=[i], patch_artist=True,
                        boxprops=boxprops,
                        whiskerprops=genprops,
                        capprops=genprops,
                        whis=1.7239,#shows 3 sigma within whiskers
                        medianprops=genprops,
                        flierprops=flierprops)
            index += 1
        if mode == 'residual':
            ax.set_ylabel('Residual (true - fitted value of c)')
        elif mode == 'uncertainty':
            ax.set_ylabel('Uncertainty in fitted value of c')
        plt.xticks(range(0,index), label)
        ax.set_xlim(-0.5, index+0.5)
        plt.xticks(rotation=45)
        fig.savefig(folder + mode + '_colour.png', dpi=200, bbox_inches='tight')
        plt.close(fig)
    except RuntimeError:
        print('An error occured - probably couldn\'t fit a Gaussian.')
        plt.close()

    try:
        fig, ax = plt.subplots()
        index = 0
        for i in range(0, len(scopes)):
            abs_c = scopes[i][1]
            boxprops = dict(linestyle='-', linewidth=0.7, color=colour[i], facecolor='white')
            flierprops = dict(marker='s', markerfacecolor=colour[i], markersize=3,
                              linestyle='none', alpha=0.3)
            genprops=dict(linestyle='-', linewidth=0.7, color=colour[i])

            a = plt.boxplot(abs_c, positions=[i], patch_artist=True,
                        boxprops=boxprops,
                        whiskerprops=genprops,
                        capprops=genprops,
                        whis=1.7239,#shows 3 sigma within whiskers
                        medianprops=genprops,
                        flierprops=flierprops)
            index += 1
        if mode == 'residual':
            ax.set_ylabel('Residual (true - fitted value of t0)')
        elif mode == 'uncertainty':
            ax.set_ylabel('Uncertainty in fitted value of t0')
        else:
            raise ValueError
        plt.xticks(range(0,index), label)
        ax.set_xlim(-0.5, index+0.5)
        plt.xticks(rotation=45)
        fig.savefig(folder + mode + '_t0.png', dpi=200, bbox_inches='tight')
        plt.close(fig)
    except RuntimeError:
        print('An error occured - probably couldn\'t fit a Gaussian.')
        plt.close()

    try:
        fig, ax = plt.subplots()
        index = 0
        for i in range(0, len(scopes)):
            abs_c = scopes[i][3]
            boxprops = dict(linestyle='-', linewidth=0.7, color=colour[i], facecolor='white')
            flierprops = dict(marker='s', markerfacecolor=colour[i], markersize=3,
                              linestyle='none', alpha=0.3)
            genprops=dict(linestyle='-', linewidth=0.7, color=colour[i])

            plt.boxplot(abs_c, positions=[i], patch_artist=True,
                        boxprops=boxprops,
                        whiskerprops=genprops,
                        capprops=genprops,
                        whis=1.7239,#shows 3 sigma within whiskers
                        medianprops=genprops,
                        flierprops=flierprops)
            index += 1
        if mode == 'residual':
            ax.set_ylabel('Residual (true - fitted value of x1)')
        elif mode == 'uncertainty':
            ax.set_ylabel('Uncertainty in fitted value of x1')
        plt.xticks(range(0,index), label)
        ax.set_xlim(-0.5, index+0.5)
        plt.xticks(rotation=45)
        fig.savefig(folder + mode +'_stretch.png', dpi=200, bbox_inches='tight')
        plt.close(fig)
    except RuntimeError:
        print('An error occured - probably couldn\'t fit a Gaussian.')
        plt.close()
    return

def plot_trimmed(scopes, label, colour, folder, mode):

    if not os.path.isdir(folder):
        os.makedirs(folder)

    plt.style.use('seaborn-paper')

    try:
        fig, ax = plt.subplots()
        fig_2, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True)

        index = 0
        for i in range(0, len(scopes)):
            data = scopes[i][0]

            # Apply 3 sigma data cut
            upper_quartile = np.percentile(data, 75)
            lower_quartile = np.percentile(data, 25)
            iqr = upper_quartile - lower_quartile
            trimmed_data = [x for x in data if x <= upper_quartile + 1.7239 * iqr and x >= lower_quartile - 1.7239 * iqr]

            # Create boxplot
            boxprops = dict(linestyle='-', linewidth=0.7, color=colour[i], facecolor='white')
            flierprops = dict(marker='s', markerfacecolor=colour[i], markersize=3,
                              linestyle='none', alpha=0.3)
            genprops = dict(linestyle='-', linewidth=0.7, color=colour[i])

            ax.boxplot(trimmed_data, positions=[i], patch_artist=True,
                        boxprops=boxprops,
                        whiskerprops=genprops,
                        capprops=genprops,
                        whis=1.7239,  # shows 3 sigma within whiskers
                        medianprops=genprops,
                        flierprops=flierprops)
            index += 1

            if mode == 'residual':
                # Create Gaussians
                bins = np.arange(min(trimmed_data), max(trimmed_data) + 0.01, 0.01)
                hist_data = ax1.hist(trimmed_data, bins, histtype='step', density=True, color=colour[i],
                                label=label[i] + '(%s fits)' % str(len(trimmed_data)))
                # Generate data from bins as a set of points
                x = [0.5 * (hist_data[1][t] + hist_data[1][t + 1]) for t in xrange(len(hist_data[1]) - 1)]
                y = hist_data[0]

                popt, pcov = curve_fit(f, x, y)

                x_fit = py.linspace(min(trimmed_data), max(trimmed_data), 200)
                y_fit = f(x_fit, *popt)

                ax2.plot(x_fit, y_fit, lw=1, color=colour[i],
                         label=label[i] + '\n $\mu = %s$' % round(popt[1], 5) + '\n $\sigma = %s$' % round(popt[2], 5))

                ax3.hist(trimmed_data, bins, histtype='step', density=True, color=colour[i],
                         label=label[i] + '(%s fits)' % str(len(trimmed_data) / 2))
                ax3.plot(x_fit, y_fit, lw=1, color=colour[i])

        if mode == 'residual':
            ax.set_ylabel('Residual (true - fitted value of c)')
            ax3.set_xlabel('Residual (true - fitted value of c)')
            ax2.set_ylabel('PDF')
            leg = ax2.legend(fontsize='large', loc='center left', bbox_to_anchor=(1, 0.5), ncol=1,
                             borderaxespad=0, frameon=False, labelspacing=1)
            for line in leg.get_lines():
                line.set_linewidth(3.0)
            fig_2.savefig(folder + mode + '_colour_normed.png', dpi=200, bbox_extra_artists=(leg,), bbox_inches='tight')
        elif mode == 'uncertainty':
            ax.set_ylabel('Uncertainty in fitted value of c')

        ax.set_xticks(range(0, index), label)
        ax.set_xlim(-0.5, index + 0.5)
        fig.savefig(folder + mode + '_colour.png', dpi=200, bbox_inches='tight')

        plt.close(fig)
        plt.close(fig_2)
    except RuntimeError:
        print('An error occured - probably couldn\'t fit a Gaussian.')
        plt.close()

    try:
        fig, ax = plt.subplots()
        fig_2, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True)

        index = 0
        for i in range(0, len(scopes)):
            data = scopes[i][1]

            # Apply 3 sigma data cut
            upper_quartile = np.percentile(data, 75)
            lower_quartile = np.percentile(data, 25)
            iqr = upper_quartile - lower_quartile
            trimmed_data = [x for x in data if x <= upper_quartile + 1.7239 * iqr and x >= lower_quartile - 1.7239 * iqr]

            # Create boxplot
            boxprops = dict(linestyle='-', linewidth=0.7, color=colour[i], facecolor='white')
            flierprops = dict(marker='s', markerfacecolor=colour[i], markersize=3,
                              linestyle='none', alpha=0.3)
            genprops = dict(linestyle='-', linewidth=0.7, color=colour[i])

            ax.boxplot(trimmed_data, positions=[i], patch_artist=True,
                        boxprops=boxprops,
                        whiskerprops=genprops,
                        capprops=genprops,
                        whis=1.7239,  # shows 3 sigma within whiskers
                        medianprops=genprops,
                        flierprops=flierprops)
            index += 1

            if mode == 'residual':
                # Create Gaussians
                bins = np.arange(min(trimmed_data), max(trimmed_data) + 0.05, 0.05)
                hist_data = ax1.hist(trimmed_data, bins, histtype='step', density=True, color=colour[i],
                                label=label[i] + '(%s fits)' % str(len(trimmed_data)))
                # Generate data from bins as a set of points
                x = [0.5 * (hist_data[1][t] + hist_data[1][t + 1]) for t in xrange(len(hist_data[1]) - 1)]
                y = hist_data[0]

                popt, pcov = curve_fit(f, x, y)

                x_fit = py.linspace(min(trimmed_data), max(trimmed_data), 200)
                y_fit = f(x_fit, *popt)

                ax2.plot(x_fit, y_fit, lw=1, color=colour[i],
                         label=label[i] + '\n $\mu = %s$' % round(popt[1], 5) + '\n $\sigma = %s$' % round(popt[2], 5))

                ax3.hist(trimmed_data, bins, histtype='step', density=True, color=colour[i],
                         label=label[i] + '(%s fits)' % str(len(trimmed_data) / 2))
                ax3.plot(x_fit, y_fit, lw=1, color=colour[i])

        if mode == 'residual':
            ax.set_ylabel('Residual (true - fitted value of t0)')
            ax3.set_xlabel('Residual (true - fitted value of t0)')
            ax2.set_ylabel('PDF')
            leg = ax2.legend(fontsize='large', loc='center left', bbox_to_anchor=(1, 0.5), ncol=1,
                             borderaxespad=0, frameon=False, labelspacing=1)
            for line in leg.get_lines():
                line.set_linewidth(3.0)
            fig_2.savefig(folder + mode + '_t0_normed.png', dpi=200, bbox_extra_artists=(leg,), bbox_inches='tight')
        elif mode == 'uncertainty':
            ax.set_ylabel('Uncertainty in fitted value of t0')

        ax.set_xticks(range(0, index), label)
        ax.set_xlim(-0.5, index + 0.5)
        fig.savefig(folder + mode + '_t0.png', dpi=200, bbox_inches='tight')

        plt.close(fig)
        plt.close(fig_2)
    except RuntimeError:
        print('An error occured - probably couldn\'t fit a Gaussian.')
        plt.close()

    try:
        fig, ax = plt.subplots()
        fig_2, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True)

        index = 0
        for i in range(0, len(scopes)):
            data = scopes[i][3]

            # Apply 3 sigma data cut
            upper_quartile = np.percentile(data, 75)
            lower_quartile = np.percentile(data, 25)
            iqr = upper_quartile - lower_quartile
            trimmed_data = [x for x in data if x <= upper_quartile + 1.7239 * iqr and x >= lower_quartile - 1.7239 * iqr]

            # Create boxplot
            boxprops = dict(linestyle='-', linewidth=0.7, color=colour[i], facecolor='white')
            flierprops = dict(marker='s', markerfacecolor=colour[i], markersize=3,
                              linestyle='none', alpha=0.3)
            genprops = dict(linestyle='-', linewidth=0.7, color=colour[i])

            ax.boxplot(trimmed_data, positions=[i], patch_artist=True,
                        boxprops=boxprops,
                        whiskerprops=genprops,
                        capprops=genprops,
                        whis=1.7239,  # shows 3 sigma within whiskers
                        medianprops=genprops,
                        flierprops=flierprops)
            index += 1

            if mode == 'residual':
                # Create Gaussians
                bins = np.arange(min(trimmed_data), max(trimmed_data) + 0.01, 0.01)
                hist_data = ax1.hist(trimmed_data, bins, histtype='step', density=True, color=colour[i],
                                label=label[i] + '(%s fits)' % str(len(trimmed_data)))
                # Generate data from bins as a set of points
                x = [0.5 * (hist_data[1][t] + hist_data[1][t + 1]) for t in xrange(len(hist_data[1]) - 1)]
                y = hist_data[0]

                popt, pcov = curve_fit(f, x, y)

                x_fit = py.linspace(min(trimmed_data), max(trimmed_data), 200)
                y_fit = f(x_fit, *popt)

                ax2.plot(x_fit, y_fit, lw=1, color=colour[i],
                         label=label[i] + '\n $\mu = %s$' % round(popt[1], 5) + '\n $\sigma = %s$' % round(popt[2], 5))

                ax3.hist(trimmed_data, bins, histtype='step', density=True, color=colour[i],
                         label=label[i] + '(%s fits)' % str(len(trimmed_data) / 2))
                ax3.plot(x_fit, y_fit, lw=1, color=colour[i])

        if mode == 'residual':
            ax.set_ylabel('Residual (true - fitted value of x1)')
            ax3.set_xlabel('Residual (true - fitted value of x1)')
            ax2.set_ylabel('PDF')
            leg = ax2.legend(fontsize='large', loc='center left', bbox_to_anchor=(1, 0.5), ncol=1,
                             borderaxespad=0, frameon=False, labelspacing=1)
            for line in leg.get_lines():
                line.set_linewidth(3.0)
            fig_2.savefig(folder + mode + '_stretch_normed.png', dpi=200, bbox_extra_artists=(leg,), bbox_inches='tight')
        elif mode == 'uncertainty':
            ax.set_ylabel('Uncertainty in fitted value of x1')

        ax.set_xticks(range(0, index), label)
        ax.set_xlim(-0.5, index + 0.5)
        fig.savefig(folder + mode + '_stretch.png', dpi=200, bbox_inches='tight')

        plt.close(fig)
        plt.close(fig_2)
    except RuntimeError:
        print('An error occured - probably couldn\'t fit a Gaussian.')
        plt.close()

    return



def plot_diffs_norm(scopes, labels, colour, folder):
    """
    Plots histograms and fitted gaussians of residuals for a set of filters, for each parameter.
    """

    crange = 100
    trange = 100
    x1range = 100

    plt.style.use('seaborn-paper')
    # C - set to [-0.2,0.2] with bins of 0.01
    try:
        figa, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True)
        for i in range(0, len(scopes)):
            abs_c = scopes[i][0]
            trimmed_c = [x for x in abs_c if x >= -crange and x <= crange]
            bins_c = np.arange(min(trimmed_c), max(trimmed_c) + 0.01, 0.01)
            data = ax1.hist(trimmed_c, bins_c, histtype='step', density=True, color=colour[i], label=labels[i]+'(%s fits)'%str(len(trimmed_c)))

            # Generate data from bins as a set of points
            x = [0.5 * (data[1][t] + data[1][t + 1]) for t in xrange(len(data[1]) - 1)]
            y = data[0]

            popt, pcov = curve_fit(f, x, y)

            x_fit = py.linspace(min(trimmed_c), max(trimmed_c), 200)
            y_fit = f(x_fit, *popt)

            ax2.plot(x_fit, y_fit, lw=2, color=colour[i], label=labels[i]+'\n $\mu = %s$'%round(popt[1], 5)+'\n $\sigma = %s$'%round(popt[2], 5))

            ax3.hist(trimmed_c, bins_c, histtype='step', density=True, color=colour[i], label=labels[i]+'(%s fits)'%str(len(trimmed_c)/2))
            ax3.plot(x_fit, y_fit, lw=2, color=colour[i])

        if not os.path.isdir(folder):
            os.makedirs(folder)

        #figa.subplots_adjust(hspace=0)
        plt.figure(dpi=1200)
        #ax2.set_xlim(-1, 1)
        #plt.setp([a.get_xticklabels() for a in figa.axes[:-1]], visible=False)
        ax3.set_xlabel('Residual (true - fitted value of c)')
        ax2.set_ylabel('PDF')
        leg = ax2.legend(fontsize='large', loc='center left', bbox_to_anchor= (1, 0.5), ncol=1,
            borderaxespad=0, frameon=False, labelspacing=1)
        for line in leg.get_lines():
            line.set_linewidth(3.0)
        figa.savefig(folder + 'colour_normed.png', dpi=200, bbox_extra_artists=(leg,), bbox_inches='tight')
        plt.close(figa)
    except RuntimeError:
        print('An error occured - probably couldn\'t fit a Gaussian.')
        plt.close()

    # t_0 - set to [-0.25,0.25] with bins of 0.0.02
    try:
        figa2, (ax4, ax5, ax6) = plt.subplots(3, sharex=True, sharey=True)

        for i in range(0, len(scopes)):
            abs_t0 = scopes[i][1]
            trimmed_t0 = [x for x in abs_t0 if x >= -trange and x <= trange]
            bins_t0 = np.arange(min(trimmed_t0), max(trimmed_t0) + 0.02, 0.02)
            data = ax4.hist(trimmed_t0, bins_t0, histtype='step', density=True, color=colour[i], label=labels[i]+'(%s fits)'%str(len(trimmed_t0)))

            # Generate data from bins as a set of points
            x = [0.5 * (data[1][t] + data[1][t + 1]) for t in xrange(len(data[1]) - 1)]
            y = data[0]

            popt, pcov = curve_fit(f, x, y)

            x_fit = py.linspace(-trange, trange, 200)
            y_fit = f(x_fit, *popt)

            ax5.plot(x_fit, y_fit, lw=2, color=colour[i], label=labels[i]+'\n $\mu = %s$'%round(popt[1], 5)+'\n $\sigma = %s$'%round(popt[2], 5))

            ax6.hist(trimmed_t0, bins_t0, histtype='step', density=True, color=colour[i], label=labels[i])
            ax6.plot(x_fit, y_fit, lw=2, color=colour[i])

        if not os.path.isdir(folder):
            os.makedirs(folder)

        figa2.subplots_adjust(hspace=0)
        plt.figure(dpi=1200)
        plt.setp([a.get_xticklabels() for a in figa2.axes[:-1]], visible=False)
        ax6.set_xlabel(r'Residual (true - fitted value of $t_0$)')
        ax5.set_ylabel('PDF')
        ax5.set_xlim(-0.1, 0.1)
        leg = ax5.legend(fontsize='large', loc='center left', bbox_to_anchor=(1, 0.5), ncol=1,
                         borderaxespad=0, frameon=False, labelspacing=1)
        for line in leg.get_lines():
            line.set_linewidth(4.0)
        figa2.savefig(folder + 't0_normed.png', dpi=200, bbox_extra_artists=(leg,), bbox_inches='tight')
        plt.close(figa2)
    except RuntimeError:
        print('An error occured.')
        plt.close()

    # x_1 - set to [-0.2, 0.2] with bins of 0.01
    try:
        figa4, (ax10, ax11, ax12) = plt.subplots(3, sharex=True, sharey=True)

        for i in range(0, len(scopes)):
            abs_x1 = scopes[i][3]
            trimmed_x1 = [x for x in abs_x1 if x >= -x1range and x <= x1range]
            bins_x1 = np.arange(min(trimmed_x1), max(trimmed_x1) + 0.01, 0.01)
            data = ax10.hist(trimmed_x1, bins_x1, histtype='step', density=True, color=colour[i], label=labels[i]+'(%s fits)'%str(len(trimmed_x1)))

            # Generate data from bins as a set of points
            x = [0.5 * (data[1][t] + data[1][t + 1]) for t in xrange(len(data[1]) - 1)]
            y = data[0]

            popt, pcov = curve_fit(f, x, y)

            x_fit = py.linspace(-x1range, x1range, 200)
            y_fit = f(x_fit, *popt)

            ax11.plot(x_fit, y_fit, lw=2, color=colour[i], label=labels[i]+'\n $\mu = %s$'%round(popt[1], 5)+'\n $\sigma = %s$'%round(popt[2], 5))

            ax12.hist(trimmed_x1, bins_x1, histtype='step', density=True, color=colour[i], label=labels[i])
            ax12.plot(x_fit, y_fit, lw=2, color=colour[i])

        if not os.path.isdir(folder):
            os.makedirs(folder)

        figa4.subplots_adjust(hspace=0)
        plt.figure(dpi=1200)
        plt.setp([a.get_xticklabels() for a in figa4.axes[:-1]], visible=False)
        ax12.set_xlabel(r'Residual (true - fitted value of $x_1$)')
        ax11.set_ylabel('PDF')
        ax11.set_xlim(-0.1, 0.1)
        leg = ax11.legend(fontsize='large', loc='center left', bbox_to_anchor=(1, 0.5), ncol=1,
                         borderaxespad=0, frameon=False, labelspacing=1)
        for line in leg.get_lines():
            line.set_linewidth(4.0)
        figa4.savefig(folder + 'x1_normed.png', dpi=200, bbox_extra_artists=(leg,), bbox_inches='tight')
        plt.close(figa4)
    except RuntimeError:
        print('An error occured.')
        plt.close()

    return

# TWO OLD METHODS THAT USE BIMODAL DISTRIBUTIONS - CLEAR THEM OUT WHEN YOU CAN
# def plot_errors_old(scopes, labels, colour, folder):
#     # """
#     # Plots histograms and fitted gaussians of residuals for a set of filters, for each parameter.
#     # """
#     # crange = 0.5
#     # trange = 0.5
#     # x1range = 0.5
#     #
#     # # C - set to [-0.2,0.2] with bins of 0.01
#     # try:
#     #     figa, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True)
#     #     for i in range(0, len(scopes)):
#     #         abs_c = copy.deepcopy(scopes[i][0])
#     #         c_flips = [-x for x in abs_c]
#     #         abs_c.extend(c_flips)
#     #         trimmed_c = [x for x in abs_c if x >= -crange and x <= crange]
#     #         bins_c = np.arange(min(trimmed_c), max(trimmed_c) + 0.01, 0.01)
#     #         data = ax1.hist(trimmed_c, bins_c, histtype='step', color=colour[i], label=labels[i]+'(%s fits)'%str(len(trimmed_c)/2))
#     #
#     #         # Generate data from bins as a set of points
#     #         x = [0.5 * (data[1][t] + data[1][t + 1]) for t in xrange(len(data[1]) - 1)]
#     #         y = data[0]
#     #
#     #         popt, pcov = curve_fit(f, x, y)
#     #
#     #         x_fit = py.linspace(-crange, crange, 200)#(min(abs_c), max(abs_c), 200)
#     #         y_fit = f(x_fit, *popt)
#     #
#     #         ax2.plot(x_fit, y_fit, lw=2, color=colour[i])
#     #
#     #         ax3.hist(trimmed_c, bins_c, histtype='step', color=colour[i], label=labels[i])
#     #         ax3.plot(x_fit, y_fit, lw=2, color=colour[i])
#     #
#     #     if not os.path.isdir(folder):
#     #         os.makedirs(folder)
#     #
#     #     figa.subplots_adjust(hspace=0)
#     #     plt.figure(dpi=1200)
#     #     plt.setp([a.get_xticklabels() for a in figa.axes[:-1]], visible=False)
#     #     ax1.set_title(r'Colour (c) Errors')
#     #     plt.xlabel('c error')
#     #     plt.ylabel('Frequency')
#     #     ax1.legend(fontsize = 'x-small')
#     #     figa.savefig(folder + 'colour.png')
#     #     plt.close()
#     # except RuntimeError as e:
#     #     print(e.message)
#     #     plt.close()
#     #
#     # # t_0 - set to [-1, 1] with bins of 0.1
#     # try:
#     #     figa2, (ax4, ax5, ax6) = plt.subplots(3, sharex=True, sharey=True)
#     #
#     #     for i in range(0, len(scopes)):
#     #         abs_t0 = copy.deepcopy(scopes[i][1])
#     #         t0_flips = [-x for x in abs_t0]
#     #         abs_t0.extend(t0_flips)
#     #         trimmed_t0 = [x for x in abs_t0 if x >= -trange and x <= trange]
#     #         bins_t0 = np.arange(min(trimmed_t0), max(trimmed_t0) + 0.02, 0.02)
#     #         data = ax4.hist(trimmed_t0, bins_t0, histtype='step', color=colour[i], label=labels[i]+'(%s fits)'%str(len(trimmed_t0)/2))
#     #
#     #         # Generate data from bins as a set of points
#     #         x = [0.5 * (data[1][t] + data[1][t + 1]) for t in xrange(len(data[1]) - 1)]
#     #         y = data[0]
#     #
#     #         popt, pcov = curve_fit(f, x, y)
#     #
#     #         x_fit = py.linspace(-trange, trange, 200)
#     #         y_fit = f(x_fit, *popt)
#     #
#     #         ax5.plot(x_fit, y_fit, lw=2, color=colour[i])
#     #
#     #         ax6.hist(trimmed_t0, bins_t0, histtype='step', color=colour[i], label=labels[i])
#     #         ax6.plot(x_fit, y_fit, lw=2, color=colour[i])
#     #
#     #     if not os.path.isdir(folder):
#     #         os.makedirs(folder)
#     #
#     #     figa2.subplots_adjust(hspace=0)
#     #     plt.figure(dpi=1200)
#     #     plt.setp([a.get_xticklabels() for a in figa2.axes[:-1]], visible=False)
#     #     ax4.set_title(r'Explosion Time ($t_0$) Errors')
#     #     plt.xlabel(r'$t_0$ error')
#     #     plt.ylabel('Frequency')
#     #     ax4.legend(fontsize = 'x-small')
#     #     figa2.savefig(folder + 't0.png')
#     #     plt.close()
#     # except RuntimeError:
#     #     print('An error occured.')
#     #     plt.close()
#     #
#     # # x_0 - unimportant.
#     # # try:
#     # #     figa3, (ax7, ax8, ax9) = plt.subplots(3, sharex=True, sharey=True)
#     # #
#     # #     for i in range(0, len(scopes)):
#     # #         abs_x0 = scopes[i][2]
#     # #         x0_flips = [-x for x in abs_x0]
#     # #         abs_x0.extend(x0_flips)
#     # #         trimmed_x0 = [x for x in abs_x0 if x >= -5 and x <= 5]
#     # #         bins_x0 = 100
#     # #         data = ax7.hist(abs_x0, bins_x0, histtype='step', color=colour[i], label=labels[i]+'(%s fits)'%str(len(abs_x0)/2))
#     # #
#     # #         # Generate data from bins as a set of points
#     # #         x = [0.5 * (data[1][t] + data[1][t + 1]) for t in xrange(len(data[1]) - 1)]
#     # #         y = data[0]
#     # #
#     # #         popt, pcov = curve_fit(f, x, y)
#     # #
#     # #         x_fit = py.linspace(min(abs_x0), max(abs_x0), 200)
#     # #         y_fit = f(x_fit, *popt)
#     # #
#     # #         ax8.plot(x_fit, y_fit, lw=2, color=colour[i])
#     # #
#     # #         ax9.hist(abs_x0, bins_x0, histtype='step', color=colour[i], label=labels[i])
#     # #         ax9.plot(x_fit, y_fit, lw=2, color=colour[i])
#     # #
#     # #     if not os.path.isdir(folder):
#     # #         os.makedirs(folder)
#     # #
#     # #
#     # #     figa3.subplots_adjust(hspace=0)
#     # #     plt.setp([a.get_xticklabels() for a in figa3.axes[:-1]], visible=False)
#     # #     ax7.set_title(r'$x_0$ Errors')
#     # #     plt.xlabel(r'$x_0$ Error')
#     # #     plt.ylabel('Frequency')
#     # #     ax7.legend(fontsize = 'x-small')
#     # #     figa3.savefig(folder + 'x0.png')
#     # #     plt.close()
#     # #     # plt.show()
#     # # except RuntimeError:
#     # #     print('An error occured.')
#     # #     plt.close()
#     #
#     #
#     #
#     # # x_1 - set to [-0.5, 0.5] with bins of 0.01 (nbins=100)
#     # try:
#     #     figa4, (ax10, ax11, ax12) = plt.subplots(3, sharex=True, sharey=True)
#     #
#     #     for i in range(0, len(scopes)):
#     #         abs_x1 = copy.deepcopy(scopes[i][3])
#     #         x1_flips = [-x for x in abs_x1]
#     #         abs_x1.extend(x1_flips)
#     #         trimmed_x1 = [x for x in abs_x1 if x >= -x1range and x <= x1range]
#     #         bins_x1 = np.arange(min(trimmed_x1), max(trimmed_x1) + 0.01, 0.01)
#     #         data = ax10.hist(trimmed_x1, bins_x1, histtype='step', color=colour[i], label=labels[i]+'(%s fits)'%str(len(trimmed_x1)/2))
#     #
#     #         # Generate data from bins as a set of points
#     #         x = [0.5 * (data[1][t] + data[1][t + 1]) for t in xrange(len(data[1]) - 1)]
#     #         y = data[0]
#     #
#     #         popt, pcov = curve_fit(f, x, y)
#     #
#     #         x_fit = py.linspace(-x1range, x1range, 200)
#     #         y_fit = f(x_fit, *popt)
#     #
#     #         ax11.plot(x_fit, y_fit, lw=2, color=colour[i])
#     #
#     #         ax12.hist(trimmed_x1, bins_x1, histtype='step', color=colour[i], label=labels[i])
#     #         ax12.plot(x_fit, y_fit, lw=2, color=colour[i])
#     #
#     #     if not os.path.isdir(folder):
#     #         os.makedirs(folder)
#     #
#     #     figa4.subplots_adjust(hspace=0)
#     #     plt.figure(dpi=1200)
#     #     plt.setp([a.get_xticklabels() for a in figa4.axes[:-1]], visible=False)
#     #     ax10.set_title(r'Fitted $x_1$ Errors')
#     #     plt.xlabel(r'Fitted $x_1$ Error')
#     #     plt.ylabel('Frequency')
#     #     ax10.legend(fontsize = 'x-small')
#     #     figa4.savefig(folder + 'x1.png')
#     #     plt.close()
#     # except RuntimeError:
#     #     print('An error occured.')
#     #     plt.close()
#     #
#     # return
#     return
#
# def plot_errors_norm_old(scopes, labels, colour, folder):
#     """
#     Plots histograms and fitted gaussians of residuals for a set of filters, for each parameter.
#     """
#
#     crange = 0.5
#     trange = 0.5
#     x1range = 0.5
#     # C - set to [-0.2,0.2] with bins of 0.01
#     try:
#         figa, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True)
#         for i in range(0, len(scopes)):
#             abs_c = copy.deepcopy(scopes[i][0])
#             c_flips = [-x for x in abs_c]
#             abs_c.extend(c_flips)
#             trimmed_c = [x for x in abs_c if x >= -crange and x <= crange]
#             bins_c = np.arange(min(trimmed_c), max(trimmed_c) + 0.01, 0.01)
#             data = ax1.hist(trimmed_c, bins_c, histtype='step', density=True, color=colour[i], label=labels[i]+'(%s fits)'%str(len(trimmed_c)/2))
#
#             # Generate data from bins as a set of points
#             x = [0.5 * (data[1][t] + data[1][t + 1]) for t in xrange(len(data[1]) - 1)]
#             y = data[0]
#
#             popt, pcov = curve_fit(f, x, y)
#
#             x_fit = py.linspace(-crange, crange, 200)#(min(abs_c), max(abs_c), 200)
#             y_fit = f(x_fit, *popt)
#
#             ax2.plot(x_fit, y_fit, lw=2, color=colour[i], label=labels[i]+'\n $\mu = %s$'%round(popt[1], 5)+'\n $\sigma = %s$'%round(popt[2], 5))
#
#             ax3.hist(trimmed_c, bins_c, histtype='step', density=True, color=colour[i], label=labels[i])
#             ax3.plot(x_fit, y_fit, lw=2, color=colour[i])
#
#         if not os.path.isdir(folder):
#             os.makedirs(folder)
#
#         figa.subplots_adjust(hspace=0)
#         plt.figure(dpi=1200)
#         plt.setp([a.get_xticklabels() for a in figa.axes[:-1]], visible=False)
#         ax3.set_xlabel('c error')
#         ax2.set_ylabel('PDF')
#         ax2.set_xlim(-0.2,0.2)
#         leg = ax2.legend(fontsize='large', loc='center left', bbox_to_anchor=(1, 0.5), ncol=1,
#                          borderaxespad=0, frameon=False, labelspacing=1)
#         for line in leg.get_lines():
#             line.set_linewidth(4.0)
#         figa.savefig(folder + 'colour_normed.png', dpi=200, bbox_extra_artists=(leg,), bbox_inches='tight')
#         plt.close(figa)
#     except RuntimeError as e:
#         print(e.message)
#         plt.close()
#
#     # t_0 - set to [-1, 1] with bins of 0.1
#     try:
#         figa2, (ax4, ax5, ax6) = plt.subplots(3, sharex=True, sharey=True)
#
#         for i in range(0, len(scopes)):
#             abs_t0 = copy.deepcopy(scopes[i][1])
#             t0_flips = [-x for x in abs_t0]
#             abs_t0.extend(t0_flips)
#             trimmed_t0 = [x for x in abs_t0 if x >= -trange and x <= trange]
#             bins_t0 = np.arange(min(trimmed_t0), max(trimmed_t0) + 0.02, 0.02)
#             data = ax4.hist(trimmed_t0, bins_t0, histtype='step', density=True, color=colour[i], label=labels[i]+'(%s fits)'%str(len(trimmed_t0)/2))
#
#             # Generate data from bins as a set of points
#             x = [0.5 * (data[1][t] + data[1][t + 1]) for t in xrange(len(data[1]) - 1)]
#             y = data[0]
#
#             popt, pcov = curve_fit(f, x, y)
#
#             x_fit = py.linspace(-trange, trange, 200)
#             y_fit = f(x_fit, *popt)
#
#             ax5.plot(x_fit, y_fit, lw=2, color=colour[i],  label=labels[i]+'\n $\mu = %s$'%round(popt[1], 5)+'\n $\sigma = %s$'%round(popt[2], 5))
#
#             ax6.hist(trimmed_t0, bins_t0, histtype='step', density=True, color=colour[i], label=labels[i])
#             ax6.plot(x_fit, y_fit, lw=2, color=colour[i])
#
#         if not os.path.isdir(folder):
#             os.makedirs(folder)
#
#         figa2.subplots_adjust(hspace=0)
#         plt.figure(dpi=1200)
#         plt.setp([a.get_xticklabels() for a in figa2.axes[:-1]], visible=False)
#         ax6.set_xlabel(r'$t_0$ error')
#         ax5.set_ylabel('PDF')
#         ax5.set_xlim(-0.2,0.2)
#         leg = ax5.legend(fontsize='large', loc='center left', bbox_to_anchor=(1, 0.5), ncol=1,
#                          borderaxespad=0, frameon=False, labelspacing=1)
#         for line in leg.get_lines():
#             line.set_linewidth(4.0)
#         figa2.savefig(folder + 't0_normed.png', dpi=200, bbox_extra_artists=(leg,), bbox_inches='tight')
#         plt.close(figa2)
#
#     except RuntimeError:
#         print('An error occured.')
#         plt.close()
#
#     # x_1 - set to [-0.5, 0.5] with bins of 0.01 (nbins=100)
#     try:
#         figa4, (ax10, ax11, ax12) = plt.subplots(3, sharex=True, sharey=True)
#
#         for i in range(0, len(scopes)):
#             abs_x1 = copy.deepcopy(scopes[i][3])
#             x1_flips = [-x for x in abs_x1]
#             abs_x1.extend(x1_flips)
#             trimmed_x1 = [x for x in abs_x1 if x >= -x1range and x <= x1range]
#             bins_x1 = np.arange(min(trimmed_x1), max(trimmed_x1) + 0.01, 0.01)
#             data = ax10.hist(trimmed_x1, bins_x1, histtype='step', density=True, color=colour[i], label=labels[i]+'(%s fits)'%str(len(trimmed_x1)/2))
#
#             # Generate data from bins as a set of points
#             x = [0.5 * (data[1][t] + data[1][t + 1]) for t in xrange(len(data[1]) - 1)]
#             y = data[0]
#
#             popt, pcov = curve_fit(f, x, y)
#
#             x_fit = py.linspace(-x1range, x1range, 200)
#             y_fit = f(x_fit, *popt)
#
#             ax11.plot(x_fit, y_fit, lw=2, color=colour[i], label=labels[i]+'\n $\mu = %s$'%round(popt[1], 5)+'\n $\sigma = %s$'%round(popt[2], 5))
#
#             ax12.hist(trimmed_x1, bins_x1, histtype='step',density=True, color=colour[i], label=labels[i])
#             ax12.plot(x_fit, y_fit, lw=2, color=colour[i])
#
#         if not os.path.isdir(folder):
#             os.makedirs(folder)
#
#         figa4.subplots_adjust(hspace=0)
#         plt.setp([a.get_xticklabels() for a in figa4.axes[:-1]], visible=False)
#         ax12.set_xlabel(r'Fitted $x_1$ Error')
#         ax11.set_ylabel('PDF')
#         ax11.set_xlim(-0.2,0.2)
#         leg = ax11.legend(fontsize='large', loc='center left', bbox_to_anchor= (1, 0.5), ncol=1,
#             borderaxespad=0, frameon=False, labelspacing=1)
#         for line in leg.get_lines():
#             line.set_linewidth(4.0)
#         figa4.savefig(folder + 'x1_normed.png', dpi=200, bbox_extra_artists=(leg,), bbox_inches='tight')
#         plt.close(figa4)
#     except RuntimeError:
#         print('An error occured.')
#         plt.close()
#
#     return


# Work in progress to fix error distributions

def plot_errors_norm(scopes, labels, colour, folder):
    """
    Plots histograms and fitted gaussians of residuals for a set of filters, for each parameter.
    """

    # Parameter ranges and resolutions - only consider errors within range
    crange = 0.5
    cres = 0.01
    trange = 0.5
    tres= 0.02
    x1range = 0.5
    x1res = 0.01

    # quick recode for label formatting

    for g in range(len(labels)):
        d = labels[g].replace('_', '\_')
        y = '$\\bf{%s}$'%d
        labels[g] = y

    # C - set to [0,crange] with bins of cres
    try:
        figa, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=False)
        # Iterate through telescope sets, plotting the errors.
        for i in range(0, len(scopes)):
            abs_c = copy.deepcopy(scopes[i][0])
            trimmed_c = [x for x in abs_c if x >= 0 and x <= crange]

            # Set bin intervals
            bins_c = np.arange(0, crange + cres, cres)
            # Plot error histograms on top plot
            data = ax1.hist(trimmed_c, bins_c, histtype='stepfilled', density=True, color=colour[i],
                            label=labels[i]+'(%s fits)'%str(len(trimmed_c)), alpha=0.4)

            # Generate data from bins as a set of points
            x = [0.5 * (data[1][t] + data[1][t + 1]) for t in xrange(len(data[1]) - 1)]
            y = data[0]
            # Fit a gaussian, plot on middle plot
            try:
                popt, pcov = curve_fit(f, x, y)
                x_fit = py.linspace(0, crange, 200)
                y_fit = f(x_fit, *popt)
            except:
                popt = [99.999, 99.999, 99.999]
                x_fit = py.linspace(0, crange, 200)
                y_fit = [1]*200
            ax2.plot(x_fit, y_fit, lw=2, color=colour[i], label=labels[i]+'\n %s fits'%str(len(trimmed_c))+'\n $\mu = %s$'%round(popt[1], 5)+'\n $\sigma = %s$'%np.abs(round(popt[2], 5)))

            # Plot data and fit overlaid on bottom plot
            ax3.hist(trimmed_c, bins_c, histtype='stepfilled', density=True, color=colour[i], label=labels[i], alpha=0.4)
            ax3.plot(x_fit, y_fit, lw=2, color=colour[i])

        if not os.path.isdir(folder):
            os.makedirs(folder)

        figa.subplots_adjust(hspace=0)
        plt.figure(dpi=1200)
        ax3.set_xlabel('c error')
        ax2.set_ylabel('PDF')
        ax2.set_xlim(0, crange)
        # Adjust tickmarks to prevent overlap and keep consistency
        ax2.set_ylim(ax1.get_ylim())
        ax3.set_ylim(ax1.get_ylim())
        ax2.set_yticks(ax1.get_yticks())
        ax3.set_yticks(ax1.get_yticks())
        nbins = len(ax1.get_yticklabels())  # added
        ax1.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune=None))  # added
        ax2.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))  # added
        ax3.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))  # added
        #
        leg = ax2.legend(fontsize='large', loc='center left', bbox_to_anchor=(1, 0.5), ncol=1,
                         borderaxespad=0, frameon=False, labelspacing=1)
        for line in leg.get_lines():
            line.set_linewidth(4.0)
        figa.savefig(folder + 'colour_normed.png', dpi=200, bbox_extra_artists=(leg,), bbox_inches='tight')
        plt.close(figa)
    except RuntimeError as e:
        print(e.message)
        plt.close()

    # # t_0
    try:
        figa2, (ax4, ax5, ax6) = plt.subplots(3, sharex=True, sharey=False)

        for i in range(0, len(scopes)):
            abs_t0 = copy.deepcopy(scopes[i][1])
            trimmed_t0 = [x for x in abs_t0 if x >= 0 and x <= trange]

            # Set bin intervals
            bins_t0 = np.arange(0, trange + tres, tres)
            # Plot error histograms on top plot
            data = ax4.hist(trimmed_t0, bins_t0, histtype='stepfilled', density=True, color=colour[i],
                            label=labels[i] + '(%s fits)' % str(len(trimmed_t0)), alpha=0.4)

            # Generate data from bins as a set of points
            x = [0.5 * (data[1][t] + data[1][t + 1]) for t in xrange(len(data[1]) - 1)]
            y = data[0]
            # Fit a gaussian, plot on middle plot
            try:
                popt, pcov = curve_fit(f, x, y)
                x_fit = py.linspace(0, trange, 200)
                y_fit = f(x_fit, *popt)
            except:
                popt = [99.999, 99.999, 99.999]
                x_fit = py.linspace(0, trange, 200)
                y_fit = [1]*200
            ax5.plot(x_fit, y_fit, lw=2, color=colour[i], label=labels[i]+'\n %s fits'%str(len(trimmed_c))+'\n $\mu = %s$'%round(popt[1], 5)+'\n $\sigma = %s$'%np.abs(round(popt[2], 5)))

            # Plot data and fit overlaid on bottom plot
            ax6.hist(trimmed_t0, bins_t0, histtype='stepfilled', density=True, color=colour[i], label=labels[i], alpha=0.4)
            ax6.plot(x_fit, y_fit, lw=2, color=colour[i])

        if not os.path.isdir(folder):
            os.makedirs(folder)

        figa2.subplots_adjust(hspace=0)
        plt.figure(dpi=1200)

        # Adjust tickmarks to prevent overlap and keep consistency
        ax5.set_ylim(ax4.get_ylim())
        ax6.set_ylim(ax4.get_ylim())
        ax5.set_yticks(ax4.get_yticks())
        ax6.set_yticks(ax4.get_yticks())
        nbins = len(ax4.get_yticklabels())  # added
        ax4.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune=None))  # added
        ax5.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))  # added
        ax6.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))  # added
        #
        ax6.set_xlabel(r'$t_0$ error')
        ax5.set_ylabel('PDF')
        ax5.set_xlim(0,trange)
        leg = ax5.legend(fontsize='large', loc='center left', bbox_to_anchor=(1, 0.5), ncol=1,
                         borderaxespad=0, frameon=False, labelspacing=1)
        for line in leg.get_lines():
            line.set_linewidth(4.0)
        figa2.savefig(folder + 't0_normed.png', dpi=200, bbox_extra_artists=(leg,), bbox_inches='tight')
        plt.close(figa2)

    except RuntimeError:
        print('An error occured.')
        plt.close()

    # # x_1
    try:
        figa4, (ax10, ax11, ax12) = plt.subplots(3, sharex=True, sharey=False)

        for i in range(0, len(scopes)):
            abs_x1 = copy.deepcopy(scopes[i][3])
            trimmed_x1 = [x for x in abs_x1 if x >= 0 and x <= x1range]

            # Set bin intervals
            bins_x1 = np.arange(0, x1range + x1res, x1res)
            # Plot error histograms on top plot
            data = ax10.hist(trimmed_x1, bins_x1, histtype='stepfilled', density=True, color=colour[i],
                            label=labels[i] + '(%s fits)' % str(len(trimmed_t0)), alpha=0.4)

            # Generate data from bins as a set of points
            x = [0.5 * (data[1][t] + data[1][t + 1]) for t in xrange(len(data[1]) - 1)]
            y = data[0]
            # Fit a gaussian, plot on middle plot
            try:
                popt, pcov = curve_fit(f, x, y)
                x_fit = py.linspace(0, x1range, 200)
                y_fit = f(x_fit, *popt)
            except:
                popt = [99.999, 99.999, 99.999]
                x_fit = py.linspace(0, x1range, 200)
                y_fit = [1] * 200
            ax11.plot(x_fit, y_fit, lw=2, color=colour[i],
                     label=labels[i]+'\n %s fits'%str(len(trimmed_c)) + '\n $\mu = %s$' % round(popt[1], 5) + '\n $\sigma = %s$' % np.abs(round(popt[2], 5)))

            # Plot data and fit overlaid on bottom plot
            ax12.hist(trimmed_x1, bins_x1, histtype='stepfilled', density=True, color=colour[i], label=labels[i],
                     alpha=0.4)
            ax12.plot(x_fit, y_fit, lw=2, color=colour[i])
        if not os.path.isdir(folder):
            os.makedirs(folder)

        figa4.subplots_adjust(hspace=0)
        # Adjust tickmarks to prevent overlap and keep consistency
        ax11.set_ylim(ax10.get_ylim())
        ax12.set_ylim(ax10.get_ylim())
        ax11.set_yticks(ax10.get_yticks())
        ax12.set_yticks(ax10.get_yticks())
        nbins = len(ax10.get_yticklabels())  # added
        ax10.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune=None))  # added
        ax11.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))  # added
        ax12.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))  # added
        #
        ax12.set_xlabel(r'Fitted $x_1$ Error')
        ax11.set_ylabel('PDF')
        ax11.set_xlim(0,x1range)
        leg = ax11.legend(fontsize='large', loc='center left', bbox_to_anchor= (1, 0.5), ncol=1,
            borderaxespad=0, frameon=False, labelspacing=1)
        for line in leg.get_lines():
            line.set_linewidth(4.0)
        figa4.savefig(folder + 'x1_normed.png', dpi=200, bbox_extra_artists=(leg,), bbox_inches='tight')
        plt.close(figa4)
    except RuntimeError:
        print('An error occured.')
        plt.close()

    return

def summary_helper(e):
    """ Returns summry satistics for a list of errors, e """
    mean =  np.mean(e)
    median = np.median(e)
    sigma = np.std(e)
    min = np.min(e)
    max = np.max(e)
    out = [mean, median, sigma, min, max]
    rounded_out = [np.round(elem, 3) for elem in out]
    return rounded_out

def uncertainty_summary(errors, lengths, scopes, colours, folder):
    """ Calculates and plots summary statistics of uncertainties - a modification of plot_errors?"""

    # Plot success rates

    # Create index array
    bar_indices = np.arange(len(errors)) * 2
    # Create array of lengths, arbitrarily using c [0] list.
    bar_lengths = [len(x[:][0]) / 10. for x in errors]
    # Text to go below scope names
    bar_strings = []#['\n %s fits \n %s cuts' % (len(x[:][0]), y) for x, y in errors, lengths]
    counter = 0
    for i in errors:
        a_string = '\n %s fits \n %s cuts' % (len(i[:][0]), lengths[counter]-len(i[:][0]))
        counter += 1
        bar_strings.append(a_string)
    bar_labels = [m + n for m, n in zip(scopes, bar_strings)]
    # Plot bars
    fig = plt.figure()
    bars_uncut = plt.bar(bar_indices, np.array(lengths) / 10., width=1.0, color=colours, alpha=0.2)
    bars = plt.bar(bar_indices, bar_lengths, width=1.0, color=colours, alpha=0.4)
    # Add counts above the bar graphs
    for rect in bars:
        height = rect.get_height()
        plt.text(rect.get_x() + rect.get_width() / 2.0, height - 5, str(height) + '%', ha='center', va='bottom')
    plt.xticks(bar_indices, bar_labels)
    plt.ylabel('Fit success rate (%)')
    fig.savefig(folder + 'success_rate.png', dpi=200, bbox_inches='tight')
    plt.close(fig)

    # For each telescope in scopes
    for i in range(len(scopes)):
        c_errors = copy.deepcopy(errors[i][0])
        t0_errors = copy.deepcopy(errors[i][1])
        x1_errors_uncut = copy.deepcopy(errors[i][3])

        x1_errors = [x for x in x1_errors_uncut if x <= 100]
        x1_cut = len([x for x in x1_errors_uncut if x > 100])

        c_summary = summary_helper(c_errors)
        t0_summary = summary_helper(t0_errors)
        x1_summary = summary_helper(x1_errors)
        try:
            fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=False, sharey=True)

            #----------COLOUR------------
            # Restrict range to 3 sigma
            trimmed_c = [x for x in c_errors if x >= 0 and x <= c_summary[2]*3.]
            # Set bin intervals
            bins_c = np.arange(0, max(trimmed_c) + 0.01, 0.01)
            # Plot c error histograms on top plot
            data = ax1.hist(trimmed_c, bins_c, histtype='stepfilled', density=True, color=colours[i],
                            label=scopes[i] + ' \n $\\bf{c}$ (%s fits within 3$\sigma$) \n $\it{mean:}$ %s \n $\it{median:}$ %s \n $\it{\sigma:}$ %s \n $\it{full\,range:}$ [%s,%s] '
                                  % (len(trimmed_c), c_summary[0], c_summary[1], c_summary[2], c_summary[3], c_summary[4]), alpha=0.4)
            ax1.set_xlabel('c error')
            ax1.set_ylabel('Frequency')
            #
            leg = ax1.legend(fontsize='large', loc='center left', bbox_to_anchor=(1, 0.5), ncol=1,
                             borderaxespad=0, frameon=False, labelspacing=1)
            for line in leg.get_lines():
                line.set_linewidth(4.0)

            # ----------EXPLOSION TIME------------
            # Restrict range to 3 sigma
            trimmed_t0 = [x for x in t0_errors if x >= 0 and x <= t0_summary[2] * 3.]
            # Set bin intervals
            bins_t0 = np.arange(0, max(trimmed_t0) + 0.02, 0.02)
            # Plot c error histograms on top plot
            data = ax2.hist(trimmed_t0, bins_t0, histtype='stepfilled', density=True, color=colours[i],
                            label='$\\bf{t_0}$ (%s fits within 3$\sigma$) \n $\it{mean:}$ %s \n $\it{median:}$ %s \n $\it{\sigma:}$ %s \n $\it{full\,range:}$ [%s,%s] '
                                  % (len(trimmed_t0), t0_summary[0], t0_summary[1], t0_summary[2], t0_summary[3],
                                     t0_summary[4]), alpha=0.4)
            ax2.set_xlabel('t0 error')
            ax2.set_ylabel('Frequency')
            #
            leg2 = ax2.legend(fontsize='large', loc='center left', bbox_to_anchor=(1, 0.5), ncol=1,
                             borderaxespad=0, frameon=False, labelspacing=1)
            for line in leg2.get_lines():
                line.set_linewidth(4.0)

            # ----------STRETCH------------
            # Restrict range to 3 sigma
            trimmed_x1 = [x for x in x1_errors if x >= 0 and x <= x1_summary[2] * 3.]
            # Set bin intervals
            bins_x1 = np.arange(0, max(trimmed_x1) + 0.02, 0.02)
            # Plot c error histograms on top plot
            data = ax3.hist(trimmed_x1, bins_x1, histtype='stepfilled', density=True, color=colours[i],
                            label='$\\bf{x_1}$ (%s fits within 3$\sigma$) \n (%s values over 100 not considered) \n $\it{mean:}$ %s \n $\it{median:}$ %s \n $\it{\sigma:}$ %s \n $\it{full\,range:}$ [%s,%s] '
                                  % (len(trimmed_x1), x1_cut, x1_summary[0], x1_summary[1], x1_summary[2], x1_summary[3],
                                     x1_summary[4]), alpha=0.4)
            ax3.set_xlabel('Absolute Error')
            ax3.set_ylabel('Frequency')
            #
            leg3 = ax3.legend(fontsize='large', loc='center left', bbox_to_anchor=(1, 0.5), ncol=1,
                              borderaxespad=0, frameon=False, labelspacing=1)
            for line in leg3.get_lines():
                line.set_linewidth(4.0)

            fig.savefig(folder + scopes[i] + '.png', dpi=200, bbox_extra_artists=(leg,), bbox_inches='tight')
            plt.close(fig)
        except RuntimeError as e:
            print(e.message)
            plt.close()
    return


def plot_wrap(smg_diffs2, smb_diffs2, kst_diffs2, bothg_diffs2, bothb_diffs2, parent):
    """
    Produces a bunch of plots of normed residuals (non-normed currently supressed)
    """
    parent = parent + 'stats/errors_removed/residuals_full/'


    # plot_diffs([smg_diffs2],
    #            ['SM_ws'],
    #            ['mediumblue'],
    #            parent + 'SM_ws/')
    #
    # plot_diffs([smb_diffs2],
    #            ['SM_ps'],
    #            ['lightskyblue'],
    #            parent + 'SM_ps/')
    #
    # plot_diffs([kst_diffs2],
    #            ['KST'],
    #            ['g'],
    #            parent + 'KST/')
    #
    # plot_diffs([bothg_diffs2],
    #            ['Combined_ws'],
    #            ['crimson'],
    #            parent + 'Combined_ws/')
    #
    # plot_diffs([bothb_diffs2],
    #            ['Combined_ps'],
    #            ['coral'],
    #            parent + 'Combined_ps/')
    #
    # plot_diffs([smg_diffs2, smb_diffs2, kst_diffs2, bothg_diffs2, bothb_diffs2],
    #            ['SM_ws', 'SM_ps', 'KST', 'Combined_ws', 'Combined_ps'],
    #            ['mediumblue', 'lightskyblue', 'g', 'crimson', 'coral'],
    #            parent + 'All/')

    plot_diffs_norm([smg_diffs2],
               ['SM_ws'],
               ['mediumblue'],
               parent + 'SM_ws/')

    plot_diffs_norm([smb_diffs2],
               ['SM_ps'],
               ['lightskyblue'],
               parent + 'SM_ps/')

    plot_diffs_norm([kst_diffs2],
               ['KST'],
               ['g'],
               parent + 'KST/')

    plot_diffs_norm([bothg_diffs2],
               ['Combined_ws'],
               ['crimson'],
               parent + 'Combined_ws/')

    plot_diffs_norm([bothb_diffs2],
                    ['Combined_ps'],
                    ['coral'],
                    parent + 'Combined_ps/')

    plot_diffs_norm([smg_diffs2, smb_diffs2, kst_diffs2, bothg_diffs2, bothb_diffs2],
               ['SM_ws', 'SM_ps', 'KST', 'Combined_ws', 'Combined_ps'],
               ['mediumblue', 'lightskyblue', 'g', 'crimson', 'coral'],
               parent + 'All/')
    return

def plot_wrap_er(smg_diffs2, smb_diffs2, kst_diffs2, bothg_diffs2, bothb_diffs2, parent):
    """
    Produces a bunch of plots of normed errors (non-normed currently suppressed)
    """
    parent = parent + 'stats/errors_removed/errors_trimmed/'


    # plot_errors([smg_diffs2],
    #            ['SM_ws'],
    #            ['mediumblue'],
    #            parent + 'SM_ws/')
    #
    # plot_errors([smb_diffs2],
    #            ['SM_ps'],
    #            ['lightskyblue'],
    #            parent + 'SM_ps/')
    #
    # plot_errors([bothb_diffs2],
    #            ['Combined_ps'],
    #            ['coral'],
    #            parent + 'Combined_ps/')
    #
    # plot_errors([smg_diffs2, smb_diffs2, kst_diffs2, bothg_diffs2, bothb_diffs2],
    #            ['SM_ws', 'SM_ps', 'KST', 'Combined_ws', 'Combined_ps'],
    #            ['mediumblue', 'lightskyblue', 'g', 'crimson', 'coral'],
    #            parent + 'All/')
    #
    # plot_errors([kst_diffs2],
    #             ['KST'],
    #             ['g'],
    #             parent + 'KST/')
    #
    # plot_errors([bothg_diffs2],
    #             ['Combined_ws'],
    #             ['crimson'],
    #             parent + 'Combined_ws/')

    plot_errors_norm([smg_diffs2],
               ['SM_ws'],
               ['mediumblue'],
               parent + 'SM_ws/')

    plot_errors_norm([smb_diffs2],
               ['SM_ps'],
               ['lightskyblue'],
               parent + 'SM_ps/')

    plot_errors_norm([kst_diffs2],
                     ['KST'],
                     ['g'],
                     parent + 'KST/')

    plot_errors_norm([bothg_diffs2],
               ['Combined_ws'],
               ['crimson'],
               parent + 'Combined_ws/')

    plot_errors_norm([bothb_diffs2],
                     ['Combined_ps'],
                     ['coral'],
                     parent + 'Combined_ps/')

    plot_errors_norm([smg_diffs2, smb_diffs2, kst_diffs2, bothg_diffs2, bothb_diffs2],
               ['SM_ws', 'SM_ps', 'KST', 'Combined_ws', 'Combined_ps'],
               ['mediumblue', 'lightskyblue', 'g', 'crimson', 'coral'],
               parent + 'All/')
    return


parent = '2019-sets/march-1000/lowz_g10/'


smg_fails, smg_diffs, smg_uncty, uncut_smg = analyse(parent + 'sm_ws/', 'skymapper_ws', wipe_fails=True)
smb_fails, smb_diffs, smb_uncty, uncut_smb = analyse(parent + 'sm_ps/', 'skymapper_ps', wipe_fails=True)

kst_fails, kst_diffs, kst_uncty, uncut_kst  = analyse(parent + 'kst/', 'kst', wipe_fails=True)

bothg_fails, bothg_diffs, bothg_uncty, uncut_bothg = analyse(parent + 'combined_ws/', 'combined_ws', fails=kst_fails, wipe_fails=True)
bothb_fails, bothb_diffs, bothb_uncty, uncut_bothb = analyse(parent + 'combined_ps/', 'combined_ps', fails=kst_fails, wipe_fails=True)

# plot_wrap(smg_diffs, smb_diffs, kst_diffs, bothg_diffs, bothb_diffs, parent)
# plt.close('all')

plot_boxplot([smg_diffs, smb_diffs, kst_diffs, bothg_diffs, bothb_diffs],
               ['SM_ws', 'SM_ps', 'KST', 'Combined_ws', 'Combined_ps'],
               ['mediumblue', 'lightskyblue', 'g', 'crimson', 'coral'],
               parent + 'stats/untrimmed/', mode='residual')

plot_boxplot([smg_uncty, smb_uncty, kst_uncty, bothg_uncty, bothb_uncty],
               ['SM_ws', 'SM_ps', 'KST', 'Combined_ws', 'Combined_ps'],
               ['mediumblue', 'lightskyblue', 'g', 'crimson', 'coral'],
               parent + 'stats/untrimmed/', mode='uncertainty')

plot_trimmed([smg_diffs, smb_diffs, kst_diffs, bothg_diffs, bothb_diffs],
               ['SM_ws', 'SM_ps', 'KST', 'Combined_ws', 'Combined_ps'],
               ['mediumblue', 'lightskyblue', 'g', 'crimson', 'coral'],
               parent + 'stats/trimmed/', mode='residual')

plot_trimmed([smg_uncty, smb_uncty, kst_uncty, bothg_uncty, bothb_uncty],
               ['SM_ws', 'SM_ps', 'KST', 'Combined_ws', 'Combined_ps'],
               ['mediumblue', 'lightskyblue', 'g', 'crimson', 'coral'],
               parent + 'stats/trimmed/', mode='uncertainty')

# plot_wrap_er(smg_diffs_er, smb_diffs_er, kst_diffs_er, bothg_diffs_er, bothb_diffs_er, parent)
# plt.close('all')

# uncertainty_summary([smg_uncty, smb_uncty, kst_uncty, bothg_uncty, bothb_uncty],
#                     [uncut_smg, uncut_smb, uncut_kst, uncut_bothg, uncut_bothb],
#                     ['SM_ws', 'SM_ps', 'KST', 'Combined_ws', 'Combined_ps'],
#                     ['mediumblue', 'lightskyblue', 'g', 'crimson', 'coral'],
#                     parent)
