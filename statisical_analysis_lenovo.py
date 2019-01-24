"""

Run a statistical anaylsis of fitted vs true parameters from model lightcurves of type 1a SN, fitted with sncosmo.
Plots distributions of residuals and errors.
G. Taylor, 2018

- For some reason my ubuntu machine won't upgrade past python 2.7.12, so this is a 'deprecated' version of the
 statistical analysis script, for running on the old machine.  - THIS IS FOR THE WINDOWS LENOVO

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


    # Check error SN and flag:
    error_file = folder + 'error_sn.txt'
    error_set = []
    f = open(error_file, "r")
    lines = f.readlines()

    for x in lines:
        y = x.translate(None, "SN :")
        error_set.append(int(y))
        # Line below doesn't work on ubuntu thinkpad, modified above.
        # error_set.append(int(x.strip('SN : \n')))
    f.close()


    # Import and strip true and fitted params
    fp = folder + 'fitted_parameters.txt'
    tp = folder + 'true_parameters.txt'

    fitted_c = []
    fitted_t0 = []
    fitted_x0 = []
    fitted_x1 = []
    fitted_z = []

    with open(fp, 'rb') as f:
        reader = csv.reader(f, delimiter=' ')
        for row in reader:
            fitted_c.append(float(row[1].replace('c:','')))
            fitted_t0.append(float(row[2].replace('t0:','')))
            fitted_x0.append(float(row[3].replace('x0:','')))
            fitted_x1.append(float(row[4].replace('x1:','')))
            fitted_z.append(float(row[5].replace('z:','')))

    # Double check for missed errors! - need to figure out why these are getting skipped
    for i in range(len(fitted_c)):
        if fitted_c[i] == 0 and fitted_t0[i] == 0 and fitted_x0[i] == 0 and fitted_x1[i] ==0 and fitted_z[i] ==0:
            if i+1 not in error_set:
                error_set.append(i+1)
    error_set.sort()

    true_c = []
    true_t0 = []
    true_x0 = []
    true_x1 = []
    true_z = []

    with open(tp, 'rb') as file:
        reader = csv.reader(file, delimiter=' ')
        for row in reader:
            true_c.append(float(row[1].replace('c:','')))
            true_t0.append(float(row[2].replace('t0:','')))
            true_x0.append(float(row[3].replace('x0:','')))
            true_x1.append(float(row[4].replace('x1:','')))
            true_z.append(float(row[5].replace('z:','')))


    # Calculate residuals

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


    # Create sn_num array
    sn_num = range(1, len(diff_c)+1)
    for i in error_set:
        sn_num[i-1] = 'fit_error'+ str(i)

    # Flag kepler errors (only applicable for combined seeing where kst t0 couldn't be passed)
    for i in fails:
        sn_num[i-1] = 'kst_error' + str(i)

    total_fails = filter(lambda x:x in error_set, fails)

    # Handles failures in both kst and current fit
    for i in total_fails:
        sn_num[i-1] = 'fit_and_kst_error' + str(i)


    # remove fails from data
    if wipe_fails:
        for i in sorted(error_set+fails, reverse=True):
            del sn_num[i-1]
            del diff_c[i-1]
            del diff_t0[i-1]
            del diff_x0[i-1]
            del diff_x1[i-1]
            del diff_z[i-1]


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


    # Writes results
    writefolder = folder + "stats/"
    if not os.path.isdir(writefolder):
        os.makedirs(writefolder)

    if wipe_fails:
        with open(writefolder + 'output_nofails.txt', 'w') as file:
            file.write(table_txt)
    else:
        with open(writefolder + 'output.txt', 'w') as file:
            file.write(table_txt)

    diffs = [diff_c, diff_t0, diff_x0, diff_x1, diff_z]

    return error_set, diffs


def analyse_errors(folder, set, fails=[], wipe_fails=False):
    """
    Extracts errors of SN parameters for one filter set.  Errors are found by an mcmc process in SNCosmo.

    :param folder: path to specific filter set data (e.g. '~\kst\')
    :param set: name of filter set, for labelling purposes (e.g. 'kst')
    :param fails: associated failed SN from another filter set.  For example, if kst has failed to fit a lightcurve to
                    SN10, the combined sm+kst fit will not be able to use kst's fitted t0 value as an initial guess,
                    and the fit will be affected.  SN10 should be counted as a 'fail' for sm+kst, as the fit has
                    underperformed.
    :param wipe_fails: if True, the 'failed' SN will not be included in the output of residuals, showing the best
                        possible performance of the code.  If False, the entire SN set is included, regardless of fit
                        success.
    :return: error set (list of 'failed' SN) and diffs, a table of fitted errors.
    """

    # Naming convention
    if wipe_fails:
        set = set + '_nofails'

    # check error SN and flag:
    error_file = folder + 'error_sn.txt'
    error_set = []
    f = open(error_file, "r")
    lines = f.readlines()

    for x in lines:
        y = x.translate(None, "SN :")
        error_set.append(int(y))
    f.close()


    # Import and strip true and fitted params
    fp = folder + 'fitted_errors.txt'

    fitted_c = []
    fitted_t0 = []
    fitted_x0 = []
    fitted_x1 = []

    with open(fp, 'rb') as f:

        reader = csv.reader(f, delimiter=' ')

        for row in reader:
            # these are actually errors, I should be less lazy
            fitted_c.append(float(row[1].replace('c:','')))
            fitted_t0.append(float(row[2].replace('t0:','')))
            fitted_x0.append(float(row[3].replace('x0:','')))
            fitted_x1.append(float(row[4].replace('x1:','')))

    for i in range(len(fitted_c)):
        if fitted_c[i] == 0 and fitted_t0[i] == 0 and fitted_x0[i] == 0 and fitted_x1[i] ==0:
            if i+1 not in error_set:
                error_set.append(i+1)
    error_set.sort()

    # Create sn_num array
    sn_num = range(1, len(fitted_c)+1)
    for i in error_set:
        sn_num[i-1] = 'fit_error'+ str(i)

    # Flag kepler errors (only applicable for combined seeing where kst t0 couldn't be passed)
    for i in fails:
        sn_num[i-1] = 'kst_error' + str(i)

    total_fails = filter(lambda x:x in error_set, fails)
    # Handles failures in both kst and current fit
    for i in total_fails:
        sn_num[i-1] = 'fit_and_kst_error' + str(i)

    # remove fails from data
    if wipe_fails:
        for i in sorted(error_set+fails, reverse=True):
            del sn_num[i-1]
            del fitted_c[i-1]
            del fitted_t0[i-1]
            del fitted_x0[i-1]
            del fitted_x1[i-1]

    t = PrettyTable()
    t.title = set
    t.add_column('SN', sn_num)
    t.add_column('c-error', fitted_c)
    t.add_column('t0-error', fitted_t0)
    t.add_column('x0-error', fitted_x0)
    t.add_column('x1-error', fitted_x1)

    table_txt = t.get_string()

    writefolder = folder + "stats/"
    if not os.path.isdir(writefolder):
        os.makedirs(writefolder)

    if wipe_fails:
        with open(writefolder + 'error_output_nofails.txt', 'w') as file:
            file.write(table_txt)
    else:
        with open(writefolder + 'error_output.txt', 'w') as file:
            file.write(table_txt)

    diffs = [fitted_c, fitted_t0, fitted_x0, fitted_x1]

    return error_set, diffs


def f(x, a, b, c):
    """
    Equation for a Gaussian
    """
    return a * py.exp((-(x - b) ** 2.0) / (2 * c ** 2))


def plot_diffs(scopes, labels, colour, folder):
#     """
#     Plots histograms and fitted gaussians of residuals for a set of filters, for each parameter.
#     """
#
#     # C - set to [-0.2,0.2] with bins of 0.01
#     crange = 0.2
#     trange = 0.2
#     x1range = 0.2
#     try:
#         figa, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True)
#
#         for i in range(0, len(scopes)):
#             abs_c = scopes[i][0]
#             trimmed_c = [x for x in abs_c if x >= -crange and x <= crange]
#             bins_c = np.arange(min(trimmed_c), max(trimmed_c) + 0.01, 0.01)
#             data = ax1.hist(trimmed_c, bins_c, histtype='step', color=colour[i], label=labels[i]+'(%s fits)'%str(len(trimmed_c)))
#
#             # Generate data from bins as a set of points
#             x = [0.5 * (data[1][t] + data[1][t + 1]) for t in xrange(len(data[1]) - 1)]
#             y = data[0]
#
#             popt, pcov = curve_fit(f, x, y)
#
#             x_fit = py.linspace(-crange, crange, 200)
#             y_fit = f(x_fit, *popt)
#
#             ax2.plot(x_fit, y_fit, lw=2, color=colour[i])
#
#             ax3.hist(trimmed_c, bins_c, histtype='step', color=colour[i], label=labels[i]+'(%s fits)'%str(len(trimmed_c)/2))
#             ax3.plot(x_fit, y_fit, lw=2, color=colour[i])
#
#         if not os.path.isdir(folder):
#             os.makedirs(folder)
#
#         figa.subplots_adjust(hspace=0)
#         plt.setp([a.get_xticklabels() for a in figa.axes[:-1]], visible=False)
#         ax1.set_title('Colour (c) Residuals')
#         plt.xlabel('Residual (true - fitted value of c)')
#         plt.ylabel('Frequency')
#         plt.xlim(-crange, crange)
#         ax1.legend(fontsize = 'x-small')
#         figa.savefig(folder + 'colour.png', dpi=1000)
#         plt.close()
#     except RuntimeError:
#         print('An error occured.')
#         plt.close()
#
#     # t_0 - set to [-0.25,0.25] with bins of 0.02
#     try:
#         figa2, (ax4, ax5, ax6) = plt.subplots(3, sharex=True, sharey=True)
#
#         for i in range(0, len(scopes)):
#             abs_t0 = scopes[i][1]
#             trimmed_t0 = [x for x in abs_t0 if x >= -trange and x <= trange]
#             bins_t0 = np.arange(min(trimmed_t0), max(trimmed_t0) + 0.02, 0.02)
#             data = ax4.hist(trimmed_t0, bins_t0, histtype='step', color=colour[i], label=labels[i]+'(%s fits)'%str(len(trimmed_t0)))
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
#             ax5.plot(x_fit, y_fit, lw=2, color=colour[i])
#
#             ax6.hist(trimmed_t0, bins_t0, histtype='step', color=colour[i], label=labels[i])
#             ax6.plot(x_fit, y_fit, lw=2, color=colour[i])
#
#         if not os.path.isdir(folder):
#             os.makedirs(folder)
#
#         figa2.subplots_adjust(hspace=0)
#         plt.figure(dpi=1200)
#         plt.setp([a.get_xticklabels() for a in figa2.axes[:-1]], visible=False)
#         ax4.set_title(r'Explosion Time ($t_0$) Residuals')
#         plt.xlabel(r'Residual (true - fitted value of $t_0$)')
#         plt.ylabel('Frequency')
#         plt.xlim(-trange, trange)
#         ax4.legend(fontsize = 'x-small')
#         figa2.savefig(folder + 't0.png')
#         plt.close()
#     except RuntimeError:
#         print('An error occured.')
#         plt.close()
#
#     # x_0 - no range because it's a bit garbage.
#     # try:
#     #     figa3, (ax7, ax8, ax9) = plt.subplots(3, sharex=True, sharey=True)
#     #
#     #     for i in range(0, len(scopes)):
#     #         abs_x0 = scopes[i][2]
#     #         bins_x0 = 100
#     #         data = ax7.hist(abs_x0, bins_x0, histtype='step', color=colour[i], label=labels[i]+'(%s fits)'%str(len(abs_x0)))
#     #
#     #         # Generate data from bins as a set of points
#     #         x = [0.5 * (data[1][t] + data[1][t + 1]) for t in xrange(len(data[1]) - 1)]
#     #         y = data[0]
#     #
#     #         popt, pcov = curve_fit(f, x, y)
#     #
#     #         x_fit = py.linspace(min(abs_x0), max(abs_x0), 200)
#     #         y_fit = f(x_fit, *popt)
#     #
#     #         ax8.plot(x_fit, y_fit, lw=2, color=colour[i])
#     #
#     #         ax9.hist(abs_x0, bins_x0, histtype='step', color=colour[i], label=labels[i])
#     #         ax9.plot(x_fit, y_fit, lw=2, color=colour[i])
#     #
#     #     if not os.path.isdir(folder):
#     #         os.makedirs(folder)
#     #
#     #     figa3.subplots_adjust(hspace=0)
#     #     plt.setp([a.get_xticklabels() for a in figa3.axes[:-1]], visible=False)
#     #     ax7.set_title(r'$x_0$ Residuals')
#     #     plt.xlabel(r'Residual (true - fitted value of $x_0$)')
#     #     plt.ylabel('Frequency')
#     #     ax7.legend(fontsize = 'x-small')
#     #     plt.xlim(-1, 1)
#     #     figa3.savefig(folder + 'x0.png')
#     #     plt.close()
#     #     # plt.show()
#     # except RuntimeError:
#     #     print('An error occured.')
#     #     plt.close()
#
#     # x_1 - set to [-0.2, 0.2] with bins of 0.01
#     try:
#         figa4, (ax10, ax11, ax12) = plt.subplots(3, sharex=True, sharey=True)
#
#         for i in range(0, len(scopes)):
#             abs_x1 = scopes[i][3]
#             trimmed_x1 = [x for x in abs_x1 if x >= -x1range and x <= x1range]
#             bins_x1 = np.arange(min(trimmed_x1), max(trimmed_x1) + 0.01, 0.01)
#             data = ax10.hist(trimmed_x1, bins_x1, histtype='step', color=colour[i], label=labels[i]+'(%s fits)'%str(len(trimmed_x1)))
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
#             ax11.plot(x_fit, y_fit, lw=2, color=colour[i])
#
#             ax12.hist(trimmed_x1, bins_x1, histtype='step', color=colour[i], label=labels[i])
#             ax12.plot(x_fit, y_fit, lw=2, color=colour[i])
#
#         if not os.path.isdir(folder):
#             os.makedirs(folder)
#
#         figa4.subplots_adjust(hspace=0)
#         plt.figure(dpi=1200)
#         plt.setp([a.get_xticklabels() for a in figa4.axes[:-1]], visible=False)
#         ax10.set_title(r' Stretch ($x_1$) Residuals')
#         plt.xlabel(r'Residual (true - fitted value of $x_1$)')
#         plt.ylabel('Frequency')
#         ax10.legend(fontsize = 'x-small')
#         plt.xlim(-x1range, x1range)
#         figa4.savefig(folder + 'x1.png')
#         figa4.clf()
#         # plt.show()
#         plt.close()
#     except RuntimeError:
#         print('An error occured.')
#         plt.close()
#
     return

def plot_diffs_norm(scopes, labels, colour, folder):
    """
    Plots histograms and fitted gaussians of residuals for a set of filters, for each parameter.
    """

    crange = 0.2
    trange = 0.2
    x1range = 0.2

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

            x_fit = py.linspace(-crange, crange, 200)
            y_fit = f(x_fit, *popt)

            ax2.plot(x_fit, y_fit, lw=2, color=colour[i], label=labels[i]+'\n $\mu = %s$'%round(popt[1], 5)+'\n $\sigma = %s$'%round(popt[2], 5))

            ax3.hist(trimmed_c, bins_c, histtype='step', density=True, color=colour[i], label=labels[i]+'(%s fits)'%str(len(trimmed_c)/2))
            ax3.plot(x_fit, y_fit, lw=2, color=colour[i])

        if not os.path.isdir(folder):
            os.makedirs(folder)

        figa.subplots_adjust(hspace=0)
        plt.figure(dpi=1200)
        ax2.set_xlim(-0.15, 0.15)
        plt.setp([a.get_xticklabels() for a in figa.axes[:-1]], visible=False)
        ax3.set_xlabel('Residual (true - fitted value of c)')
        ax2.set_ylabel('PDF')
        leg = ax2.legend(fontsize='large', loc='center left', bbox_to_anchor= (1, 0.5), ncol=1,
            borderaxespad=0, frameon=False, labelspacing=1)
        for line in leg.get_lines():
            line.set_linewidth(4.0)
        figa.savefig(folder + 'colour_normed.png', dpi=200, bbox_extra_artists=(leg,), bbox_inches='tight')
        plt.close(figa)
    except RuntimeError:
        print('An error occured.')
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


def uncertainty_summary(errors, scopes, colours, folder):
    """ Calculates and plots summary statistics of uncertainties - a modification of plot_errors?"""

    # Plot success rates

    # Create index array
    bar_indices = np.arange(len(errors))*2
    # Create array of lengths, arbitrarily using c [0] list.
    bar_lengths = [len(x[:][0])/10. for x in errors]
    # Text to go below scope names
    bar_strings = ['\n %s fits'%len(x[:][0]) for x in errors]
    bar_labels = [m+n for m,n in zip(scopes,bar_strings)]
    # Plot bars
    fig = plt.figure()
    bars = plt.bar(bar_indices, bar_lengths, width = 1.0, color=colours, alpha = 0.4)
    # Add counts above the bar graphs
    for rect in bars:
        height = rect.get_height()
        plt.text(rect.get_x() + rect.get_width() / 2.0, height-5, str(height) + '%', ha='center', va='bottom')
    plt.xticks(bar_indices, bar_labels)
    plt.ylabel('Fit success rate (%)')
    fig.savefig(folder + 'success_rate.png', dpi=200, bbox_inches='tight')
    plt.close(fig)

    return

    # For each telescope in scopes

    # calculate mean, median, mode, sigma

    # plot histograms of errors from [0,3sigma]


def plot_wrap(smg_diffs2, smb_diffs2, kst_diffs2, bothg_diffs2, bothb_diffs2, parent):
    """
    Produces a bunch of plots of normed residuals (non-normed currently supressed)
    """
    parent = parent + 'stats/errors_removed/residuals_trimmed/'


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


parent = 'Honours_data_sets/062618/1000setcopy/'


smg_fails2, smg_diffs2 = analyse(parent + 'sm_ws/',
                    'skymapper_ws', wipe_fails=True)
smb_fails2, smb_diffs2 = analyse(parent + 'sm_ps/',
                    'skymapper_ps', wipe_fails=True)

kst_fails2, kst_diffs2 = analyse(parent + 'kst/', 'kst',
                   wipe_fails=True)

bothg_fails2, bothg_diffs2 = analyse(parent + 'combined_ws/',
                   'combined_ws', fails=kst_fails2, wipe_fails=True)
bothb_fails2, bothb_diffs2 = analyse(parent + 'combined_ps/',
                   'combined_ps', fails=kst_fails2, wipe_fails=True)
# plot_wrap(smg_diffs2, smb_diffs2, kst_diffs2, bothg_diffs2, bothb_diffs2, parent)
# plt.close('all')


smg_fails_er, smg_diffs_er = analyse_errors(parent + 'sm_ws/',
                    'skymapper_ws', wipe_fails=True)
smb_fails_er, smb_diffs_er = analyse_errors(parent + 'sm_ps/',
                    'skymapper_ps', wipe_fails=True)
kst_fails_er, kst_diffs_er = analyse_errors(parent + 'kst/', 'kst',
                   wipe_fails=True)
bothg_fails_er, bothg_diffs_er = analyse_errors(parent + 'combined_ws/',
                   'combined_ws', fails=kst_fails_er, wipe_fails=True)
bothb_fails_er, bothb_diffs_er = analyse_errors(parent + 'combined_ps/',
                   'combined_ps', fails=kst_fails_er, wipe_fails=True)
# plot_wrap_er(smg_diffs_er, smb_diffs_er, kst_diffs_er, bothg_diffs_er, bothb_diffs_er, parent)
# plt.close('all')

uncertainty_summary([smg_diffs_er, smb_diffs_er, kst_diffs_er, bothg_diffs_er, bothb_diffs_er], ['SM_ws', 'SM_ps', 'KST', 'Combined_ws', 'Combined_ps'],
               ['mediumblue', 'lightskyblue', 'g', 'crimson', 'coral'], parent)