"""

Run a statistical anaylsis of fitted vs true parameters from model lightcurves of type 1a SN, fitted with sncosmo.
Plots distributions of residuals and errors.
G. Taylor, 2018

"""

import csv
from prettytable import PrettyTable
import os
import matplotlib.pyplot as plt
import pylab as py
import numpy as np
from scipy.optimize import curve_fit
import copy


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
        error_set.append(int(x.strip('SN : \n')))
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
        error_set.append(int(x.strip('SN : \n')))
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
    return a * py.exp(-(x - b) ** 2.0 / (2 * c ** 2))


def plot_diffs(scopes, labels, colour, folder):
    """
    Plots histograms and fitted gaussians of residuals for a set of filters, for each parameter.
    """

    # C - set to [-0.2,0.2] with bins of 0.01
    try:
        figa, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True)

        for i in range(0, len(scopes)):
            abs_c = scopes[i][0]
            trimmed_c = [x for x in abs_c if x >= -0.2 and x <= 0.2]
            bins_c = np.arange(min(trimmed_c), max(trimmed_c) + 0.01, 0.01)
            data = ax1.hist(trimmed_c, bins_c, histtype='step', color=colour[i], label=labels[i]+'(%s fits)'%str(len(trimmed_c)))

            # Generate data from bins as a set of points
            x = [0.5 * (data[1][t] + data[1][t + 1]) for t in xrange(len(data[1]) - 1)]
            y = data[0]

            popt, pcov = curve_fit(f, x, y)

            x_fit = py.linspace(-2, 2, 200)
            y_fit = f(x_fit, *popt)

            ax2.plot(x_fit, y_fit, lw=2, color=colour[i])

            ax3.hist(trimmed_c, bins_c, histtype='step', color=colour[i], label=labels[i]+'(%s fits)'%str(len(trimmed_c)/2))
            ax3.plot(x_fit, y_fit, lw=2, color=colour[i])

        if not os.path.isdir(folder):
            os.makedirs(folder)

        figa.subplots_adjust(hspace=0)
        plt.setp([a.get_xticklabels() for a in figa.axes[:-1]], visible=False)
        ax1.set_title('Colour (c) Residuals')
        plt.xlabel('Residual (true - fitted value of c)')
        plt.ylabel('Frequency')
        plt.xlim(-0.2, 0.2)
        ax1.legend(fontsize = 'x-small')
        figa.savefig(folder + 'colour.png')
        plt.close()
    except RuntimeError:
        print('An error occured.')
        plt.close()

    # t_0 - set to [-0.25,0.25] with bins of 0.02
    try:
        figa2, (ax4, ax5, ax6) = plt.subplots(3, sharex=True, sharey=True)

        for i in range(0, len(scopes)):
            abs_t0 = scopes[i][1]
            trimmed_t0 = [x for x in abs_t0 if x >= -0.25 and x <= .25]
            bins_t0 = np.arange(min(trimmed_t0), max(trimmed_t0) + 0.02, 0.02)
            data = ax4.hist(trimmed_t0, bins_t0, histtype='step', color=colour[i], label=labels[i]+'(%s fits)'%str(len(trimmed_t0)))

            # Generate data from bins as a set of points
            x = [0.5 * (data[1][t] + data[1][t + 1]) for t in xrange(len(data[1]) - 1)]
            y = data[0]

            popt, pcov = curve_fit(f, x, y)

            x_fit = py.linspace(-.25, .25, 200)
            y_fit = f(x_fit, *popt)

            ax5.plot(x_fit, y_fit, lw=2, color=colour[i])

            ax6.hist(trimmed_t0, bins_t0, histtype='step', color=colour[i], label=labels[i])
            ax6.plot(x_fit, y_fit, lw=2, color=colour[i])

        if not os.path.isdir(folder):
            os.makedirs(folder)

        figa2.subplots_adjust(hspace=0)
        plt.setp([a.get_xticklabels() for a in figa2.axes[:-1]], visible=False)
        ax4.set_title(r'Explosion Time ($t_0$) Residuals')
        plt.xlabel(r'Residual (true - fitted value of $t_0$)')
        plt.ylabel('Frequency')
        plt.xlim(-.25, .25)
        ax4.legend(fontsize = 'x-small')
        figa2.savefig(folder + 't0.png')
        plt.close()
    except RuntimeError:
        print('An error occured.')
        plt.close()

    # x_0 - no range because it's a bit garbage.
    # try:
    #     figa3, (ax7, ax8, ax9) = plt.subplots(3, sharex=True, sharey=True)
    #
    #     for i in range(0, len(scopes)):
    #         abs_x0 = scopes[i][2]
    #         bins_x0 = 100
    #         data = ax7.hist(abs_x0, bins_x0, histtype='step', color=colour[i], label=labels[i]+'(%s fits)'%str(len(abs_x0)))
    #
    #         # Generate data from bins as a set of points
    #         x = [0.5 * (data[1][t] + data[1][t + 1]) for t in xrange(len(data[1]) - 1)]
    #         y = data[0]
    #
    #         popt, pcov = curve_fit(f, x, y)
    #
    #         x_fit = py.linspace(min(abs_x0), max(abs_x0), 200)
    #         y_fit = f(x_fit, *popt)
    #
    #         ax8.plot(x_fit, y_fit, lw=2, color=colour[i])
    #
    #         ax9.hist(abs_x0, bins_x0, histtype='step', color=colour[i], label=labels[i])
    #         ax9.plot(x_fit, y_fit, lw=2, color=colour[i])
    #
    #     if not os.path.isdir(folder):
    #         os.makedirs(folder)
    #
    #     figa3.subplots_adjust(hspace=0)
    #     plt.setp([a.get_xticklabels() for a in figa3.axes[:-1]], visible=False)
    #     ax7.set_title(r'$x_0$ Residuals')
    #     plt.xlabel(r'Residual (true - fitted value of $x_0$)')
    #     plt.ylabel('Frequency')
    #     ax7.legend(fontsize = 'x-small')
    #     plt.xlim(-1, 1)
    #     figa3.savefig(folder + 'x0.png')
    #     plt.close()
    #     # plt.show()
    # except RuntimeError:
    #     print('An error occured.')
    #     plt.close()

    # x_1 - set to [-0.2, 0.2] with bins of 0.01
    try:
        figa4, (ax10, ax11, ax12) = plt.subplots(3, sharex=True, sharey=True)

        for i in range(0, len(scopes)):
            abs_x1 = scopes[i][3]
            trimmed_x1 = [x for x in abs_x1 if x >= -0.2 and x <= 0.2]
            bins_x1 = np.arange(min(trimmed_x1), max(trimmed_x1) + 0.01, 0.01)
            data = ax10.hist(trimmed_x1, bins_x1, histtype='step', color=colour[i], label=labels[i]+'(%s fits)'%str(len(trimmed_x1)))

            # Generate data from bins as a set of points
            x = [0.5 * (data[1][t] + data[1][t + 1]) for t in xrange(len(data[1]) - 1)]
            y = data[0]

            popt, pcov = curve_fit(f, x, y)

            x_fit = py.linspace(-0.5, 0.5, 200)
            y_fit = f(x_fit, *popt)

            ax11.plot(x_fit, y_fit, lw=2, color=colour[i])

            ax12.hist(trimmed_x1, bins_x1, histtype='step', color=colour[i], label=labels[i])
            ax12.plot(x_fit, y_fit, lw=2, color=colour[i])

        if not os.path.isdir(folder):
            os.makedirs(folder)

        figa4.subplots_adjust(hspace=0)
        plt.setp([a.get_xticklabels() for a in figa4.axes[:-1]], visible=False)
        ax10.set_title(r' Stretch ($x_1$) Residuals')
        plt.xlabel(r'Residual (true - fitted value of $x_1$)')
        plt.ylabel('Frequency')
        ax10.legend(fontsize = 'x-small')
        plt.xlim(-0.2, 0.2)
        figa4.savefig(folder + 'x1.png')
        figa4.clf()
        # plt.show()
        plt.close()
    except RuntimeError:
        print('An error occured.')
        plt.close()

    return

def plot_diffs_norm(scopes, labels, colour, folder):
    """
    Plots histograms and fitted gaussians of residuals for a set of filters, for each parameter.
    """

    # C - set to [-0.2,0.2] with bins of 0.01
    try:
        figa, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True)

        for i in range(0, len(scopes)):
            abs_c = scopes[i][0]
            trimmed_c = [x for x in abs_c if x >= -0.2 and x <= 0.2]
            bins_c = np.arange(min(trimmed_c), max(trimmed_c) + 0.01, 0.01)
            data = ax1.hist(trimmed_c, bins_c, histtype='step', density=True, color=colour[i], label=labels[i]+'(%s fits)'%str(len(trimmed_c)))

            # Generate data from bins as a set of points
            x = [0.5 * (data[1][t] + data[1][t + 1]) for t in xrange(len(data[1]) - 1)]
            y = data[0]

            popt, pcov = curve_fit(f, x, y)

            x_fit = py.linspace(-2, 2, 200)
            y_fit = f(x_fit, *popt)

            ax2.plot(x_fit, y_fit, lw=2, color=colour[i])

            ax3.hist(trimmed_c, bins_c, histtype='step', density=True, color=colour[i], label=labels[i]+'(%s fits)'%str(len(trimmed_c)/2))
            ax3.plot(x_fit, y_fit, lw=2, color=colour[i])

        if not os.path.isdir(folder):
            os.makedirs(folder)

        figa.subplots_adjust(hspace=0)
        plt.setp([a.get_xticklabels() for a in figa.axes[:-1]], visible=False)
        ax1.set_title('NORMALIZED Colour (c) Residuals')
        plt.xlabel('Residual (true - fitted value of c)')
        plt.ylabel('PDF')
        plt.xlim(-0.2, 0.2)
        ax1.legend(fontsize = 'x-small')
        figa.savefig(folder + 'colour_normed.png')
        plt.close()
    except RuntimeError:
        print('An error occured.')
        plt.close()

    # t_0 - set to [-0.25,0.25] with bins of 0.0.02
    try:
        figa2, (ax4, ax5, ax6) = plt.subplots(3, sharex=True, sharey=True)

        for i in range(0, len(scopes)):
            abs_t0 = scopes[i][1]
            trimmed_t0 = [x for x in abs_t0 if x >= -0.25 and x <= .25]
            bins_t0 = np.arange(min(trimmed_t0), max(trimmed_t0) + 0.02, 0.02)
            data = ax4.hist(trimmed_t0, bins_t0, histtype='step', density=True, color=colour[i], label=labels[i]+'(%s fits)'%str(len(trimmed_t0)))

            # Generate data from bins as a set of points
            x = [0.5 * (data[1][t] + data[1][t + 1]) for t in xrange(len(data[1]) - 1)]
            y = data[0]

            popt, pcov = curve_fit(f, x, y)

            x_fit = py.linspace(-.25, .25, 200)
            y_fit = f(x_fit, *popt)

            ax5.plot(x_fit, y_fit, lw=2, color=colour[i])

            ax6.hist(trimmed_t0, bins_t0, histtype='step', density=True, color=colour[i], label=labels[i])
            ax6.plot(x_fit, y_fit, lw=2, color=colour[i])

        if not os.path.isdir(folder):
            os.makedirs(folder)

        figa2.subplots_adjust(hspace=0)
        plt.setp([a.get_xticklabels() for a in figa2.axes[:-1]], visible=False)
        ax4.set_title(r'NORMALIZED Explosion Time ($t_0$) Residuals')
        plt.xlabel(r'Residual (true - fitted value of $t_0$)')
        plt.ylabel('PDF')
        plt.xlim(-.25, .25)
        ax4.legend(fontsize = 'x-small')
        figa2.savefig(folder + 't0_normed.png')
        plt.close()
    except RuntimeError:
        print('An error occured.')
        plt.close()

    # x_1 - set to [-0.2, 0.2] with bins of 0.01
    try:
        figa4, (ax10, ax11, ax12) = plt.subplots(3, sharex=True, sharey=True)

        minn = 10000
        maxx = -10000
        for i in range(0, len(scopes)):
            abs_x1 = scopes[i][3]
            trimmed_x1 = [x for x in abs_x1 if x >= -0.2 and x <= 0.2]
            bins_x1 = np.arange(min(trimmed_x1), max(trimmed_x1) + 0.01, 0.01)
            data = ax10.hist(trimmed_x1, bins_x1, histtype='step', density=True, color=colour[i], label=labels[i]+'(%s fits)'%str(len(trimmed_x1)))

            # Generate data from bins as a set of points
            x = [0.5 * (data[1][t] + data[1][t + 1]) for t in xrange(len(data[1]) - 1)]
            y = data[0]

            popt, pcov = curve_fit(f, x, y)

            x_fit = py.linspace(-0.5, 0.5, 200)
            y_fit = f(x_fit, *popt)

            ax11.plot(x_fit, y_fit, lw=2, color=colour[i])

            ax12.hist(trimmed_x1, bins_x1, histtype='step', density=True, color=colour[i], label=labels[i])
            ax12.plot(x_fit, y_fit, lw=2, color=colour[i])

        if not os.path.isdir(folder):
            os.makedirs(folder)

        figa4.subplots_adjust(hspace=0)
        plt.setp([a.get_xticklabels() for a in figa4.axes[:-1]], visible=False)
        ax10.set_title(r'NORMALIZED Stretch ($x_1$) Residuals')
        plt.xlabel(r'Residual (true - fitted value of $x_1$)')
        plt.ylabel('PDF')
        ax10.legend(fontsize = 'x-small')
        plt.xlim(-0.2, 0.2)
        figa4.savefig(folder + 'x1_normed.png')
        figa4.clf()
        # plt.show()
        plt.close()
    except RuntimeError:
        print('An error occured.')
        plt.close()

    return

def plot_errors(scopes, labels, colour, folder):
    """
    Plots histograms and fitted gaussians of residuals for a set of filters, for each parameter.
    """

    # C - set to [-0.2,0.2] with bins of 0.01
    try:
        figa, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True)
        for i in range(0, len(scopes)):
            abs_c = copy.deepcopy(scopes[i][0])
            c_flips = [-x for x in abs_c]
            abs_c.extend(c_flips)
            trimmed_c = [x for x in abs_c if x >= -0.2 and x <= 0.2]
            bins_c = np.arange(min(trimmed_c), max(trimmed_c) + 0.01, 0.01)
            data = ax1.hist(trimmed_c, bins_c, histtype='step', color=colour[i], label=labels[i]+'(%s fits)'%str(len(trimmed_c)/2))

            # Generate data from bins as a set of points
            x = [0.5 * (data[1][t] + data[1][t + 1]) for t in xrange(len(data[1]) - 1)]
            y = data[0]

            popt, pcov = curve_fit(f, x, y)

            x_fit = py.linspace(-.2, .2, 200)#(min(abs_c), max(abs_c), 200)
            y_fit = f(x_fit, *popt)

            ax2.plot(x_fit, y_fit, lw=2, color=colour[i])

            ax3.hist(trimmed_c, bins_c, histtype='step', color=colour[i], label=labels[i])
            ax3.plot(x_fit, y_fit, lw=2, color=colour[i])

        if not os.path.isdir(folder):
            os.makedirs(folder)

        figa.subplots_adjust(hspace=0)
        plt.setp([a.get_xticklabels() for a in figa.axes[:-1]], visible=False)
        ax1.set_title(r'Colour (c) Errors')
        plt.xlabel('c error')
        plt.ylabel('Frequency')
        ax1.legend(fontsize = 'x-small')
        figa.savefig(folder + 'colour.png')
        plt.close()
    except RuntimeError as e:
        print(e.message)
        plt.close()

    # t_0 - set to [-1, 1] with bins of 0.1
    try:
        figa2, (ax4, ax5, ax6) = plt.subplots(3, sharex=True, sharey=True)

        for i in range(0, len(scopes)):
            abs_t0 = copy.deepcopy(scopes[i][1])
            t0_flips = [-x for x in abs_t0]
            abs_t0.extend(t0_flips)
            trimmed_t0 = [x for x in abs_t0 if x >= -0.25 and x <= 0.25]
            bins_t0 = np.arange(min(trimmed_t0), max(trimmed_t0) + 0.02, 0.02)
            data = ax4.hist(trimmed_t0, bins_t0, histtype='step', color=colour[i], label=labels[i]+'(%s fits)'%str(len(trimmed_t0)/2))

            # Generate data from bins as a set of points
            x = [0.5 * (data[1][t] + data[1][t + 1]) for t in xrange(len(data[1]) - 1)]
            y = data[0]

            popt, pcov = curve_fit(f, x, y)

            x_fit = py.linspace(min(trimmed_t0), max(trimmed_t0), 200)
            y_fit = f(x_fit, *popt)

            ax5.plot(x_fit, y_fit, lw=2, color=colour[i])

            ax6.hist(trimmed_t0, bins_t0, histtype='step', color=colour[i], label=labels[i])
            ax6.plot(x_fit, y_fit, lw=2, color=colour[i])

        if not os.path.isdir(folder):
            os.makedirs(folder)

        figa2.subplots_adjust(hspace=0)
        plt.setp([a.get_xticklabels() for a in figa2.axes[:-1]], visible=False)
        ax4.set_title(r'Explosion Time ($t_0$) Errors')
        plt.xlabel(r'$t_0$ error')
        plt.ylabel('Frequency')
        ax4.legend(fontsize = 'x-small')
        figa2.savefig(folder + 't0.png')
        plt.close()
    except RuntimeError:
        print('An error occured.')
        plt.close()

    # x_0 - unimportant.
    # try:
    #     figa3, (ax7, ax8, ax9) = plt.subplots(3, sharex=True, sharey=True)
    #
    #     for i in range(0, len(scopes)):
    #         abs_x0 = scopes[i][2]
    #         x0_flips = [-x for x in abs_x0]
    #         abs_x0.extend(x0_flips)
    #         trimmed_x0 = [x for x in abs_x0 if x >= -5 and x <= 5]
    #         bins_x0 = 100
    #         data = ax7.hist(abs_x0, bins_x0, histtype='step', color=colour[i], label=labels[i]+'(%s fits)'%str(len(abs_x0)/2))
    #
    #         # Generate data from bins as a set of points
    #         x = [0.5 * (data[1][t] + data[1][t + 1]) for t in xrange(len(data[1]) - 1)]
    #         y = data[0]
    #
    #         popt, pcov = curve_fit(f, x, y)
    #
    #         x_fit = py.linspace(min(abs_x0), max(abs_x0), 200)
    #         y_fit = f(x_fit, *popt)
    #
    #         ax8.plot(x_fit, y_fit, lw=2, color=colour[i])
    #
    #         ax9.hist(abs_x0, bins_x0, histtype='step', color=colour[i], label=labels[i])
    #         ax9.plot(x_fit, y_fit, lw=2, color=colour[i])
    #
    #     if not os.path.isdir(folder):
    #         os.makedirs(folder)
    #
    #
    #     figa3.subplots_adjust(hspace=0)
    #     plt.setp([a.get_xticklabels() for a in figa3.axes[:-1]], visible=False)
    #     ax7.set_title(r'$x_0$ Errors')
    #     plt.xlabel(r'$x_0$ Error')
    #     plt.ylabel('Frequency')
    #     ax7.legend(fontsize = 'x-small')
    #     figa3.savefig(folder + 'x0.png')
    #     plt.close()
    #     # plt.show()
    # except RuntimeError:
    #     print('An error occured.')
    #     plt.close()



    # x_1 - set to [-0.5, 0.5] with bins of 0.01 (nbins=100)
    try:
        figa4, (ax10, ax11, ax12) = plt.subplots(3, sharex=True, sharey=True)

        for i in range(0, len(scopes)):
            abs_x1 = copy.deepcopy(scopes[i][3])
            x1_flips = [-x for x in abs_x1]
            abs_x1.extend(x1_flips)
            trimmed_x1 = [x for x in abs_x1 if x >= -0.5 and x <= 0.5]
            bins_x1 = np.arange(min(trimmed_x1), max(trimmed_x1) + 0.01, 0.01)
            data = ax10.hist(trimmed_x1, bins_x1, histtype='step', color=colour[i], label=labels[i]+'(%s fits)'%str(len(trimmed_x1)/2))

            # Generate data from bins as a set of points
            x = [0.5 * (data[1][t] + data[1][t + 1]) for t in xrange(len(data[1]) - 1)]
            y = data[0]

            popt, pcov = curve_fit(f, x, y)

            x_fit = py.linspace(min(trimmed_x1), max(trimmed_x1), 200)
            y_fit = f(x_fit, *popt)

            ax11.plot(x_fit, y_fit, lw=2, color=colour[i])

            ax12.hist(trimmed_x1, bins_x1, histtype='step', color=colour[i], label=labels[i])
            ax12.plot(x_fit, y_fit, lw=2, color=colour[i])

        if not os.path.isdir(folder):
            os.makedirs(folder)

        figa4.subplots_adjust(hspace=0)
        plt.setp([a.get_xticklabels() for a in figa4.axes[:-1]], visible=False)
        ax10.set_title(r'Fitted $x_1$ Errors')
        plt.xlabel(r'Fitted $x_1$ Error')
        plt.ylabel('Frequency')
        ax10.legend(fontsize = 'x-small')
        figa4.savefig(folder + 'x1.png')
        figa4.clf()
        # plt.show()
        plt.close()
    except RuntimeError:
        print('An error occured.')
        plt.close()

    return

def plot_errors_norm(scopes, labels, colour, folder):
    """
    Plots histograms and fitted gaussians of residuals for a set of filters, for each parameter.
    """

    # C - set to [-0.2,0.2] with bins of 0.01
    try:
        figa, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True)
        for i in range(0, len(scopes)):
            abs_c = copy.deepcopy(scopes[i][0])
            c_flips = [-x for x in abs_c]
            abs_c.extend(c_flips)
            trimmed_c = [x for x in abs_c if x >= -0.2 and x <= 0.2]
            bins_c = np.arange(min(trimmed_c), max(trimmed_c) + 0.01, 0.01)
            data = ax1.hist(trimmed_c, bins_c, histtype='step', density=True, color=colour[i], label=labels[i]+'(%s fits)'%str(len(trimmed_c)/2))

            # Generate data from bins as a set of points
            x = [0.5 * (data[1][t] + data[1][t + 1]) for t in xrange(len(data[1]) - 1)]
            y = data[0]

            popt, pcov = curve_fit(f, x, y)

            x_fit = py.linspace(-.2, .2, 200)#(min(abs_c), max(abs_c), 200)
            y_fit = f(x_fit, *popt)

            ax2.plot(x_fit, y_fit, lw=2, color=colour[i])

            ax3.hist(trimmed_c, bins_c, histtype='step', density=True, color=colour[i], label=labels[i])
            ax3.plot(x_fit, y_fit, lw=2, color=colour[i])

        if not os.path.isdir(folder):
            os.makedirs(folder)

        figa.subplots_adjust(hspace=0)
        plt.setp([a.get_xticklabels() for a in figa.axes[:-1]], visible=False)
        ax1.set_title(r'NORMALIZED Colour (c) Errors')
        plt.xlabel('c error')
        plt.ylabel('PDF')
        ax1.legend(fontsize = 'x-small')
        figa.savefig(folder + 'colour_normed.png')
        plt.close()
    except RuntimeError as e:
        print(e.message)
        plt.close()

    # t_0 - set to [-1, 1] with bins of 0.1
    try:
        figa2, (ax4, ax5, ax6) = plt.subplots(3, sharex=True, sharey=True)

        for i in range(0, len(scopes)):
            abs_t0 = copy.deepcopy(scopes[i][1])
            t0_flips = [-x for x in abs_t0]
            abs_t0.extend(t0_flips)
            trimmed_t0 = [x for x in abs_t0 if x >= -0.25 and x <= 0.25]
            bins_t0 = np.arange(min(trimmed_t0), max(trimmed_t0) + 0.02, 0.02)
            data = ax4.hist(trimmed_t0, bins_t0, histtype='step', density=True, color=colour[i], label=labels[i]+'(%s fits)'%str(len(trimmed_t0)/2))

            # Generate data from bins as a set of points
            x = [0.5 * (data[1][t] + data[1][t + 1]) for t in xrange(len(data[1]) - 1)]
            y = data[0]

            popt, pcov = curve_fit(f, x, y)

            x_fit = py.linspace(min(trimmed_t0), max(trimmed_t0), 200)
            y_fit = f(x_fit, *popt)

            ax5.plot(x_fit, y_fit, lw=2, color=colour[i])

            ax6.hist(trimmed_t0, bins_t0, histtype='step', density=True, color=colour[i], label=labels[i])
            ax6.plot(x_fit, y_fit, lw=2, color=colour[i])

        if not os.path.isdir(folder):
            os.makedirs(folder)

        figa2.subplots_adjust(hspace=0)
        plt.setp([a.get_xticklabels() for a in figa2.axes[:-1]], visible=False)
        ax4.set_title(r'NORMALIZED Explosion Time ($t_0$) Errors')
        plt.xlabel(r'$t_0$ error')
        plt.ylabel('PDF')
        ax4.legend(fontsize = 'x-small')
        figa2.savefig(folder + 't0_normed.png')
        plt.close()
    except RuntimeError:
        print('An error occured.')
        plt.close()

    # x_1 - set to [-0.5, 0.5] with bins of 0.01 (nbins=100)
    try:
        figa4, (ax10, ax11, ax12) = plt.subplots(3, sharex=True, sharey=True)

        for i in range(0, len(scopes)):
            abs_x1 = copy.deepcopy(scopes[i][3])
            x1_flips = [-x for x in abs_x1]
            abs_x1.extend(x1_flips)
            trimmed_x1 = [x for x in abs_x1 if x >= -0.5 and x <= 0.5]
            bins_x1 = np.arange(min(trimmed_x1), max(trimmed_x1) + 0.01, 0.01)
            data = ax10.hist(trimmed_x1, bins_x1, histtype='step', density=True, color=colour[i], label=labels[i]+'(%s fits)'%str(len(trimmed_x1)/2))

            # Generate data from bins as a set of points
            x = [0.5 * (data[1][t] + data[1][t + 1]) for t in xrange(len(data[1]) - 1)]
            y = data[0]

            popt, pcov = curve_fit(f, x, y)

            x_fit = py.linspace(min(trimmed_x1), max(trimmed_x1), 200)
            y_fit = f(x_fit, *popt)

            ax11.plot(x_fit, y_fit, lw=2, color=colour[i])

            ax12.hist(trimmed_x1, bins_x1, histtype='step',density=True, color=colour[i], label=labels[i])
            ax12.plot(x_fit, y_fit, lw=2, color=colour[i])

        if not os.path.isdir(folder):
            os.makedirs(folder)

        figa4.subplots_adjust(hspace=0)
        plt.setp([a.get_xticklabels() for a in figa4.axes[:-1]], visible=False)
        ax10.set_title(r'NORMALIZED Fitted $x_1$ Errors')
        plt.xlabel(r'Fitted $x_1$ Error')
        plt.ylabel('PDF')
        ax10.legend(fontsize = 'x-small')
        figa4.savefig(folder + 'x1_normed.png')
        figa4.clf()
        # plt.show()
        plt.close()
    except RuntimeError:
        print('An error occured.')
        plt.close()

    return

def plot_wrap(smb_diffs2, kst_diffs2, bothb_diffs2, parent):
    parent = parent + 'stats/errors_removed/residuals/'

    plot_diffs([smb_diffs2],
               ['SM'],
               ['b'],
               parent + 'SM/')

    plot_diffs_norm([smb_diffs2],
               ['SM'],
               ['b'],
               parent + 'SM/')

    plot_diffs([kst_diffs2],
               ['KST'],
               ['g'],
               parent + 'KST/')

    plot_diffs_norm([kst_diffs2],
               ['KST'],
               ['g'],
               parent + 'KST/')

    plot_diffs([bothb_diffs2],
               ['Combined'],
               ['y'],
               parent + 'Combined/')

    plot_diffs_norm([bothb_diffs2],
               ['Combined'],
               ['y'],
               parent + 'Combined/')

    plot_diffs([smb_diffs2, kst_diffs2, bothb_diffs2],
               ['SM', 'KST', 'Combined'],
               ['b', 'g', 'y'],
               parent + 'All/')

    plot_diffs_norm([smb_diffs2, kst_diffs2, bothb_diffs2],
               ['SM', 'KST', 'Combined'],
               ['b', 'g', 'y'],
               parent + 'All/')
    return

def plot_wrap_er(smb_diffs2, kst_diffs2, bothb_diffs2, parent):
    parent = parent + 'stats/errors_removed/errors/'

    plot_errors([smb_diffs2],
               ['SM'],
               ['b'],
               parent + 'SM/')

    plot_errors_norm([smb_diffs2],
                ['SM'],
                ['b'],
                parent + 'SM/')

    plot_errors([kst_diffs2],
               ['KST'],
               ['g'],
               parent + 'KST/')

    plot_errors_norm([kst_diffs2],
                ['KST'],
                ['g'],
                parent + 'KST/')

    plot_errors([bothb_diffs2],
               ['Combined'],
               ['y'],
               parent + 'Combined/')

    plot_errors_norm([bothb_diffs2],
               ['Combined'],
               ['y'],
               parent + 'Combined/')

    plot_errors([smb_diffs2, kst_diffs2, bothb_diffs2],
               ['SM', 'KST', 'Combined'],
               ['b', 'g', 'y'],
               parent + 'All/')

    plot_errors_norm([smb_diffs2, kst_diffs2, bothb_diffs2],
                ['SM', 'KST', 'Combined'],
                ['b', 'g', 'y'],
                parent + 'All/')
    return


parent = 'Honours_data_sets/060818/ws_sm/sn300/'


smb_fails2, smb_diffs2 = analyse(parent + 'sm-ws/',
                    'skymapper', wipe_fails=True)
kst_fails2, kst_diffs2 = analyse(parent + 'kst/', 'kst',
                   wipe_fails=True)
bothb_fails2, bothb_diffs2 = analyse(parent + 'combined-ws/',
                   'combined', fails=kst_fails2, wipe_fails=True)
# plot_wrap(smb_diffs2, kst_diffs2, bothb_diffs2, parent)


smb_fails_er, smb_diffs_er = analyse_errors(parent + 'sm-ws/',
                    'skymapper', wipe_fails=True)
kst_fails_er, kst_diffs_er = analyse_errors(parent + 'kst/', 'kst',
                   wipe_fails=True)
bothb_fails_er, bothb_diffs_er = analyse_errors(parent + 'combined-ws/',
                   'combined', fails=kst_fails_er, wipe_fails=True)
# plot_wrap_er(smb_diffs_er, kst_diffs_er, bothb_diffs_er,  parent)