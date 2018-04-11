""" Run a statistical anaylsis of fitted vs true parameters from model lightcurves of type 1a SN, fitted with sncosmo."""
import csv
from prettytable import PrettyTable
import os
import matplotlib.pyplot as plt
import pylab as py
from scipy.stats import norm
import matplotlib.mlab as mlab
import numpy as np
import win32api
from scipy.optimize import curve_fit


def analyse(folder, set, fails=[], wipe_fails=False):

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

    # For percentage difference not using right now)
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

    t = PrettyTable()
    t.title = set
    t.add_column('SN', sn_num)
    t.add_column('c-diff-percent', diff_c)
    t.add_column('t0-diff-days', diff_t0)
    t.add_column('x0-diff-percent', diff_x0)
    t.add_column('x1-diff-percent', diff_x1)
    t.add_column('z-diff-percent', diff_z)

    table_txt = t.get_string()

    writefolder = folder + "stats\\"
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

def abs_and_weight(list):
    # Not currently using absolute values
    #abs_list = map(abs, list)
    abs_list = list
    weights = np.ones_like(abs_list) / float(len(list))
    return abs_list, weights

# Equation for Gaussian
def f(x, a, b, c):
    return a * py.exp(-(x - b) ** 2.0 / (2 * c ** 2))

def plot_diffs_2(scopes, labels, colour, folder):

    figa, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True)
    min_c = 10000
    max_c = -10000
    for i in range(0, len(scopes)):
        abs_c = scopes[i][0]
        trimmed_c = [x for x in abs_c if x >= -2 and x <= 2]
        bins_c = 40
        data = ax1.hist(trimmed_c, bins_c, histtype='step', color=colour[i], label=labels[i])

        # Generate data from bins as a set of points
        x = [0.5 * (data[1][t] + data[1][t + 1]) for t in xrange(len(data[1]) - 1)]
        y = data[0]

        popt, pcov = curve_fit(f, x, y)

        x_fit = py.linspace(-2, 2, 200)
        y_fit = f(x_fit, *popt)

        ax2.plot(x_fit, y_fit, lw=2, color=colour[i])

        ax3.hist(trimmed_c, bins_c, histtype='step', color=colour[i], label=labels[i])
        ax3.plot(x_fit, y_fit, lw=2, color=colour[i])

        if min_c >= min(trimmed_c):
            min_c = min(trimmed_c)
        if max_c <= max(trimmed_c):
            max_c = max(trimmed_c)

    if not os.path.isdir(folder):
        os.makedirs(folder)

    figa.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in figa.axes[:-1]], visible=False)
    ax1.set_title('Colour Residuals')
    plt.xlabel('Residual (true - fitted value of c)')
    plt.ylabel('Frequency')
    plt.xlim(min_c,max_c)
    ax1.legend()
    figa.savefig(folder + 'colour.png')
    plt.close()


    figa2, (ax4, ax5, ax6) = plt.subplots(3, sharex=True, sharey=True)

    minn = 10000
    maxx = -10000
    for i in range(0, len(scopes)):
        abs_t0 = scopes[i][1]
        trimmed_t0 = [x for x in abs_t0 if x >= -5 and x <= 5]
        bins_t0 = 40
        data = ax4.hist(trimmed_t0, bins_t0, histtype='step', color=colour[i], label=labels[i])

        # Generate data from bins as a set of points
        x = [0.5 * (data[1][t] + data[1][t + 1]) for t in xrange(len(data[1]) - 1)]
        y = data[0]

        popt, pcov = curve_fit(f, x, y)

        x_fit = py.linspace(-5, 5, 200)
        y_fit = f(x_fit, *popt)

        ax5.plot(x_fit, y_fit, lw=2, color=colour[i])

        ax6.hist(trimmed_t0, bins_t0, histtype='step', color=colour[i], label=labels[i])
        ax6.plot(x_fit, y_fit, lw=2, color=colour[i])

        if minn >= min(trimmed_t0):
            minn = min(trimmed_t0)
        if maxx <= max(trimmed_t0):
            maxx = max(trimmed_t0)

    if not os.path.isdir(folder):
        os.makedirs(folder)

    figa2.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in figa.axes[:-1]], visible=False)
    ax4.set_title('Explosion Time Residuals')
    plt.xlabel('Residual (true - fitted value of t0)')
    plt.ylabel('Probability')
    plt.xlim(minn, maxx)
    ax3.legend()
    figa2.savefig(folder + 't0.png')
    plt.close()

    #
    figa3, (ax7, ax8, ax9) = plt.subplots(3, sharex=True, sharey=True)

    minn = 10000
    maxx = -10000
    for i in range(0, len(scopes)):
        abs_x0 = scopes[i][2]
        trimmed_x0 = [x for x in abs_x0 if x >= -5 and x <= 5]
        bins_x0 = 100
        data = ax7.hist(abs_x0, bins_x0, histtype='step', color=colour[i], label=labels[i])

        # Generate data from bins as a set of points
        x = [0.5 * (data[1][t] + data[1][t + 1]) for t in xrange(len(data[1]) - 1)]
        y = data[0]

        popt, pcov = curve_fit(f, x, y)

        x_fit = py.linspace(min(abs_x0), max(abs_x0), 200)
        y_fit = f(x_fit, *popt)

        ax8.plot(x_fit, y_fit, lw=2, color=colour[i])

        ax9.hist(abs_x0, bins_x0, histtype='step', color=colour[i], label=labels[i])
        ax9.plot(x_fit, y_fit, lw=2, color=colour[i])

    if not os.path.isdir(folder):
        os.makedirs(folder)


    figa3.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in figa.axes[:-1]], visible=False)
    ax7.set_title('x0 Residuals')
    plt.xlabel('Residual (true - fitted value of x0)')
    plt.ylabel('Probability')
    ax7.legend()
    plt.xlim(-1, 1)
    figa3.savefig(folder + 'x0.png')
    plt.close()
    # plt.show()


    figa4, (ax10, ax11, ax12) = plt.subplots(3, sharex=True, sharey=True)

    minn = 10000
    maxx = -10000
    for i in range(0, len(scopes)):
        abs_x1 = scopes[i][3]
        trimmed_x1 = [x for x in abs_x1 if x >= -2 and x <= 2]
        bins_x1 = 16
        data = ax10.hist(trimmed_x1, bins_x1, histtype='step', color=colour[i], label=labels[i])

        # Generate data from bins as a set of points
        x = [0.5 * (data[1][t] + data[1][t + 1]) for t in xrange(len(data[1]) - 1)]
        y = data[0]

        popt, pcov = curve_fit(f, x, y)

        x_fit = py.linspace(-2, 2, 200)
        y_fit = f(x_fit, *popt)

        ax11.plot(x_fit, y_fit, lw=2, color=colour[i])

        ax12.hist(trimmed_x1, bins_x1, histtype='step', color=colour[i], label=labels[i])
        ax12.plot(x_fit, y_fit, lw=2, color=colour[i])

        if minn >= min(trimmed_x1):
            minn = min(trimmed_x1)
        if maxx <= max(trimmed_x1):
            maxx = max(trimmed_x1)

    if not os.path.isdir(folder):
        os.makedirs(folder)

    figa4.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in figa.axes[:-1]], visible=False)
    ax10.set_title('x1 Residuals')
    plt.xlabel('Residual (true - fitted value of x1)')
    plt.ylabel('Probability')
    ax10.legend()
    plt.xlim(minn, maxx)
    figa4.savefig(folder + 'x1.png')
    figa4.clf()
    # plt.show()

    return


smb_fails2, smb_diffs2 = analyse('Honours_data_sets/5_270218/1_stat_sample/Kepler_6hours/SM_5day/vObs_2/100SN_3/sm_bad_seeing/',
                    'sm_bad_seeing', wipe_fails=True)
smg_fails2, smg_diffs2 = analyse('Honours_data_sets/5_270218/1_stat_sample/Kepler_6hours/SM_5day/vObs_2/100SN_3/sm_good_seeing/',
                    'sm_good_seeing', wipe_fails=True)
kst_fails2, kst_diffs2 = analyse('Honours_data_sets/5_270218/1_stat_sample/Kepler_6hours/SM_5day/vObs_2/100SN_3/kst/', 'kst',
                   wipe_fails=True)
bothb_fails2, bothb_diffs2 = analyse('Honours_data_sets/5_270218/1_stat_sample/Kepler_6hours/SM_5day/vObs_2/100SN_3/both_bad_seeing/',
                   'both_bad_seeing', fails=kst_fails2, wipe_fails=True)
bothg_fails2, bothg_diffs2 = analyse(u"\\\\?\\c:\\Users\\gltay\\AppData\\Local\\Packages\\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\\LocalState\\rootfs\\home\\gtaylor\\lightcurves\\Honours_data_sets\\5_270218\\1_stat_sample\\Kepler_6hours\\SM_5day\\vObs_2\\100SN_3\\both_good_seeing\\",
                   'both_good_seeing', fails=kst_fails2, wipe_fails=True)


# Plot different scope combinations

plot_diffs_2([smb_diffs2, smg_diffs2, kst_diffs2, bothb_diffs2, bothg_diffs2],
             ['SM bad seeing', 'SM good seeing', 'KST', 'Both bad seeing', 'Both good seeing'],
             ['b', 'r', 'g', 'y', 'k'],
             u"\\\\?\\c:\\Users\\gltay\\AppData\\Local\\Packages\\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\\LocalState\\rootfs\\home\\gtaylor\\lightcurves\\Honours_data_sets\\5_270218\\1_stat_sample\\Kepler_6hours\\SM_5day\\vObs_2\\100SN_3\\stats"
             u"\\errors_removed"
             u"\\all"
             u"\\")

plot_diffs_2([smb_diffs2],
             ['SM bad seeing'],
             ['b'],
             u"\\\\?\\c:\\Users\\gltay\\AppData\\Local\\Packages\\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\\LocalState\\rootfs\\home\\gtaylor\\lightcurves\\Honours_data_sets\\5_270218\\1_stat_sample\\Kepler_6hours\\SM_5day\\vObs_2\\100SN_3\\stats"
             u"\\errors_removed"
             u"\\SM_b\\")

plot_diffs_2([smg_diffs2],
             ['SM good seeing'],
             ['r'],
             u"\\\\?\\c:\\Users\\gltay\\AppData\\Local\\Packages\\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\\LocalState\\rootfs\\home\\gtaylor\\lightcurves\\Honours_data_sets\\5_270218\\1_stat_sample\\Kepler_6hours\\SM_5day\\vObs_2\\100SN_3\\stats"
             u"\\errors_removed"
             u"\\SM_g\\")

plot_diffs_2([kst_diffs2],
             ['KST'],
             ['g'],
             u"\\\\?\\c:\\Users\\gltay\\AppData\\Local\\Packages\\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\\LocalState\\rootfs\\home\\gtaylor\\lightcurves\\Honours_data_sets\\5_270218\\1_stat_sample\\Kepler_6hours\\SM_5day\\vObs_2\\100SN_3\\stats"
             u"\\errors_removed"
             u"\\KST"
             u"\\")

plot_diffs_2([bothb_diffs2],
             ['Both bad seeing'],
             ['y'],
             u"\\\\?\\c:\\Users\\gltay\\AppData\\Local\\Packages\\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\\LocalState\\rootfs\\home\\gtaylor\\lightcurves\\Honours_data_sets\\5_270218\\1_stat_sample\\Kepler_6hours\\SM_5day\\vObs_2\\100SN_3\\stats"
             u"\\errors_removed"
             u"\\Both_b\\")

plot_diffs_2([bothg_diffs2],
             ['Both good seeing'],
             ['k'],
             u"\\\\?\\c:\\Users\\gltay\\AppData\\Local\\Packages\\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\\LocalState\\rootfs\\home\\gtaylor\\lightcurves\\Honours_data_sets\\5_270218\\1_stat_sample\\Kepler_6hours\\SM_5day\\vObs_2\\100SN_3\\stats\\errors_removed\\Both_g\\")

plot_diffs_2([smb_diffs2, smg_diffs2],
             ['SM bad seeing', 'SM good seeing'],
             ['b', 'r'],
             u"\\\\?\\c:\\Users\\gltay\\AppData\\Local\\Packages\\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\\LocalState\\rootfs\\home\\gtaylor\\lightcurves\\Honours_data_sets\\5_270218\\1_stat_sample\\Kepler_6hours\\SM_5day\\vObs_2\\100SN_3\\stats\\errors_removed\\SM_g_SM_b\\")

plot_diffs_2([bothb_diffs2, bothg_diffs2],
             ['Both bad seeing', 'Both good seeing'],
             ['y', 'k'],
             u"\\\\?\\c:\\Users\\gltay\\AppData\\Local\\Packages\\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\\LocalState\\rootfs\\home\\gtaylor\\lightcurves\\Honours_data_sets\\5_270218\\1_stat_sample\\Kepler_6hours\\SM_5day\\vObs_2\\100SN_3\\stats\\errors_removed\\Both_g_Both_b\\")

plot_diffs_2([smg_diffs2, kst_diffs2, bothg_diffs2],
             ['SM good seeing', 'KST', 'Both good seeing'],
             ['r', 'g', 'k'],
             u"\\\\?\\c:\\Users\\gltay\\AppData\\Local\\Packages\\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\\LocalState\\rootfs\\home\\gtaylor\\lightcurves\\Honours_data_sets\\5_270218\\1_stat_sample\\Kepler_6hours\\SM_5day\\vObs_2\\100SN_3\\stats\\errors_removed\\SM_g_KST_Both_g\\"
            )

plot_diffs_2([smb_diffs2, kst_diffs2, bothb_diffs2],
             ['SM bad seeing', 'KST', 'Both bad seeing'],
             ['b', 'g', 'y'],
             u"\\\\?\\c:\\Users\\gltay\\AppData\\Local\\Packages\\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\\LocalState\\rootfs\\home\\gtaylor\\lightcurves\\Honours_data_sets\\5_270218\\1_stat_sample\\Kepler_6hours\\SM_5day\\vObs_2\\100SN_3\\stats\\errors_removed\\SM_b_KST_Both_b\\"
             )
