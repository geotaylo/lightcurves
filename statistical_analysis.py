""" Run a statistical anaylsis of fitted vs true parameters from model lightcurves of type 1a SN, fitted with sncosmo."""
import csv
from prettytable import PrettyTable
import os
import matplotlib.pyplot as plt
import numpy as np


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









    # Calculate differences
    diff_c = []
    diff_t0 = []
    diff_x0 = []
    diff_x1 = []
    diff_z = []


    for i in range(len(true_c)):
        diff_c.append(true_c[i] - fitted_c[i])
        diff_t0.append(true_t0[i] - fitted_t0[i])
        diff_x0.append(true_x0[i] - fitted_x0[i])
        diff_x1.append(true_x1[i] - fitted_x1[i])
        diff_z.append(true_z[i] - fitted_z[i])



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
    t.add_column('c-diff', diff_c)
    t.add_column('t0-diff', diff_t0)
    t.add_column('x0-diff', diff_x0)
    t.add_column('x1-diff', diff_x1)
    t.add_column('z-diff', diff_z)

    print t

    table_txt = t.get_string()

    writefolder = folder + 'stats/'
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
    abs_list = map(abs, list)
    weights = np.ones_like(abs_list) / float(len(list))
    return abs_list, weights

def plot_diffs(smb_diffs, smg_diffs, kst_diffs, bothb_diffs, bothg_diffs):

    abs_c_smb, weights_c_smb = abs_and_weight(smb_diffs[0])
    abs_c_smg, weights_c_smg = abs_and_weight(smg_diffs[0])
    abs_c_kst, weights_c_kst = abs_and_weight(kst_diffs[0])
    abs_c_bothb, weights_c_bothb = abs_and_weight(bothb_diffs[0])
    abs_c_bothg, weights_c_bothg = abs_and_weight(bothg_diffs[0])

    # Plot all telescopes on same graph
    plt.hist(abs_c_smb, bins=50, range =[0,1], histtype='step', weights=weights_c_smb, color='b', label='SM bad seeing')
    plt.hist(abs_c_smg, bins=50, range =[0,1], histtype='step', weights=weights_c_smg, color='r', label='SM good seeing')
    plt.hist(abs_c_kst, bins=50, range =[0,1], histtype='step', weights=weights_c_kst, color='g', label='KST')
    plt.hist(abs_c_bothb, bins=50, range =[0,1], histtype='step', weights=weights_c_bothb, color='y', label='Both bad seeing')
    plt.hist(abs_c_bothg, bins=50, range =[0,1], histtype='step', weights=weights_c_bothg, color='k', label='Both good seeing')
    plt.title('Colour Residuals')
    plt.xlabel('Residual (true - fitted value of c)')
    plt.ylabel('Probability')
    plt.legend()
    plt.show()


smb_fails, smb_diffs = analyse('Honours_data_sets/5_270218/1_stat_sample/Kepler_6hours/SM_5day/vObs_2/100SN_3/sm_bad_seeing/',
                    'sm_bad_seeing')
smg_fails, smg_diffs = analyse('Honours_data_sets/5_270218/1_stat_sample/Kepler_6hours/SM_5day/vObs_2/100SN_3/sm_good_seeing/',
                    'sm_good_seeing')
kst_fails, kst_diffs = analyse('Honours_data_sets/5_270218/1_stat_sample/Kepler_6hours/SM_5day/vObs_2/100SN_3/kst/', 'kst')
bothb_fails, bothb_diffs = analyse('Honours_data_sets/5_270218/1_stat_sample/Kepler_6hours/SM_5day/vObs_2/100SN_3/both_bad_seeing/',
                   'both_bad_seeing', fails=kst_fails)
bothg_fails, bothg_diffs = analyse('Honours_data_sets/5_270218/1_stat_sample/Kepler_6hours/SM_5day/vObs_2/100SN_3/both_good_seeing/',
                   'both_good_seeing', fails=kst_fails)

plot_diffs(smb_diffs, smg_diffs, kst_diffs, bothb_diffs, bothg_diffs)





smb_fails2, smb_diffs2 = analyse('Honours_data_sets/5_270218/1_stat_sample/Kepler_6hours/SM_5day/vObs_2/100SN_3/sm_bad_seeing/',
                    'sm_bad_seeing', wipe_fails=True)
smg_fails2, smg_diffs2 = analyse('Honours_data_sets/5_270218/1_stat_sample/Kepler_6hours/SM_5day/vObs_2/100SN_3/sm_good_seeing/',
                    'sm_good_seeing', wipe_fails=True)
kst_fails2, kst_diffs2 = analyse('Honours_data_sets/5_270218/1_stat_sample/Kepler_6hours/SM_5day/vObs_2/100SN_3/kst/', 'kst',
                   wipe_fails=True)
bothb_fails2, bothb_diffs2 = analyse('Honours_data_sets/5_270218/1_stat_sample/Kepler_6hours/SM_5day/vObs_2/100SN_3/both_bad_seeing/',
                   'both_bad_seeing', fails=kst_fails2, wipe_fails=True)
bothg_fails2, bothg_diffs2 = analyse('Honours_data_sets/5_270218/1_stat_sample/Kepler_6hours/SM_5day/vObs_2/100SN_3/both_good_seeing/',
                   'both_good_seeing', fails=kst_fails2, wipe_fails=True)
plot_diffs(smb_diffs2, smg_diffs2, kst_diffs2, bothb_diffs2, bothg_diffs2)