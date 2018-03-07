""" Run a statistical anaylsis of fitted vs true parameters from model lightcurves of type 1a SN, fitted with sncosmo."""
import csv
from prettytable import PrettyTable
from prettytable import MSWORD_FRIENDLY


def analyse(folder, set):

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
        sn_num[i-1] = 'error'+ str(i)



    t = PrettyTable()
    t.title = set
    t.add_column('SN', sn_num)
    t.add_column('c-diff', diff_c)
    t.add_column('t0-diff', diff_t0)
    t.add_column('x0-diff', diff_x0)
    t.add_column('x1-diff', diff_x1)
    t.add_column('z-diff', diff_z)

    t.set_style(MSWORD_FRIENDLY)
    print t

    table_txt = t.get_string()
    with open(folder + 'output.txt', 'w') as file:
        file.write(table_txt)

    return

analyse('Honours_data_sets/5_270218/1_stat_sample/Kepler_6hours/SM_5day/vObs_2/100SN_2/sm_bad_seeing/', 'sm_bad_seeing')
analyse('Honours_data_sets/5_270218/1_stat_sample/Kepler_6hours/SM_5day/vObs_2/100SN_2/kst/', 'kst')
analyse('Honours_data_sets/5_270218/1_stat_sample/Kepler_6hours/SM_5day/vObs_2/100SN_2/both_bad_seeing/', 'both_bad_seeing')