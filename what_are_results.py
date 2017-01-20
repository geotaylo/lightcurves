""" Run code to get histograms of colour and stretch discrepancies between
real and fitted values, for sncosmo sims."""
import matplotlib.pyplot as plt
import csv

parent = 'cad_sm_4d/cad_k_6h/30sn/'

folder_1 = 'cad_sm_4d/cad_k_6h/30sn/sm_bad_seeing/'
name_1 = 'SkyMapper (bad seeing) 30SN sample - '

folder_2 = 'cad_sm_4d/cad_k_6h/30sn/sm_good_seeing/'
name_2 = 'SkyMapper (good seeing) 30SN sample - '

folder_3 = 'cad_sm_4d/cad_k_6h/30sn/k_bad_seeing/'
name_3 = 'Kepler 30SN sample - '

folder_4 = 'cad_sm_4d/cad_k_6h/30sn/both_bad_seeing/'
name_4 = 'SkyMapper (bad seeing) and Kepler 30SN sample - '

folder_5 = 'cad_sm_4d/cad_k_6h/30sn/both_good_seeing/'
name_5 = 'SkyMapper (good seeing) and Kepler 30SN sample - '

def get_diff_util(folder, name):

    """ Reads text files of fitted and true parameters, and calculates
        diffences for each SN parameter """

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

            fitted_c.append(float(row[1].strip('c:')))
            fitted_t0.append(float(row[2].strip('t0:')))
            fitted_x0.append(float(row[3].strip('x0:')))
            fitted_x1.append(float(row[4].strip('x1:')))
            fitted_z.append(float(row[5].strip('z:')))

    true_c = []
    true_t0 = []
    true_x0 = []
    true_x1 = []
    true_z = []

    with open(tp, 'rb') as file:

        reader = csv.reader(file, delimiter=' ')

        for row in reader:

            true_c.append(float(row[1].strip('c:')))
            true_t0.append(float(row[2].strip('t0:')))
            true_x0.append(float(row[3].strip('x0:')))
            true_x1.append(float(row[4].strip('x1:')))
            true_z.append(float(row[5].strip('z:')))


    c_diff = []
    x0_diff = []

    c_string = 'Excluded: \n'
    x_string = 'Excluded: \n'
    for i in range(len(true_c)):
        if abs(true_c[i] - fitted_c[i]) > 1:
            print 'Colour value %s excluded: absolute difference exceeds 1'%(true_c[i] - fitted_c[i])
            c_string = c_string + str(true_c[i] - fitted_c[i]) + '\n'
        if abs(true_x0[i] - fitted_x0[i]) > 1:
            print 'Stretch value %s excluded: absolute difference exceeds ' \
                  '1'%abs(true_x0[i] - fitted_x0[i])
            x_string = x_string + str(true_x0[i] - fitted_x0[i]) + '\n'
        c_diff.append(true_c[i] - fitted_c[i])
        x0_diff.append(true_x0[i] - fitted_x0[i])


    plt.hist(c_diff, 15, range=(-0.5,0.5))
    plt.title(name + 'colour')
    plt.xlabel('True - Fitted')
    plt.ylabel('Frequency')
    plt.figtext(0.2, 0.7, c_string)
    plt.savefig(parent + name + 'colour.png')
    plt.show()

    plt.hist(x0_diff, 15, range=(-0.2,0.2))
    plt.title(name + 'stretch')
    plt.xlabel('True - Fitted')
    plt.ylabel('Frequency')
    plt.figtext(0.2, 0.7, x_string)
    plt.savefig(parent + name + 'stretch.png')
    plt.show()

    return

get_diff_util(folder_1, name_1)
get_diff_util(folder_2, name_2)
get_diff_util(folder_3, name_3)
get_diff_util(folder_4, name_4)
get_diff_util(folder_5, name_5)