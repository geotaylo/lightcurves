import luminosity_distance_fitter as ldf
from prettytable import PrettyTable
import os
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pylab as py
import numpy as np
from scipy.optimize import curve_fit
import heaven_and_earth as hae


# Where the data comes from (simulation set)
parent = 'Honours_data_sets/062618/1000set/'


# Load parameters and SN IDs of failed fits, so only 'perfect' processes are included in diagram
smg_fails, smg_params = ldf.get_params(parent + 'sm_ws/',
                    'skymapper_ws')
smb_fails, smb_params = ldf.get_params(parent + 'sm_ps/',
                    'skymapper_ps')

kst_fails, kst_params = ldf.get_params(parent + 'kst/', 'kst')

# If the kepler fitting didn't work, the combined fitting won't be optimal as it won't have been passed Kepler's t0.
bothg_fails, bothg_params = ldf.get_params(parent + 'combined_ws/', 'combined_ws', fails=kst_fails)

bothb_fails, bothb_params = ldf.get_params(parent + 'combined_ps/', 'combined_ps', fails=kst_fails)


# Load true parameters, with which to compare to other fits
true_params = ldf.get_true(parent + 'kst/')

def mu_residuals(focus_params, true_params):
    broken_sn = []
    residual = []
    for i in range(len(focus_params[0])):
        c_f = focus_params[0][i]
        x0_f = focus_params[2][i]
        x1_f = focus_params[3][i]
        c_t = true_params[0][i]
        x0_t = true_params[2][i]
        x1_t = true_params[3][i]
        if x0_f <= 0.:
            broken_sn.append(i)
        else:
            mu_f = ldf.distance_modulus(x0_f, x1_f, c_f, parent + 'sm_ws/', i + 1)
            mu_t = ldf.distance_modulus(x0_t, x1_t, c_t, parent + 'sm_ws/', i + 1)
            x = mu_t - mu_f
            residual.append(x)
    return residual

def f(x, a, b, c):
    """
    Equation for a Gaussian
    """
    return a * py.exp(-(x - b) ** 2.0 / (2 * c ** 2))


def plot_residuals(sg, kst, bg, true, parent):
    sg_res = mu_residuals(sg, true)
    kst_res = mu_residuals(kst, true)
    bg_res = mu_residuals(bg, true)

    plt.hist(sg_res, range=[-10,10], bins=50, histtype='step',  color='mediumblue', label='SM_ws')
    plt.hist(kst_res, range=[-10,10], bins=50, histtype='step', color='g', label='KST')
    plt.hist(bg_res, range=[-10,10], bins=50, histtype='step',  color='crimson', label='Combined_ws')
    plt.legend()
    plt.savefig(parent+'cosmology/wide_res.png', dpi=200)
    plt.close()

    plt.hist(sg_res, range=[-2, 2], bins=20, histtype='step', color='mediumblue', label='SM_ws')
    plt.hist(kst_res, range=[-2, 2], bins=20, histtype='step', color='g', label='KST')
    plt.hist(bg_res, range=[-2, 2], bins=20, histtype='step', color='crimson', label='Combined_ws')
    plt.legend()
    plt.savefig(parent + 'cosmology/narrow_res.png', dpi=200)
    plt.close()

    plt.hist(sg_res, range=[-2, 2], bins=10, histtype='step', color='mediumblue', label='SM_ws')
    plt.hist(kst_res, range=[-2, 2], bins=10, histtype='step', color='g', label='KST')
    plt.hist(bg_res, range=[-2, 2], bins=10, histtype='step', color='crimson', label='Combined_ws')
    plt.legend()
    plt.savefig(parent + 'cosmology/narrow_res_lessbins.png', dpi=200)
    plt.close()

    figa, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True)
    data_sm = ax1.hist(sg_res, range=[-2, 2], bins=10, histtype='step', density=True, color='mediumblue', label='SM_ws')
    data_kst = ax1.hist(kst_res, range=[-2, 2], bins=40, histtype='step', density=True, color='g', label='KST')
    data_combo = ax1.hist(bg_res, range=[-2, 2], bins=10, histtype='step', density=True, color='crimson', label='Combined_ws')


    data = [data_sm, data_kst, data_combo]
    colour = ['mediumblue', 'g', 'crimson']
    labels = ['SM_ws', 'KST', 'Combined_ws']
    for i in range(len(data)):
        # Generate data from bins as a set of points
        x = [0.5 * (data[i][1][t] + data[i][1][t + 1]) for t in xrange(len(data[i][1]) - 1)]
        y = data[i][0]

        popt, pcov = curve_fit(f, x, y)
        print popt

        x_fit = py.linspace(-2, 2, 200)
        y_fit = f(x_fit, *popt)

        ax2.plot(x_fit, y_fit, lw=2, color=colour[i], label=labels[i]+'\n $\mu = %s$'%round(popt[1], 3)+'\n $\sigma = %s$'%round(popt[2], 3))

    figa.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in figa.axes[:-1]], visible=False)
    ax2.set_xlabel('Residual ($\mu_{true}- \mu_{fitted}$)')
    ax2.set_ylabel('PDF')
    leg = plt.legend(fontsize='large', loc='center left', bbox_to_anchor=(1, 1), ncol=1,
                     borderaxespad=0, frameon=False, labelspacing=1)
    for line in leg.get_lines():
        line.set_linewidth(4.0)
    figa.savefig(parent + 'cosmology/narrow_res_morebins_normed.png', dpi=200, bbox_extra_artists=(leg,), bbox_inches='tight')
    plt.close(figa)

    plt.hist(sg_res, range=[-0.5, 0.5], bins=20, histtype='step', color='mediumblue', label='SM_ws')
    plt.hist(kst_res, range=[-0.5, 0.5], bins=20, histtype='step', color='g', label='KST')
    plt.hist(bg_res, range=[-0.5, 0.5], bins=20, histtype='step', color='crimson', label='Combined_ws')
    plt.legend()
    plt.savefig(parent + 'cosmology/vnarrow_res.png', dpi=200)
    plt.close()

    plt.hist(sg_res, range=[-0.5, 0.5], bins=10, histtype='step', color='mediumblue', label='SM_ws')
    plt.hist(kst_res, range=[-0.5, 0.5], bins=10, histtype='step', color='g', label='KST')
    plt.hist(bg_res, range=[-0.5, 0.5], bins=10, histtype='step', color='crimson', label='Combined_ws')
    plt.legend()
    plt.savefig(parent + 'cosmology/vnarrow_res_lessbins.png', dpi=200)
    plt.close()

def plot_single(focus_params, set):
    distance_modulus = []
    luminosity_distance = []
    # indexes to delete when errors encountered
    del_i = []

    # Get distance for each SN
    for i in range(len(focus_params[0])):
        # Order of params: c, t0, x0, x1, z, SN
        c = focus_params[0][i]
        x0 = focus_params[2][i]
        x1 = focus_params[3][i]
        z = focus_params[4][i]

        # Can't take the log of a negative number - these should be counted as fails
        #  (need to go back through and ammend code at some point)
        if x0 <= 0.:
            del_i.append(i)

        else:
            mu = ldf.distance_modulus(x0, x1, c, parent+'sm_ws/', i+1)
            d_l = ldf.distance(mu)
            distance_modulus.append(mu)
            luminosity_distance.append(d_l)

    # Remove SN with negative x0 (flagged in previous loop)
    for i in sorted(del_i, reverse=True):
        del focus_params[0][i - 1]
        del focus_params[1][i - 1]
        del focus_params[2][i - 1]
        del focus_params[3][i - 1]
        del focus_params[4][i - 1]
        del focus_params[5][i - 1]

    # Creates residuals table
    t = PrettyTable()
    t.title = set
    t.add_column('SN', focus_params[5])
    t.add_column('c', focus_params[0])
    t.add_column('t0', focus_params[1])
    t.add_column('x0', focus_params[2])
    t.add_column('x1', focus_params[3])
    t.add_column('z', focus_params[4])
    t.add_column('mu', distance_modulus)
    t.add_column('distance(Mpc)', luminosity_distance)

    table_txt = t.get_string()

    # Writes results
    writefolder = parent + "cosmology/" + set + "/"
    if not os.path.isdir(writefolder):
        os.makedirs(writefolder)

    with open(writefolder + 'fitted_distances_' + set + '.txt', 'w') as file:
        file.write(table_txt)

    # speed of light (km/s)
    c = 299792

    # recessional velocity = z*c is this the right equation?
    vel = [i * c for i in focus_params[4]]

    plt.plot(focus_params[4], distance_modulus, 'o', markersize=2)
    plt.ylim(0, 50)
    plt.xlabel('Redshift')
    plt.ylabel(r'Distance Modulus, $\mu = -2.5\ln(x_0)-10.095 - M_b + (\alpha x_1)-(\beta c)$')
    plt.savefig(writefolder + 'z_vs_dm.png')
    plt.close()

    # plt.plot(focus_params[4], vel)
    # plt.show()
    # creates hubble diagram
    plt.plot(luminosity_distance, vel, 'o', markersize=2)
    plt.xlabel('Distance (Mpc)')
    plt.ylabel('Redshift velocity (km/s)')
    plt.xlim(0, 0.5)
    plt.savefig(writefolder + 'd_vs_v.png')
    return

# Which subset to focus on e.g SkyMapper poorly sampled (smb), Kepler (kst)
focus_params = kst_params

setname = 'kst'

def plot_multi(param_list, set_list):

    # fig, ax1 = plt.subplots(1)
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True, sharey=True)


    # fig2, (ax12, ax22, ax32, ax42) = plt.subplots(4, sharex=True, sharey=True)


    col = ['mediumblue', 'g', 'crimson', 'black']
    for i in range(len(param_list)):
        focus_params = param_list[i]
        set = set_list[i]
        distance_modulus = []
        luminosity_distance = []
        h_0 = []
        # indexes to delete when errors encountered
        del_i = []

        # Get distance for each SN
        for a in range(len(focus_params[0])):
            # Order of params: c, t0, x0, x1, z, SN
            c = focus_params[0][a]
            x0 = focus_params[2][a]
            x1 = focus_params[3][a]
            z = focus_params[4][a]

            # Can't take the log of a negative number - these should be counted as fails
            #  (need to go back through and ammend code at some point)
            if x0 <= 0.:
                del_i.append(a)

            else:
                mu = ldf.distance_modulus(x0, x1, c, parent + 'sm_ws/', a+1)
                d_l = ldf.distance(mu)
                hubest = ldf.h0_est(z, d_l)
                distance_modulus.append(mu)
                luminosity_distance.append(d_l)
                h_0.append(hubest)


        # Remove SN with negative x0 (flagged in previous loop)
        for q in sorted(del_i, reverse=True):
            del focus_params[0][q - 1]
            del focus_params[1][q - 1]
            del focus_params[2][q - 1]
            del focus_params[3][q - 1]
            del focus_params[4][q - 1]
            del focus_params[5][q - 1]


        # Writes results
        writefolder = parent + "cosmology/"
        if not os.path.isdir(writefolder):
            os.makedirs(writefolder)

        # speed of light (km/s)
        c = 299792

        # recessional velocity = z*c is this the right equation?
        vel = [g * c for g in focus_params[4]]

        axes = [ax1, ax2, ax3, ax4]
        # axes2 = [ax12, ax22, ax32, ax42]

        axes[i].plot(focus_params[4], distance_modulus, 'o', markersize=2, alpha=0.3, label = set, color=col[i])
        # leg = axes[i].legend(fontsize='large', loc='center left', bbox_to_anchor=(1, 0.5), ncol=1,
        #              borderaxespad=0, frameon=False, labelspacing=1)
        # leg.legendHandles[0]._sizes = [30]
        # ax1.plot(np.log10(focus_params[4]), distance_modulus, 'o', markersize=5, alpha=0.5, label = set, color=col[i])

        # axes2[i].hist(h_0, bins=50, range = [0,100], label=set, color=col[i])
        # axes2[i].legend()

    legend_elements = [Line2D([0], [0], marker='o', color='w', label=set_list[0],
                              markerfacecolor=col[0], markersize=7),
                       Line2D([0], [0], marker='o', color='w', label=set_list[1],
                              markerfacecolor=col[1], markersize=7),
                       Line2D([0], [0], marker='o', color='w', label=set_list[2],
                              markerfacecolor=col[2], markersize=7),
                       Line2D([0], [0], marker='o', color='w', label=set_list[3],
                              markerfacecolor=col[3], markersize=7)]

    leg = ax2.legend(handles=legend_elements, fontsize='large', loc='center left', bbox_to_anchor=(1, 0), ncol=1,
                                borderaxespad=0, frameon=False, labelspacing=1)

    ax1.set_ylim(27, 42)
    ax4.set_xlabel(r'$Z$')
    ax3.set_ylabel(r'     Distance Modulus, $\mu$')
    plt.yticks(np.arange(27, 47, 5.0))

    # plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)

    fig.savefig(writefolder + 'loghub.png', dpi=200, bbox_extra_artists=(leg,),
                 bbox_inches='tight')


#plot_multi([true_params], ['true'])

# plot_residuals(smg_params, kst_params, bothg_params, true_params, parent)

plot_multi([smg_params, kst_params, bothg_params, true_params], ['sm_ws', 'kst', 'combined_ws', 'true'])
