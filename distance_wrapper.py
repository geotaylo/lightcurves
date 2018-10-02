import luminosity_distance_fitter as ldf
from prettytable import PrettyTable
import os
import matplotlib.pyplot as plt
import run_sncosmo_k2fields as hae


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
            mu = ldf.distance_modulus(x0, x1, c)
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
    fig, ax1 = plt.subplots(1)
    for i in range(len(param_list)):
        focus_params = param_list[i]
        set = set_list[i]
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
                mu = ldf.distance_modulus(x0, x1, c)
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


        # Writes results
        writefolder = parent + "cosmology/"
        if not os.path.isdir(writefolder):
            os.makedirs(writefolder)

        # speed of light (km/s)
        c = 299792

        # recessional velocity = z*c is this the right equation?
        vel = [i * c for i in focus_params[4]]

        ax1.plot(focus_params[4], distance_modulus, 'o', markersize=1, label = set)

    ax1.set_ylim(0, 30)
    ax1.set_xlabel('Redshift')
    ax1.set_ylabel(r'Distance Modulus, $\mu$')
    ax1.legend()

    plt.savefig(writefolder + 'all.png')

plot_multi([smg_params, smb_params, kst_params, bothg_params, bothb_params, true_params], ['sm_ws', 'sm_ps', 'kst', 'combined_ws', 'combined_ps', 'true'])

print ldf.distance_modulus(0.00082198, 2.0, 0.3)
