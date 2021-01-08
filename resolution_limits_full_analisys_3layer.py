from TDSA import *
import csv


def epsilon(e_s, e_inf, tau, freq):  # Debye model
    omg = 2 * pi * freq
    e_w = e_inf + (e_s - e_inf) / (1 + 1j * omg * tau)
    return e_w


def nk_from_eps(e_s, e_inf, tau, freq):
    e_w = epsilon(e_s, e_inf, tau, freq)
    n = sqrt((abs(e_w) + real(e_w)) / 2)
    k = sqrt((abs(e_w) - real(e_w)) / 2)
    return n, k


# fh = open('resolution_limit.csv', 'r')
# csvfile = csv.reader('resolution_limit.csv')
test_dirs = list()
for i in range(15):
    test_dirs.append('test_' + str(i+1))
error_dirs = ['1%', '2%', '5%', '10%']  # , '20%']

for test_dir in test_dirs:
    for error_dir in error_dirs:
        d_mat_i = list()
        d_mat_mean_i = list()
        d_mat_std_i = list()

        e_s_i = list()
        e_s_mean_i = list()
        e_s_std_i = list()

        e_inf_i = list()
        e_inf_mean_i = list()
        e_inf_std_i = list()

        tau_i = list()
        tau_mean_i = list()
        tau_std_i = list()

        d_mat_m = list()
        d_mat_mean_m = list()
        d_mat_std_m = list()

        e_s_m = list()
        e_s_mean_m = list()
        e_s_std_m = list()

        e_inf_m = list()
        e_inf_mean_m = list()
        e_inf_std_m = list()

        tau_m = list()
        tau_mean_m = list()
        tau_std_m = list()

        d_mat_o = list()
        d_mat_mean_o = list()
        d_mat_std_o = list()

        e_s_o = list()
        e_s_mean_o = list()
        e_s_std_o = list()

        e_inf_o = list()
        e_inf_mean_o = list()
        e_inf_std_o = list()

        tau_o = list()
        tau_mean_o = list()
        tau_std_o = list()

        d_air = list()
        d_air_mean = list()
        d_air_std = list()

        f_cutoff = list()
        pDr = list()

        sample_dir = test_dir + '/' + error_dir

        print('Processing', sample_dir)

        with open('./output/simulation_results_3_layer/' + sample_dir + '/resolution_limit.csv', 'r') as f:
            reader = csv.reader(f)
            for row in reader:

                d_mat_i.append(float(row[0]))
                d_mat_mean_i.append(float(row[1]))
                d_mat_std_i.append(float(row[2]))
                e_s_i.append(float(row[3]))
                e_s_mean_i.append(float(row[4]))
                e_s_std_i.append(float(row[5]))
                e_inf_i.append(float(row[6]))
                e_inf_mean_i.append(float(row[7]))
                e_inf_std_i.append(float(row[8]))
                tau_i.append(float(row[9]))
                tau_mean_i.append(float(row[10]))
                tau_std_i.append(float(row[11]))
                d_mat_m.append(float(row[12]))
                d_mat_mean_m.append(float(row[13]))
                d_mat_std_m.append(float(row[14]))
                e_s_m.append(float(row[15]))
                e_s_mean_m.append(float(row[16]))
                e_s_std_m.append(float(row[17]))
                e_inf_m.append(float(row[18]))
                e_inf_mean_m.append(float(row[19]))
                e_inf_std_m.append(float(row[20]))
                tau_m.append(float(row[21]))
                tau_mean_m.append(float(row[22]))
                tau_std_m.append(float(row[23]))
                d_mat_o.append(float(row[24]))
                d_mat_mean_o.append(float(row[25]))
                d_mat_std_o.append(float(row[26]))
                e_s_o.append(float(row[27]))
                e_s_mean_o.append(float(row[28]))
                e_s_std_o.append(float(row[29]))
                e_inf_o.append(float(row[30]))
                e_inf_mean_o.append(float(row[31]))
                e_inf_std_o.append(float(row[32]))
                tau_o.append(float(row[33]))
                tau_mean_o.append(float(row[34]))
                tau_std_o.append(float(row[35]))
                d_air.append(float(row[36]))
                d_air_mean.append(float(row[37]))
                d_air_std.append(float(row[38]))
                f_cutoff.append(float(row[39]))
                pDr.append(float(row[40]))

        d_mat_i = array(d_mat_i)  # * 1e6
        d_mat_mean_i = array(d_mat_mean_i)  # * 1e6
        d_mat_std_i = array(d_mat_std_i)  # * 1e6
        d_mat_pDr_i = array(pDr)
        e_s_i = array(e_s_i)
        e_s_mean_i = array(e_s_mean_i)
        e_s_std_i = array(e_s_std_i)
        e_inf_i = array(e_inf_i)
        e_inf_mean_i = array(e_inf_mean_i)
        e_inf_std_i = array(e_inf_std_i)
        tau_i = array(tau_i)
        tau_mean_i = array(tau_mean_i)
        tau_std_i = array(tau_std_i)
        d_mat_m = array(d_mat_m)  # * 1e6
        d_mat_mean_m = array(d_mat_mean_m)  # * 1e6
        d_mat_std_m = array(d_mat_std_m)  # * 1e6
        d_mat_pDr_m = array(pDr)
        e_s_m = array(e_s_m)
        e_s_mean_m = array(e_s_mean_m)
        e_s_std_m = array(e_s_std_m)
        e_inf_m = array(e_inf_m)
        e_inf_mean_m = array(e_inf_mean_m)
        e_inf_std_m = array(e_inf_std_m)
        tau_m = array(tau_m)
        tau_mean_m = array(tau_mean_m)
        tau_std_m = array(tau_std_m)
        d_mat_o = array(d_mat_o)  # * 1e6
        d_mat_mean_o = array(d_mat_mean_o)  # * 1e6
        d_mat_std_o = array(d_mat_std_o)  # * 1e6
        d_mat_pDr_o = array(pDr)
        e_s_o = array(e_s_o)
        e_s_mean_o = array(e_s_mean_o)
        e_s_std_o = array(e_s_std_o)
        e_inf_o = array(e_inf_o)
        e_inf_mean_o = array(e_inf_mean_o)
        e_inf_std_o = array(e_inf_std_o)
        tau_o = array(tau_o)
        tau_mean_o = array(tau_mean_o)
        tau_std_o = array(tau_std_o)

        d_air = array(d_air)
        d_air_mean = array(d_air_mean)
        d_air_std = array(d_air_std)
        f_cutoff = array(f_cutoff)
        pDr = array(pDr)
        lmda = (c_0 / (f_cutoff * 1e12)) * 1e6  # um

        d_mat_t = (d_mat_i + d_mat_m + d_mat_o) / 3
        d_mat_mean_t = (d_mat_mean_i + d_mat_mean_m + d_mat_mean_o) / 3
        d_mat_std_t = sqrt(d_mat_std_i**2 + d_mat_std_m**2 + d_mat_std_o**2) / 3

        figure(error_dir)
        ax = axes()
        ax.set_xscale('log')
        ax.set_yscale('log')

        # if test_dir in ['test_1', 'test_2', 'test_3']:
        #     if test_dir == 'test_1':
        #         ax.plot(d_mat_i, abs(d_mat_mean_i - d_mat_i) / d_mat_i, 'b.')
        #         ax.plot(d_mat_m, abs(d_mat_mean_m - d_mat_m) / d_mat_m, 'b.')
        #         ax.plot(d_mat_o, abs(d_mat_mean_o - d_mat_o) / d_mat_o, 'b.')
        #     elif test_dir == 'test_2':
        #         ax.plot(d_mat_i, abs(d_mat_mean_i - d_mat_i) / d_mat_i, 'b+')
        #         ax.plot(d_mat_m, abs(d_mat_mean_m - d_mat_m) / d_mat_m, 'b+')
        #         ax.plot(d_mat_o, abs(d_mat_mean_o - d_mat_o) / d_mat_o, 'b+')
        #     elif test_dir == 'test_3':
        #         ax.plot(d_mat_i, abs(d_mat_mean_i - d_mat_i) / d_mat_i, 'bx')
        #         ax.plot(d_mat_m, abs(d_mat_mean_m - d_mat_m) / d_mat_m, 'bx')
        #         ax.plot(d_mat_o, abs(d_mat_mean_o - d_mat_o) / d_mat_o, 'bx')
        # elif test_dir in ['test_4', 'test_5', 'test_6']:
        #     if test_dir == 'test_4':
        #         ax.plot(d_mat_i, abs(d_mat_mean_i - d_mat_i) / d_mat_i, 'g.')
        #         ax.plot(d_mat_m, abs(d_mat_mean_m - d_mat_m) / d_mat_m, 'g.')
        #         ax.plot(d_mat_o, abs(d_mat_mean_o - d_mat_o) / d_mat_o, 'g.')
        #     elif test_dir == 'test_5':
        #         ax.plot(d_mat_i, abs(d_mat_mean_i - d_mat_i) / d_mat_i, 'g+')
        #         ax.plot(d_mat_m, abs(d_mat_mean_m - d_mat_m) / d_mat_m, 'g+')
        #         ax.plot(d_mat_o, abs(d_mat_mean_o - d_mat_o) / d_mat_o, 'g+')
        #     elif test_dir == 'test_6':
        #         ax.plot(d_mat_i, abs(d_mat_mean_i - d_mat_i) / d_mat_i, 'gx')
        #         ax.plot(d_mat_m, abs(d_mat_mean_m - d_mat_m) / d_mat_m, 'gx')
        #         ax.plot(d_mat_o, abs(d_mat_mean_o - d_mat_o) / d_mat_o, 'gx')
        # elif test_dir in ['test_13', 'test_14', 'test_15']:
        #     if test_dir == 'test_13':
        #         ax.plot(d_mat_i, abs(d_mat_mean_i - d_mat_i) / d_mat_i, 'r.')
        #         ax.plot(d_mat_m, abs(d_mat_mean_m - d_mat_m) / d_mat_m, 'r.')
        #         ax.plot(d_mat_o, abs(d_mat_mean_o - d_mat_o) / d_mat_o, 'r.')
        #     elif test_dir == 'test_14':
        #         ax.plot(d_mat_i, abs(d_mat_mean_i - d_mat_i) / d_mat_i, 'r+')
        #         ax.plot(d_mat_m, abs(d_mat_mean_m - d_mat_m) / d_mat_m, 'r+')
        #         ax.plot(d_mat_o, abs(d_mat_mean_o - d_mat_o) / d_mat_o, 'r+')
        #     elif test_dir == 'test_15':
        #         ax.plot(d_mat_i, abs(d_mat_mean_i - d_mat_i) / d_mat_i, 'rx')
        #         ax.plot(d_mat_m, abs(d_mat_mean_m - d_mat_m) / d_mat_m, 'rx')
        #         ax.plot(d_mat_o, abs(d_mat_mean_o - d_mat_o) / d_mat_o, 'rx')
        # else:
        #     if test_dir in ['test_7', 'test_10']:
        #         ax.plot(d_mat_i, abs(d_mat_mean_i - d_mat_i) / d_mat_i, 'y.')
        #         ax.plot(d_mat_m, abs(d_mat_mean_m - d_mat_m) / d_mat_m, 'y.')
        #         ax.plot(d_mat_o, abs(d_mat_mean_o - d_mat_o) / d_mat_o, 'y.')
        #     elif test_dir in ['test_8', 'test_11']:
        #         ax.plot(d_mat_i, abs(d_mat_mean_i - d_mat_i) / d_mat_i, 'y+')
        #         ax.plot(d_mat_m, abs(d_mat_mean_m - d_mat_m) / d_mat_m, 'y+')
        #         ax.plot(d_mat_o, abs(d_mat_mean_o - d_mat_o) / d_mat_o, 'y+')
        #     elif test_dir in ['test_9', 'test_12']:
        #         ax.plot(d_mat_i, abs(d_mat_mean_i - d_mat_i) / d_mat_i, 'yx')
        #         ax.plot(d_mat_m, abs(d_mat_mean_m - d_mat_m) / d_mat_m, 'yx')
        #         ax.plot(d_mat_o, abs(d_mat_mean_o - d_mat_o) / d_mat_o, 'yx')
        # if test_dir in ['test_1', 'test_2', 'test_3']:
        #     if test_dir == 'test_1':
        #         ax.plot(d_mat_t, abs(d_mat_mean_t - d_mat_t) / d_mat_t, 'b.')
        #     elif test_dir == 'test_2':
        #         ax.plot(d_mat_t, abs(d_mat_mean_t - d_mat_t) / d_mat_t, 'b+')
        #     elif test_dir == 'test_3':
        #         ax.plot(d_mat_t, abs(d_mat_mean_t - d_mat_t) / d_mat_t, 'bx')
        # elif test_dir in ['test_4', 'test_5', 'test_6']:
        #     if test_dir == 'test_4':
        #         ax.plot(d_mat_t, abs(d_mat_mean_t - d_mat_t) / d_mat_t, 'g.')
        #     elif test_dir == 'test_5':
        #         ax.plot(d_mat_t, abs(d_mat_mean_t - d_mat_t) / d_mat_t, 'g+')
        #     elif test_dir == 'test_6':
        #         ax.plot(d_mat_t, abs(d_mat_mean_t - d_mat_t) / d_mat_t, 'gx')
        # elif test_dir in ['test_13', 'test_14', 'test_15']:
        #     if test_dir == 'test_13':
        #         ax.plot(d_mat_t, abs(d_mat_mean_t - d_mat_t) / d_mat_t, 'r.')
        #     elif test_dir == 'test_14':
        #         ax.plot(d_mat_t, abs(d_mat_mean_t - d_mat_t) / d_mat_t, 'r+')
        #     elif test_dir == 'test_15':
        #         ax.plot(d_mat_t, abs(d_mat_mean_t - d_mat_t) / d_mat_t, 'rx')
        # else:
        #     if test_dir in ['test_7', 'test_10']:
        #         ax.plot(d_mat_t, abs(d_mat_mean_t - d_mat_t) / d_mat_t, 'y.')
        #     elif test_dir in ['test_8', 'test_11']:
        #         ax.plot(d_mat_t, abs(d_mat_mean_t - d_mat_t) / d_mat_t, 'y+')
        #     elif test_dir in ['test_9', 'test_12']:
        #         ax.plot(d_mat_t, abs(d_mat_mean_t - d_mat_t) / d_mat_t, 'yx')
        if test_dir in ['test_1', 'test_2', 'test_3']:
            if test_dir == 'test_1':
                ax.plot(d_mat_t, d_mat_mean_t, 'b.')
            elif test_dir == 'test_2':
                ax.plot(d_mat_t, d_mat_mean_t, 'b+')
            elif test_dir == 'test_3':
                ax.plot(d_mat_t, d_mat_mean_t, 'bx')
        elif test_dir in ['test_4', 'test_5', 'test_6']:
            if test_dir == 'test_4':
                ax.plot(d_mat_t, d_mat_mean_t, 'g.')
            elif test_dir == 'test_5':
                ax.plot(d_mat_t, d_mat_mean_t, 'g+')
            elif test_dir == 'test_6':
                ax.plot(d_mat_t, d_mat_mean_t, 'gx')
        elif test_dir in ['test_13', 'test_14', 'test_15']:
            if test_dir == 'test_13':
                ax.plot(d_mat_t, d_mat_mean_t, 'r.')
            elif test_dir == 'test_14':
                ax.plot(d_mat_t, d_mat_mean_t, 'r+')
            elif test_dir == 'test_15':
                ax.plot(d_mat_t, d_mat_mean_t, 'rx')
        else:
            if test_dir in ['test_7', 'test_10']:
                ax.plot(d_mat_t, d_mat_mean_t, 'y.')
            elif test_dir in ['test_8', 'test_11']:
                ax.plot(d_mat_t, d_mat_mean_t, 'y+')
            elif test_dir in ['test_9', 'test_12']:
                ax.plot(d_mat_t, d_mat_mean_t, 'yx')

        # if test_dir in ['test_1', 'test_2', 'test_3']:
        #     if test_dir == 'test_1':
        #         ax.errorbar(d_mat_t, abs(d_mat_mean_t - d_mat_t) / d_mat_t, yerr=d_mat_std_t / d_mat_t, fmt='b.')
        #     elif test_dir == 'test_2':
        #         ax.errorbar(d_mat_t, abs(d_mat_mean_t - d_mat_t) / d_mat_t, yerr=d_mat_std_t / d_mat_t, fmt='b+')
        #     elif test_dir == 'test_3':
        #         ax.errorbar(d_mat_t, abs(d_mat_mean_t - d_mat_t) / d_mat_t, yerr=d_mat_std_t / d_mat_t, fmt='bx')
        # elif test_dir in ['test_4', 'test_5', 'test_6']:
        #     if test_dir == 'test_4':
        #         ax.errorbar(d_mat_t, abs(d_mat_mean_t - d_mat_t) / d_mat_t, yerr=d_mat_std_t / d_mat_t, fmt='g.')
        #     elif test_dir == 'test_5':
        #         ax.errorbar(d_mat_t, abs(d_mat_mean_t - d_mat_t) / d_mat_t, yerr=d_mat_std_t / d_mat_t, fmt='g+')
        #     elif test_dir == 'test_6':
        #         ax.errorbar(d_mat_t, abs(d_mat_mean_t - d_mat_t) / d_mat_t, yerr=d_mat_std_t / d_mat_t, fmt='gx')
        # elif test_dir in ['test_13',   'test_14', 'test_15']:
        #     if test_dir == 'test_13':
        #         ax.errorbar(d_mat_t, abs(d_mat_mean_t - d_mat_t) / d_mat_t, yerr=d_mat_std_t / d_mat_t, fmt='r.')
        #     elif test_dir == 'test_14':
        #         ax.errorbar(d_mat_t, abs(d_mat_mean_t - d_mat_t) / d_mat_t, yerr=d_mat_std_t / d_mat_t, fmt='r+')
        #     elif test_dir == 'test_15':
        #         ax.errorbar(d_mat_t, abs(d_mat_mean_t - d_mat_t) / d_mat_t, yerr=d_mat_std_t / d_mat_t, fmt='rx')
        # else:
        #     if test_dir in ['test_7', 'test_10']:
        #         ax.errorbar(d_mat_t, abs(d_mat_mean_t - d_mat_t) / d_mat_t, yerr=d_mat_std_t / d_mat_t, fmt='y.')
        #     elif test_dir in ['test_8', 'test_11']:
        #         ax.errorbar(d_mat_t, abs(d_mat_mean_t - d_mat_t) / d_mat_t, yerr=d_mat_std_t / d_mat_t, fmt='y+')
        #     elif test_dir in ['test_9', 'test_12']:
        #         ax.errorbar(d_mat_t, abs(d_mat_mean_t - d_mat_t) / d_mat_t, yerr=d_mat_std_t / d_mat_t, fmt='yx')

        legend_disp = [
            Line2D([0], [0], color='b', marker='.', label='No disp', linestyle='None'),
            Line2D([0], [0], color='g', marker='.', label='Low disp', linestyle='None'),
            Line2D([0], [0], color='y', marker='.', label='Mid disp', linestyle='None'),
            Line2D([0], [0], color='r', marker='.', label='High disp', linestyle='None')
        ]
        legend_contr = [
            Line2D([0], [0], color='b', marker='.', label='High contr', linestyle='None'),
            Line2D([0], [0], color='b', marker='+', label='Mid contr', linestyle='None'),
            Line2D([0], [0], color='b', marker='x', label='Low contr', linestyle='None')
        ]
        legend1 = legend(handles=legend_disp, loc=1)
        legend2 = legend(handles=legend_contr, loc=3)
        ax.add_artist(legend1)
        ax.add_artist(legend2)
        xlabel(r'$d_{sim}$')
        ylabel('Dispersion')
        # ax.plot(d_mat_t, ones(d_mat_t.size) / 100, 'k-', lw=0.5)
        ax.plot(d_mat_t, d_mat_t, 'k-', lw=0.5)
        ylim([1e-5, 2e4])

        if test_dir == 'test_15':
            savefig('./output/simulation_results_3_layer/' + error_dir.replace('%','_100_3') + '.pdf')
# show()
