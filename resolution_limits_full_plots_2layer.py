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

        with open('./output/simulation_results/' + sample_dir + '/resolution_limit.csv', 'r') as f:
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
                d_mat_o.append(float(row[12]))
                d_mat_mean_o.append(float(row[13]))
                d_mat_std_o.append(float(row[14]))
                e_s_o.append(float(row[15]))
                e_s_mean_o.append(float(row[16]))
                e_s_std_o.append(float(row[17]))
                e_inf_o.append(float(row[18]))
                e_inf_mean_o.append(float(row[19]))
                e_inf_std_o.append(float(row[20]))
                tau_o.append(float(row[21]))
                tau_mean_o.append(float(row[22]))
                tau_std_o.append(float(row[23]))
                d_air.append(float(row[24]))
                d_air_mean.append(float(row[25]))
                d_air_std.append(float(row[26]))
                f_cutoff.append(float(row[27]))
                pDr.append(float(row[28]))


        # print(pDr)
        # quit()
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

        figure(test_dir + '1')
        ax = axes()
        ax.set_xscale('log')
        ax.set_yscale('log')
        # if error_dir == '1%':
        #     ax.plot(d_mat_i, d_mat_std_i / d_mat_mean_i, '.b')
        # elif error_dir == '2%':
        #     ax.plot(d_mat_i, d_mat_std_i / d_mat_mean_i, '.k')
        # elif error_dir == '5%':
        #     ax.plot(d_mat_i, d_mat_std_i / d_mat_mean_i, '.g')
        # elif error_dir == '10%':
        #     ax.plot(d_mat_i, d_mat_std_i / d_mat_mean_i, '.y')
        # elif error_dir == '20%':
        #     ax.plot(d_mat_i, d_mat_std_i / d_mat_mean_i, '.r')
        if error_dir == '1%':
            ax.errorbar(d_mat_i, d_mat_mean_i, yerr=d_mat_std_i, fmt='.b')
        elif error_dir == '2%':
            ax.errorbar(d_mat_i, d_mat_mean_i, yerr=d_mat_std_i, fmt='.c')
        elif error_dir == '5%':
            ax.errorbar(d_mat_i, d_mat_mean_i, yerr=d_mat_std_i, fmt='.g')
        elif error_dir == '10%':
            ax.errorbar(d_mat_i, d_mat_mean_i, yerr=d_mat_std_i, fmt='.r')
        elif error_dir == '20%':
            ax.errorbar(d_mat_i, d_mat_mean_i, yerr=d_mat_std_i, fmt='.y')

        figure(test_dir + '2')
        ax = axes()
        ax.set_xscale('log')
        ax.set_yscale('log')
        if error_dir == '1%':
            ax.errorbar(d_mat_o, d_mat_mean_o, yerr=d_mat_std_o, fmt='.b')
        elif error_dir == '2%':
            ax.errorbar(d_mat_o, d_mat_mean_o, yerr=d_mat_std_o, fmt='.c')
        elif error_dir == '5%':
            ax.errorbar(d_mat_o, d_mat_mean_o, yerr=d_mat_std_o, fmt='.g')
        elif error_dir == '10%':
            ax.errorbar(d_mat_o, d_mat_mean_o, yerr=d_mat_std_o, fmt='.r')
        elif error_dir == '20%':
            ax.errorbar(d_mat_o, d_mat_mean_o, yerr=d_mat_std_o, fmt='.y')
        # plot(d_mat_m, d_mat_std_m / d_mat_mean_m, '.r')
        # plot(d_mat_o, d_mat_std_o / d_mat_mean_o, '.g')

    legend_elements = [
        Line2D([0], [0], color='b', marker ='.'),
        Line2D([0], [0], color='c', marker='.'),
        Line2D([0], [0], color='g', marker='.'),
        Line2D([0], [0], color='r', marker='.'),
        Line2D([0], [0], color='y', marker='.')
    ]

    figure(test_dir + '1')
    plot(d_mat_i, d_mat_i, '-', c='black', label='expected')
    # legend(handles=legend_elements)
    savefig('./output/simulation_results/' + test_dir + '/' + test_dir + '_i.png')
    figure(test_dir + '3')
    plot(d_mat_o, d_mat_o, '-', c='black', label='expected')
    # legend(handles=legend_elements)
    savefig('./output/simulation_results/' + test_dir + '/' + test_dir + '_o.png')

# show()
