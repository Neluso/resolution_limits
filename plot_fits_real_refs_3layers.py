from TDSA import *
import scipy.interpolate as intplt
from scipy.optimize import differential_evolution
from matplotlib import gridspec


deg_in = 0  # incidence angle in degrees
snell_sin = n_air * sin(deg_in * pi / 180)
# n_subs = 1.17 - 0.0 * 1j  # substrate refractive index -- cork
n_subs = 1e20 - 0.0 * 1j  # substrate refractive index -- metal
eps = 1e-20


# function definitions
def theta(n):
    return arcsin(snell_sin / real(n))


def ct2(n_l, n_l_1):
    n_l *= cos(theta(n_l))
    n_l_1 *= cos(theta(n_l_1))
    return 4 * n_l * n_l_1 / (n_l + n_l_1) ** 2


def cr_l_1_l(n_l, n_l_1):  # from n_l-1 to n_l
    n_l_1 *= cos(theta(n_l_1))
    n_l *= cos(theta(n_l))
    return (n_l_1 - n_l) / (n_l_1 + n_l)


def phase_factor(n, k, thick, freq):  # theta in radians
    omg = 2 * pi * freq
    thick *= cos(theta(n))
    phi = 2 * omg * thick / c_0
    return exp(- 1j * n * phi) * exp(- k * phi)


def recover_material(freqs, material):  # ['ABS', 'HDPE', 'PA6', 'PAMD', 'PC', 'PET', 'PFTE', 'PMMA', 'PP']
    fn_read, n_read = read_1file('./sim_resources/polymer_database/n_' + material + '.csv')
    fa_read, alpha_read = read_1file('./sim_resources/polymer_database/alpha_' + material + '.csv')
    n_intplt = intplt.interp1d(fn_read, n_read, bounds_error=False, fill_value=(n_read[0], n_read[-1]))
    n = n_intplt(freqs)
    alpha_intplt = intplt.interp1d(fa_read, alpha_read, bounds_error=False, fill_value=(alpha_read[0], alpha_read[-1]))
    k = 1e-10 * c_0 * alpha_intplt(freqs) / (4 * pi * freqs + eps)
    return n, k


def H_sim(freq, n_i, k_i, thick_i, n_m, k_m, thick_m, n_o, k_o, thick_o, d_air):
    H_i = cr_l_1_l(n_subs, n_i - 1j * k_i)

    rlm1l = cr_l_1_l(n_i - 1j * k_i, n_m - 1j * k_m)
    tt = ct2(n_i - 1j * k_i, n_m - 1j * k_m)
    exp_phi = phase_factor(n_i, k_i, thick_i, freq)

    H_i = rlm1l + (tt * H_i * exp_phi) / (1 + rlm1l * H_i * exp_phi)

    rlm1l = cr_l_1_l(n_m - 1j * k_m, n_o - 1j * k_o)
    tt = ct2(n_m - 1j * k_m, n_o - 1j * k_o)
    exp_phi = phase_factor(n_m, k_m, thick_m, freq)

    H_i = rlm1l + (tt * H_i * exp_phi) / (1 + rlm1l * H_i * exp_phi)

    rlm1l = cr_l_1_l(n_o - 1j * k_o, n_air_cplx)
    tt = ct2(n_o - 1j * k_o, n_air_cplx)
    exp_phi = phase_factor(n_o, k_o, thick_o, freq)

    H_i = rlm1l + (tt * H_i * exp_phi) / (1 + rlm1l * H_i * exp_phi)

    return exp(- 1j * 2 * 2 * pi * freq * d_air / c_0) * H_i


def cost_function(params, *args):
    d_air, thick_i, thick_m, thick_o = params
    E_sam, E_ref_w, freqs, n_i, k_i, n_m, k_m, n_o, k_o = args
    H_teo = H_sim(freqs, n_i, k_i, thick_i, n_m, k_m, thick_m, n_o, k_o, thick_o, d_air)
    E_teo = irfft(H_teo * E_ref_w)
    return sum((E_sam - E_teo) ** 2)


in_dir = './output/simulation_real_refs/3_layer/traces/'
out_dir = './output/simulation_real_refs/3_layer/'
ref_dir = './sim_resources/refs/'
t_ref, E_ref = read_1file(ref_dir + '100k_1.txt')
f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)
# plot(f_ref, toDb_0(E_ref_w))
# xlabel(r'$f\ (THz)$')
# xlim([0, 3])
# ylim([-80, 5])
# show()
# quit()
data_base_dir = './sim_resources/polymer_database/'
dir_list = os.listdir(in_dir)
rows = len(dir_list)
columns = 10
data = list()
fh = open('./output/simulation_real_refs/3_layer/results.txt', 'r')
i = 0
for line in fh:
    line = line.replace('\n', '')
    line = line.replace('k', '')
    line = line.split(' ')
    for j in range(columns):
        line[j] = float(line[j])
    data.append(tuple(line))

data = array(data)

sort_idxs = argsort(data, axis=0)[:, 0]
data = data[sort_idxs]
for i in range(4):
    # fig = figure(i+1)
    fig, axs = subplots(2, 1, sharex=True, gridspec_kw={
                           'height_ratios': [3, 1]})
    fig.subplots_adjust(hspace=0)
    # ax = fig.add_subplot(2, 1, 1)
    d_sims = data[16 * i:16 * i + 15, 1]
    d_fits_inner = data[16 * i:16 * i + 15, 2]
    d_fits_middle = data[16 * i:16 * i + 15, 4]
    d_fint_outer = data[16 * i:16 * i + 15, 6]
    axs[0].loglog(sort(d_sims), sort(d_sims), 'k-', lw=0.8, label='sims')
    axs[0].loglog(d_sims, d_fits_inner, '.', label='inner')
    axs[0].loglog(d_sims, d_fits_middle, '.', label='middle')
    axs[0].loglog(d_sims, d_fint_outer, '.', label='outer')
    # title(data[8 * i, 0])
    axs[0].legend()
    axs[0].set_ylabel(r'$d_{fit}\ (\mu m)$')
    # axs[0].set_box_aspect(1/0.2)
    # axs[1].loglog(sort(data[16 * i:16 * i + 15, 1]), sort(data[16 * i:16 * i + 15, 1]), 'k-', lw=0.8, label='sims')
    axs[1].loglog(d_sims, abs(d_sims - d_fits_inner) / d_sims, '.')
    axs[1].loglog(d_sims, abs(d_sims - d_fits_middle) / d_sims, '.')
    axs[1].loglog(d_sims, abs(d_sims - d_fint_outer) / d_sims, '.')
    axs[1].loglog(d_sims, ones(d_sims.size), 'k-', lw=0.5, label='0')
    # axs[1].loglog(d_sims, 10 * ones(d_sims.size), 'r-', lw=0.8, label=r'$10\ \mu m$')
    # axs[1].legend()
    axs[1].set_ylabel('desv')
    # axs[0].set_box_aspect(0.999)
    # axs[1].set_box_aspect(0.001)

    xlabel(r'$d_{sim}\ (\mu m)$')
    savefig(out_dir + str(data[16 * i, 0]).replace('.0','') + '.pdf')
show()
