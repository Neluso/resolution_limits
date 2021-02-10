from TDSA import *
import scipy.interpolate as intplt


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


def H_sim(freq, n_i, k_i, thick_i, n_o, k_o, thick_o, d_air):
    H_i = cr_l_1_l(n_subs, n_i - 1j * k_i)

    rlm1l = cr_l_1_l(n_i - 1j * k_i, n_o - 1j * k_o)
    tt = ct2(n_i - 1j * k_i, n_o - 1j * k_o)
    exp_phi = phase_factor(n_i, k_i, thick_i, freq)

    H_i = rlm1l + (tt * H_i * exp_phi) / (1 + rlm1l * H_i * exp_phi)

    rlm1l = cr_l_1_l(n_o - 1j * k_o, n_air_cplx)
    tt = ct2(n_o - 1j * k_o, n_air_cplx)
    exp_phi = phase_factor(n_o, k_o, thick_o, freq)

    H_i = rlm1l + (tt * H_i * exp_phi) / (1 + rlm1l * H_i * exp_phi)

    return exp(- 1j * 2 * 2 * pi * freq * d_air / c_0) * H_i


def sim_traces(mat_i, mat_o):
    t0 = time_ns()
    out_dir = './output/simulation_real_refs/2_layer/traces/'
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    in_dir = './sim_resources/refs/'
    ref_list = os.listdir(in_dir)
    for ref_file in ref_list:
        if ref_file.split('_')[1] == '1.txt':
            t_ref, E_ref = read_1file(in_dir + ref_file)  # t_ref in ps
            t_ref *= 1e-12
            f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)  # f_ref in THz
            n_i, k_i = recover_material(f_ref, mat_i)
            n_o, k_o = recover_material(f_ref, mat_o)
            for d_mat in pow(10, arange(-1, 3, 1 / 4)):
                name_trace = str(d_mat) + '_' + mat_i + '_' + mat_o + '_' + ref_file.replace('_1', '')
                d_mat *= 1e-6  # um
                H_sim_teo = H_sim(f_ref, n_i, k_i, d_mat, n_o, k_o, d_mat, 0.0)
                E_sim_w = H_sim_teo * E_ref_w
                E_sim = irfft(E_sim_w)
                tw = open(out_dir + name_trace, 'w')
                for i in range(E_sim.size):
                    tw.write(str(t_ref[i] * 1e12) + ',' + str(E_sim[i]) + '\n')
                tw.close()
        else:
            continue
    return 0


def sim_traces_v2(n_i_eff, n_o_eff, dispersion, k_eff_prima):
    out_dir = './output/simulation_real_refs/2_layer/traces/'
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    in_dir = './sim_resources/refs/'
    ref_file = '100k_1.txt'
    t_ref, E_ref = read_1file(in_dir + ref_file)  # t_ref in ps
    t_ref *= 1e-12
    f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)  # f_ref in THz
    f_ref *= 1e-12
    n_i = n_i_eff + dispersion * f_ref
    n_o = n_o_eff + dispersion * f_ref
    k_i = k_eff_prima * f_ref
    k_o = k_eff_prima * f_ref
    figure()
    plot(f_ref, n_i, f_ref, n_o)
    figure()
    plot(f_ref, k_i, f_ref, k_o)
    show()
    quit()
    f_ref *= 1e12
    step_d_sim = 1 / 6
    for d_mat in pow(10, arange(0, 2. + step_d_sim, step_d_sim)):
    # for d_mat in [1, pow(10, )]:
        name_trace = str(d_mat) + '_' + str(n_i_eff) + '_' + str(n_o_eff) + '_' + str(dispersion) + '_' +\
                     str(k_eff_prima) + '_' + ref_file.replace('_1', '')
        d_mat *= 1e-6  # um
        H_sim_teo = H_sim(f_ref, n_i, k_i, d_mat, n_o, k_o, d_mat, 0.0)
        E_sim_w = H_sim_teo * E_ref_w
        E_sim = irfft(E_sim_w)
        tw = open(out_dir + name_trace, 'w')
        for i in range(E_sim.size):
            tw.write(str(t_ref[i] * 1e12) + ',' + str(E_sim[i]) + '\n')
        tw.close()
    return 0


# sim_traces('PFTE', 'PA6')  # ['ABS', 'HDPE', 'PA6', 'PAMD', 'PC', 'PET', 'PFTE', 'PMMA', 'PP']
# sim_traces('PET', 'PC')  # ['ABS', 'HDPE', 'PA6', 'PAMD', 'PC', 'PET', 'PFTE', 'PMMA', 'PP']
# sim_traces('PFTE', 'PP')  # ['ABS', 'HDPE', 'PA6', 'PAMD', 'PC', 'PET', 'PFTE', 'PMMA', 'PP']
# sim_traces('PET', 'PP')  # ['ABS', 'HDPE', 'PA6', 'PAMD', 'PC', 'PET', 'PFTE', 'PMMA', 'PP']


for disps in range(3, 7):
    sim_traces_v2(1.4, 1.6, - 0.005 * disps, 0.025)
