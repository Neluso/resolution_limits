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


def cost_function(params, *args):
    d_air, thick_i, thick_o = params
    E_sam, E_ref_w, freqs, n_i, k_i, n_o, k_o = args
    H_teo = H_sim(freqs, n_i, k_i, thick_i, n_o, k_o, thick_o, d_air)
    E_teo = irfft(H_teo * E_ref_w)
    return sum((E_sam - E_teo) ** 2)


