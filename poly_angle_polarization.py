from TDSA import *
import scipy.interpolate as intplt
from scipy.optimize import differential_evolution

deg_in = 0  # incidence angle in degrees
snell_sin = n_air * sin(deg_in * pi / 180)
# n_subs = 1.17 - 0.0 * 1j  # substrate refractive index -- cork
n_subs = 1e20 - 0.0 * 1j  # substrate refractive index -- metal
eps = 1e-20


# function definitions
def theta(n, deg_in_air):
    snell_sin = n_air * sin(deg_in_air * pi / 180)
    return arcsin(snell_sin / real(n))


def ct2(n_l, n_l_1, deg_in_air, pol='s'):
    if pol == 's':
        n_l *= cos(theta(n_l, deg_in_air))
        n_l_1 *= cos(theta(n_l_1, deg_in_air))
    if pol == 'p':
        n_l *= cos(theta(n_l_1, deg_in_air))
        n_l_1 *= cos(theta(n_l, deg_in_air))
    return 4 * n_l * n_l_1 / (n_l + n_l_1) ** 2


def cr_l_1_l(n_l, n_l_1, deg_in_air, pol='s'):  # from n_l-1 to n_l
    if pol == 's':
        n_l *= cos(theta(n_l, deg_in_air))
        n_l_1 *= cos(theta(n_l_1, deg_in_air))
    if pol == 'p':
        n_l *= cos(theta(n_l_1, deg_in_air))
        n_l_1 *= cos(theta(n_l, deg_in_air))
    return (n_l_1 - n_l) / (n_l_1 + n_l)


def phase_factor(n, k, thick, deg_in_air, freq):  # theta in radians
    omg = 2 * pi * freq
    thick *= cos(theta(n, deg_in_air))
    phi = 2 * omg * thick / c_0
    return exp(- 1j * n * phi) * exp(- k * phi)


def H_sim(freq, n_i, k_i, thick_i, n_o, k_o, thick_o, d_air, deg_in_air, polaritz):
    H_i = cr_l_1_l(n_subs, n_i - 1j * k_i, deg_in_air, pol=polaritz)

    rlm1l = cr_l_1_l(n_i - 1j * k_i, n_o - 1j * k_o, deg_in_air, pol=polaritz)
    tt = ct2(n_i - 1j * k_i, n_o - 1j * k_o, deg_in_air, pol=polaritz)
    exp_phi = phase_factor(n_i, k_i, thick_i, deg_in_air, freq)

    H_i = rlm1l + (tt * H_i * exp_phi) / (1 + rlm1l * H_i * exp_phi)

    rlm1l = cr_l_1_l(n_o - 1j * k_o, n_air_cplx, deg_in_air, pol=polaritz)
    tt = ct2(n_o - 1j * k_o, n_air_cplx, deg_in_air, pol=polaritz)
    exp_phi = phase_factor(n_o, k_o, thick_o, deg_in_air, freq)

    H_i = rlm1l + (tt * H_i * exp_phi) / (1 + rlm1l * H_i * exp_phi)

    return exp(- 1j * 2 * 2 * pi * freq * d_air / c_0) * H_i


