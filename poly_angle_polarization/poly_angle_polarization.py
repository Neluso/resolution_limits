from TDSA import *
import scipy.interpolate as intplt
from scipy.optimize import differential_evolution

# deg_in = 0  # incidence angle in degrees
# snell_sin = n_air * sin(deg_in * pi / 180)
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


def cost_function(params, *args):
    d_air, thick_i, thick_o = params
    E_sam1, E_sam2, E_ref_w, freqs, n_i, k_i, n_o, k_o, inc_angle1, inc_angle2, polaritz = args
    H_teo1 = H_sim(freqs, n_i, k_i, thick_i, n_o, k_o, thick_o, d_air, inc_angle1, polaritz)
    H_teo2 = H_sim(freqs, n_i, k_i, thick_i, n_o, k_o, thick_o, d_air, inc_angle2, polaritz)
    E_teo1 = irfft(H_teo1 * E_ref_w)
    E_teo2 = irfft(H_teo2 * E_ref_w)
    return sum((E_sam1 - E_teo1) ** 2 + (E_sam2 - E_teo2) ** 2)


# t_ref, E_ref = read_1file('./ref.txt')
t_ref, E_ref = read_1file('./100k_1.txt')
t_ref *= 1e-12
f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)

n_i = 1.6 * ones(f_ref.size)
n_o = 1.5 * ones(f_ref.size)
k_i = 0.025 * f_ref * 1e-12
k_o = 0.025 * f_ref * 1e-12


d_mat = 20 * 1e-6  # 35 um


H_teo_15s = H_sim(f_ref, n_i, k_i, d_mat, n_o, k_o, d_mat, 0, 15, 's')
H_teo_15p = H_sim(f_ref, n_i, k_i, d_mat, n_o, k_o, d_mat, 0, 15, 'p')
H_teo_20s = H_sim(f_ref, n_i, k_i, d_mat, n_o, k_o, d_mat, 0, 20, 's')
H_teo_20p = H_sim(f_ref, n_i, k_i, d_mat, n_o, k_o, d_mat, 0, 20, 'p')
H_teo_25s = H_sim(f_ref, n_i, k_i, d_mat, n_o, k_o, d_mat, 0, 25, 's')
H_teo_25p = H_sim(f_ref, n_i, k_i, d_mat, n_o, k_o, d_mat, 0, 25, 'p')
H_teo_30s = H_sim(f_ref, n_i, k_i, d_mat, n_o, k_o, d_mat, 0, 30, 's')
H_teo_30p = H_sim(f_ref, n_i, k_i, d_mat, n_o, k_o, d_mat, 0, 30, 'p')
H_teo_35s = H_sim(f_ref, n_i, k_i, d_mat, n_o, k_o, d_mat, 0, 35, 's')
H_teo_35p = H_sim(f_ref, n_i, k_i, d_mat, n_o, k_o, d_mat, 0, 35, 'p')
H_teo_40s = H_sim(f_ref, n_i, k_i, d_mat, n_o, k_o, d_mat, 0, 40, 's')
H_teo_40p = H_sim(f_ref, n_i, k_i, d_mat, n_o, k_o, d_mat, 0, 40, 'p')
H_teo_45s = H_sim(f_ref, n_i, k_i, d_mat, n_o, k_o, d_mat, 0, 45, 's')
H_teo_45p = H_sim(f_ref, n_i, k_i, d_mat, n_o, k_o, d_mat, 0, 45, 'p')


E_sim_15s = irfft(H_teo_15s * E_ref_w)
E_sim_15p = irfft(H_teo_15p * E_ref_w)
E_sim_20s = irfft(H_teo_20s * E_ref_w)
E_sim_20p = irfft(H_teo_20p * E_ref_w)
E_sim_25s = irfft(H_teo_25s * E_ref_w)
E_sim_25p = irfft(H_teo_25p * E_ref_w)
E_sim_30s = irfft(H_teo_30s * E_ref_w)
E_sim_30p = irfft(H_teo_30p * E_ref_w)
E_sim_35s = irfft(H_teo_35s * E_ref_w)
E_sim_35p = irfft(H_teo_35p * E_ref_w)
E_sim_40s = irfft(H_teo_40s * E_ref_w)
E_sim_40p = irfft(H_teo_40p * E_ref_w)
E_sim_45s = irfft(H_teo_45s * E_ref_w)
E_sim_45p = irfft(H_teo_45p * E_ref_w)

# plot(t_ref, E_sim_15s, label='15s')
# plot(t_ref, E_sim_15p, label='15p')
# plot(t_ref, E_sim_45s, label='45s')
# plot(t_ref, E_sim_45p, label='45p')
# legend()
# show()

# trace_list = [E_sim_15s, E_sim_15p, E_sim_30s, E_sim_30p, E_sim_45s, E_sim_45p]
# trace_ID = ['15s', '15p', '30s', '30p', '45s', '45p']
# trace_list = [E_sim_20s, E_sim_20p, E_sim_25s, E_sim_25p, E_sim_35s, E_sim_35p, E_sim_40s, E_sim_40p]
# trace_ID = ['20s', '20p', '25s', '25p', '35s', '35p', '40s', '40p']

E_sam1 = E_sim_25s
E_sam2 = E_sim_30s
inc_angle1 = 25
inc_angle2 = 30
polaritz = 's'

wh = open('./output/' + str(d_mat*1e6) + '_results.txt', 'a')

if __name__ == '__main__':
    k_bounds = [
        (-1e-9, 1e-9),  # d_air
        (0, 1e-3),  # d_mat
        (0, 1e-3)  # d_mat
    ]
    num_statistics = 10
    delta_error = 0.01
    error_mod = 1 + delta_error * (2 * random.rand(num_statistics) - 1)
    resx0 = list()
    resx1 = list()
    resx2 = list()
    for i in range(num_statistics):
        res = differential_evolution(cost_function,
                                     k_bounds,
                                     args=(E_sam1, E_sam2, E_ref_w, f_ref,
                                           n_i * error_mod[i], k_i * error_mod[i],
                                           n_o * error_mod[i], k_o * error_mod[i],
                                           inc_angle1, inc_angle2, polaritz),
                                     popsize=30,
                                     maxiter=3000,
                                     updating='deferred',
                                     workers=-1,
                                     disp=True,  # step cost_function value
                                     polish=True
                                     )
        resx0.append(res.x[0])
        resx1.append(res.x[1])
        resx2.append(res.x[2])
    resx0 = array(resx0)
    resx1 = array(resx1)
    resx2 = array(resx2)

    wh.write(str(inc_angle1) + ' ' + str(inc_angle2) + ' ' + polaritz + ' '
             + str(mean(resx1) * 1e6) + ' ' + str(std(resx1) * 1e6) + ' '
             + str(mean(resx2) * 1e6) + ' ' + str(std(resx2) * 1e6) + '\n'
             # + str(mean(resx0) * 1e6) + ' ' + str(std(resx0) * 1e6) + '\n'
             )
