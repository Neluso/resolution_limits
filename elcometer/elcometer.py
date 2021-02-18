from TDSA import *
import scipy.interpolate as intplt
from scipy.optimize import differential_evolution


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


def H_sim(freq, n, k, thick, d_air):  # n = 1.8, alpha
    H_i = cr_l_1_l(n_subs, n - 1j * k)

    rlm1l = cr_l_1_l(n - 1j * k, n_air_cplx)
    tt = ct2(n - 1j * k, n_air_cplx)
    exp_phi = phase_factor(n, k, thick, freq)

    H_i = rlm1l + (tt * H_i * exp_phi) / (1 + rlm1l * H_i * exp_phi)

    return exp(- 1j * 2 * 2 * pi * freq * d_air / c_0) * H_i


def cost_function(params, *args):
    d_air, thick = params
    E_sam, E_ref_w, freqs, n, k = args
    H_teo = H_sim(freqs, n, k, thick, d_air)
    E_teo = - irfft(H_teo * E_ref_w)
    return sum((E_sam - E_teo) ** 2)


t0 = time_ns()
name_files = ['morat', 'cian', 'marro', 'blanc']  # , 'negre', 'gris']
d_film = {'morat': 23.6e-6, 'cian': 50.0e-6, 'marro': 123.2e-6, 'blanc': 244.8e-6, 'negre': 463e-6, 'gris': 978e-6}
mat_film = {'morat': 'polyester', 'cian': 'polyester', 'marro': 'polyester', 'blanc': 'polyester', 'negre': 463e-6, 'gris': 978e-6}

if __name__ == '__main__':
    for name_file in name_files:
        cfh = open('./calibration_' + mat_film[name_file] + '.txt')
        freq_cal = list()
        n_cal = list()
        alpha_cal = list()
        for line in cfh.read().split('\n'):
            line = line.split(',')
            try:
                freq_cal.append(float(line[0]))
                n_cal.append(float(line[1]))
                alpha_cal.append(float(line[3]))
            except:
                continue

        freq_cal = array(freq_cal)
        n_cal = array(n_cal)
        alpha_cal = array(alpha_cal)
        n_intplt = intplt.interp1d(freq_cal, n_cal, bounds_error=False, fill_value=(n_cal[0], n_cal[-1]))
        alpha_intplt = intplt.interp1d(freq_cal, alpha_cal, bounds_error=False, fill_value=(alpha_cal[0], alpha_cal[-1]))

        for j in range(2):
            t_ref, E_ref = read_1file('./' + name_file + '_ref_' + str(j + 1) + '.txt')
            t_sam, E_sam = read_1file('./' + name_file + '_sam_' + str(j + 1) + '.txt')
            enlargement = 0 * t_ref.size
            delta_t_ref = mean(diff(t_ref))
            E_ref = zero_padding(E_ref, 0, enlargement)
            t_ref = concatenate((t_ref, t_ref[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
            E_sam = zero_padding(E_sam, 0, enlargement)
            t_sam = concatenate((t_sam, t_sam[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
            t_ref *= 1e-12
            t_sam *= 1e-12
            f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)
            f_sam, E_sam_w = fourier_analysis(t_ref, E_ref)

            k_bounds = [
                (-1e-3, 1e-3),  # d_air
                # (1.79, 1.81),  # n
                # (0, 1),  # k
                (1e-6, 1e-3)  # thickness
            ]

            n = n_intplt(f_ref)
            k = 1e-10 * c_0 * alpha_intplt(f_ref) / (4 * pi * f_ref + eps)

            num_statistics = 50
            delta_error = 0.01
            error_mod = 1 + delta_error * (2 * random.rand(num_statistics) - 1)
            resx0 = list()
            # resx1 = list()
            # resx2 = list()
            resx3 = list()
            for i in range(num_statistics):
                res = differential_evolution(cost_function,
                                             k_bounds,
                                             args=(E_sam, E_ref_w, f_ref, n * error_mod[i], k * error_mod[i]),
                                             # popsize=30,
                                             # maxiter=3000,
                                             updating='deferred',
                                             workers=-1,
                                             disp=False,  # step cost_function value
                                             polish=True
                                             )
                resx0.append(res.x[0])
                # resx1.append(res.x[1])
                # resx2.append(res.x[2])
                resx3.append(res.x[1])
            resx0 = array(resx0)
            # resx1 = array(resx1)
            # resx2 = array(resx2)
            resx3 = array(resx3) * 1e6
            print('Sample', name_file, j+1)
            # print('n = ', mean(resx1), '+-', std(resx1))
            # print('k = ', mean(resx2), '+-', std(resx2))
            print('d = ', mean(resx3), '+-', std(resx3))
            H_fit = H_sim(f_ref, n, k, mean(resx3) * 1e-6, - mean(resx3) * 1e-6)
            H_teo = H_sim(f_ref, n, k, d_film[name_file], - d_film[name_file])
            figure(name_file + '_' + str(j+1))
            plot(t_sam, E_sam, label='medida')
            plot(t_ref, - irfft(H_fit * E_ref_w), label='ajuste')
            plot(t_ref, - irfft(H_teo * E_ref_w), label='calculada')
            legend()
show()
