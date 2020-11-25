from TDSA import *
import os


in_dir = './output/simulation_results/test_1/traces/'
dir_list = os.listdir(in_dir)


t_ref, E_ref = read_1file('./output/simulation_results/test_1/refs/-60_ref.txt')
f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)
# plot(f_ref, unwrap(angle(E_ref_w)), lw=0.5)
# plot(f_ref, toDb_0(E_ref_w), lw=0.5)
plot(t_ref, E_ref, lw=1, label='ref')


for file_i in dir_list:
    t_sim, E_sim = read_1file(in_dir + file_i)  # t_ref in s
    plot(t_sim, E_sim, lw=1, label=file_i)
    f_sim, E_sim_w = fourier_analysis(t_sim, E_sim)  # f_ref in Hz
legend()
show()
