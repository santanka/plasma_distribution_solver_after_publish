import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import stats
import os

grid_ionosphere_middle = 14
grid_middle_magnetosphere = 109
grid_fix = 175

BC_number = 12
min_number = 149

grid_number = 242

series_1 = 10   # Magnetospheric electrons
series_2 = 12   # Trapped electrons
series_3 = 4   # Ionospheric hot electrons

l_shell = 10E0
planet_radius = 6378.1E3
series_number = 12

if (series_number == 12):
    mass_series = np.array([2.677950266103E-26, 1.67262192369E-27, 9.1093837015E-31, 9.1093837015E-31, 2.677950266103E-26, 1.67262192369E-27, 9.1093837015E-31, 9.1093837015E-31, 1.67262192369E-27, 9.1093837015E-31, 9.1093837015E-31, 9.1093837015E-31])
    charge_series = np.array([1E0, 1E0, -1E0, -1E0, 1E0, 1E0, -1E0, -1E0, 1E0, -1E0, -1E0, -1E0])
elif (series_number == 10):
    mass_series = np.array([2.677950266103E-26, 1.67262192369E-27, 9.1093837015E-31, 2.677950266103E-26, 1.67262192369E-27, 9.1093837015E-31, 1.67262192369E-27, 9.1093837015E-31, 9.1093837015E-31, 9.1093837015E-31])
    charge_series = np.array([1E0, 1E0, -1E0, 1E0, 1E0, -1E0, 1E0, -1E0, -1E0, -1E0])

dir_name = f'/mnt/j/plasma_distribution_solver_after_publish/Earth_L_10_Imajo/alpha_perp_12_parallel_12/grid_{str(grid_ionosphere_middle).zfill(3)}_{str(grid_middle_magnetosphere).zfill(3)}_{str(grid_fix).zfill(3)}/'
dir_BC_name = f'boundary_condition_{str(BC_number)}/'

file_all_name = f'all/all_{str(min_number).zfill(3)}.csv'
path_all_name = f'{dir_name}{dir_BC_name}{file_all_name}'

file_velocity_distribution_function_name_1 = f'velocity_distribution_function/min_{str(min_number).zfill(3)}_grid_{str(grid_number).zfill(3)}_series_{str(series_1).zfill(2)}.csv'
path_velocity_distribution_function_name_1 = f'{dir_name}{dir_BC_name}{file_velocity_distribution_function_name_1}'

file_velocity_distribution_function_name_2 = f'velocity_distribution_function/min_{str(min_number).zfill(3)}_grid_{str(grid_number).zfill(3)}_series_{str(series_2).zfill(2)}.csv'
path_velocity_distribution_function_name_2 = f'{dir_name}{dir_BC_name}{file_velocity_distribution_function_name_2}'

file_velocity_distribution_function_name_3 = f'velocity_distribution_function/min_{str(min_number).zfill(3)}_grid_{str(grid_number).zfill(3)}_series_{str(series_3).zfill(2)}.csv'
path_velocity_distribution_function_name_3 = f'{dir_name}{dir_BC_name}{file_velocity_distribution_function_name_3}'

print(path_all_name)
print(path_velocity_distribution_function_name_1)
print(path_velocity_distribution_function_name_2)
print(path_velocity_distribution_function_name_3)

data_all = np.genfromtxt(path_all_name, delimiter=',', unpack=True)
number_density_1 = data_all[6+series_1, grid_number-1]
number_density_2 = data_all[6+series_2, grid_number-1]
number_density_3 = data_all[6+series_3, grid_number-1]
mlat_deg = data_all[3, grid_number]
length2planet_per_R = data_all[1, grid_number-1] / planet_radius

print(mlat_deg)
print(length2planet_per_R)
print(number_density_1*1E-6)
print(number_density_2*1E-6)
print(number_density_3*1E-6)

data_velocity_distribution_function_1 = np.genfromtxt(path_velocity_distribution_function_name_1, delimiter=',', unpack=True)
data_velocity_distribution_function_2 = np.genfromtxt(path_velocity_distribution_function_name_2, delimiter=',', unpack=True)
data_velocity_distribution_function_3 = np.genfromtxt(path_velocity_distribution_function_name_3, delimiter=',', unpack=True)

mlat_deg_1_array            = data_velocity_distribution_function_1[0, :]
v_perp_i_1                  = data_velocity_distribution_function_1[1, :] * 1E-3 # [km/s]
v_para_i_1                  = data_velocity_distribution_function_1[2, :] * 1E-3 # [km/s]
pdf_1                       = data_velocity_distribution_function_1[5, :]   # Probability distribution function [m^-3 s^3]
vdf_1                       = pdf_1 * number_density_1 * 1E3                # Velocity distribution function [cm^-3 (km^-1 s)^3]

mlat_deg_2_array            = data_velocity_distribution_function_2[0, :]
v_perp_i_2                  = data_velocity_distribution_function_2[1, :] * 1E-3 # [km/s]
v_para_i_2                  = data_velocity_distribution_function_2[2, :] * 1E-3 # [km/s]
pdf_2                       = data_velocity_distribution_function_2[5, :]   # Probability distribution function [m^-3 s^3]
vdf_2                       = pdf_2 * number_density_2 * 1E3                # Velocity distribution function [cm^-3 (km^-1 s)^3]

mlat_deg_3_array            = data_velocity_distribution_function_3[0, :]
v_perp_i_3                  = data_velocity_distribution_function_3[1, :] * 1E-3 # [km/s]
v_para_i_3                  = data_velocity_distribution_function_3[2, :] * 1E-3 # [km/s]
pdf_3                       = data_velocity_distribution_function_3[5, :]   # Probability distribution function [m^-3 s^3]
vdf_3                       = pdf_3 * number_density_3 * 1E3                # Velocity distribution function [cm^-3 (km^-1 s)^3]

vdf_array = np.append(vdf_1, vdf_2)
vdf_array = np.append(vdf_array, vdf_3)
v_perp_array = np.append(v_perp_i_1, v_perp_i_2)
v_perp_array = np.append(v_perp_array, v_perp_i_3)
v_para_array = np.append(v_para_i_1, v_para_i_2)
v_para_array = np.append(v_para_array, v_para_i_3)

vdf_max = np.nanmax(vdf_array)
vdf_min = vdf_max * 1E-20
vdf_array[vdf_array < vdf_min] = np.nan
v_perp_array = np.where(np.isnan(vdf_array), np.nan, v_perp_array)
v_para_array = np.where(np.isnan(vdf_array), np.nan, v_para_array)
vdf_min = np.nanmin(vdf_array)

#vperp minus
vdf_array = np.append(vdf_array, vdf_array)
v_perp_array = np.append(v_perp_array, -v_perp_array)
v_para_array = np.append(v_para_array, v_para_array)

mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = r'\usepackage{bm}'
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = ['Computer Modern Roman']
mpl.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams["font.size"] = 40

fig = plt.figure(figsize=(25, 15), dpi=300)
ax = fig.add_subplot(111, xlabel=r'$v_{\parallel i}$ [km/s] ($+$ : South $\rightarrow$ North, $-$ : North $\rightarrow$ South)', ylabel=r'$v_{\perp i}$ [km/s]', aspect='equal')
ax.set_title(r'Electron distribution function at Altitude $=$ ' + str(round(length2planet_per_R, 2)) + r' $\mathrm{R_E}$')
mappable = ax.scatter(v_para_array, v_perp_array, c=np.log10(vdf_array), cmap='turbo', vmin=np.log10(vdf_min), vmax=np.log10(vdf_max), s=50, alpha=0.7)
fig.colorbar(mappable=mappable, ax=ax, label=r'$\log_{10} (f(r_{\parallel i}, \bm{v}_{i}) [\mathrm{cm}^{-3} (\mathrm{km}^{-1} \, \mathrm{s})^{3}])$')
ax.minorticks_on()
ax.grid(which='both', alpha=0.3)
plt.tight_layout()

fig_name = f'plot/probably_density_function_trappedelectron_min_{str(min_number).zfill(3)}_grid_{str(grid_number).zfill(3)}'
path_fig_name = f'{dir_name}{dir_BC_name}{fig_name}.png'
os.makedirs(os.path.dirname(path_fig_name), exist_ok=True)
plt.savefig(path_fig_name)
plt.close()