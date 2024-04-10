import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os
from tqdm.auto import tqdm
from multiprocessing import Pool


grid_ionosphere_middle = 29
grid_middle_magnetosphere = 109
grid_fix = 175

BC_number = 21
min_number = 110

grid_number = 277

series_1 = 10   # Magnetospheric electrons
alpha_1 = r'5.E-1'
series_2 = 12   # Trapped electrons
alpha_2 = r'5.E-1'
series_3 = 4   # Ionospheric hot electrons
alpha_3 = r'1.E-1'

pdf_min = 1E-21

l_shell = 10E0
planet_radius = 6378.1E3
series_number = 12

speed_of_light = 299792458E0    #[m s-1]
elementary_charge = 1.6021766208E-19    #[A s]

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

file_velocity_distribution_function_name_1 = f'velocity_distribution_function/min_{str(min_number).zfill(3)}_grid_{str(grid_number).zfill(3)}_series_{str(series_1).zfill(2)}_alpha_{alpha_1}.csv'
path_velocity_distribution_function_name_1 = f'{dir_name}{dir_BC_name}{file_velocity_distribution_function_name_1}'

file_velocity_distribution_function_name_2 = f'velocity_distribution_function/min_{str(min_number).zfill(3)}_grid_{str(grid_number).zfill(3)}_series_{str(series_2).zfill(2)}_alpha_{alpha_2}.csv'
path_velocity_distribution_function_name_2 = f'{dir_name}{dir_BC_name}{file_velocity_distribution_function_name_2}'

file_velocity_distribution_function_name_3 = f'velocity_distribution_function/min_{str(min_number).zfill(3)}_grid_{str(grid_number).zfill(3)}_series_{str(series_3).zfill(2)}_alpha_{alpha_3}.csv'
path_velocity_distribution_function_name_3 = f'{dir_name}{dir_BC_name}{file_velocity_distribution_function_name_3}'

print(path_all_name)
print(path_velocity_distribution_function_name_1)
print(path_velocity_distribution_function_name_2)
print(path_velocity_distribution_function_name_3)

fig_name = f'plot/probably_density_function_trappedelectron_min_{str(min_number).zfill(3)}_grid_{str(grid_number).zfill(3)}_alpha_{alpha_1}_{alpha_2}_{alpha_3}'
path_fig_name = f'{dir_name}{dir_BC_name}{fig_name}.png'
os.makedirs(os.path.dirname(path_fig_name), exist_ok=True)

data_all = np.genfromtxt(path_all_name, delimiter=',', unpack=True)
number_density = np.array([data_all[6+series_1, grid_number-1], data_all[6+series_2, grid_number-1], data_all[6+series_3, grid_number-1]])
mlat_deg = data_all[3, grid_number]
length2planet_per_R = data_all[1, grid_number-1] / planet_radius

print(mlat_deg)
print(length2planet_per_R)

def input_velocity_distribution_function_data(args):
    path_name, path_number = args
    data = np.genfromtxt(path_name, delimiter=',', unpack=True)
    mlat_deg_array = data[0, :]
    v_perp_i = data[1, :] * 1E-3    # [km/s]
    v_para_i = data[2, :] * 1E-3    # [km/s]
    pdf = data[5, :] * number_density[path_number]  # velocity distribution function [m^-3 s^3]
    pdf_max = np.nanmax(pdf)
    return mlat_deg_array, v_perp_i, v_para_i, pdf, pdf_max

if __name__ == '__main__':
    args = [(path_velocity_distribution_function_name_1, 0), (path_velocity_distribution_function_name_2, 1), (path_velocity_distribution_function_name_3, 2)]
    with Pool(3) as p:
        results = list(tqdm(p.map(input_velocity_distribution_function_data, args), total=len(args)))
    
    mlat_deg_1_array, v_perp_i_1, v_para_i_1, pdf_1, pdf_max_1 = results[0]
    mlat_deg_2_array, v_perp_i_2, v_para_i_2, pdf_2, pdf_max_2 = results[1]
    mlat_deg_3_array, v_perp_i_3, v_para_i_3, pdf_3, pdf_max_3 = results[2]

pdf_max = max(pdf_max_1, pdf_max_2, pdf_max_3)

mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = r'\usepackage{bm}'
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = ['Computer Modern Roman']
mpl.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams["font.size"] = 30

fig = plt.figure(figsize=(15, 15), dpi=100)
gs = fig.add_gridspec(1, 2, width_ratios=[1, 0.05])
ax = fig.add_subplot(gs[0, 0], xlabel=r'Perpendicular velocity [$\mathrm{km/s}$]', ylabel=r'Field-aligned velocity [$\mathrm{km/s}$]', aspect='equal')
ax.set_title(r'$\mathrm{e}^{-}$ at Altitude $=$ ' + str(round(length2planet_per_R, 2)) + r' $\mathrm{R_E}$')

pdf_max_log10 = np.log10(pdf_max)
pdf_min_log10 = np.log10(pdf_min)
print(pdf_max_log10, pdf_min_log10)

cmap_color = cm.turbo
norm = mpl.colors.Normalize(vmin=pdf_min_log10, vmax=pdf_max_log10)
sm = cm.ScalarMappable(cmap=cmap_color, norm=norm)
sm.set_array([])
cbarax = fig.add_subplot(gs[0, 1])
cbar = plt.colorbar(sm, cax=cbarax)
cbar.set_label(r'Phase space density $\mathrm{log}_{10} \left( [\mathrm{s^3/m^6}] \right)$')


def plot_distribution_function(args):
    ax, v_perp_i, v_para_i, pdf, pdf_min = args
    pdf[pdf < pdf_min] = np.nan
    v_perp_i = np.where(np.isnan(pdf), np.nan, v_perp_i)
    v_para_i = np.where(np.isnan(pdf), np.nan, v_para_i)
    ax.scatter(v_perp_i, v_para_i, c=np.log10(pdf), cmap=cmap_color, norm=norm, s=50, alpha=0.7)
    v_perp_max = np.nanmax(v_perp_i)
    v_perp_min = np.nanmin(v_perp_i)
    v_para_max = np.nanmax(v_para_i)
    v_para_min = np.nanmin(v_para_i)
    return v_perp_max, v_perp_min, v_para_max, v_para_min

args = [(ax, v_perp_i_1, v_para_i_1, pdf_1, pdf_min), (ax, v_perp_i_2, v_para_i_2, pdf_2, pdf_min), (ax, v_perp_i_3, v_para_i_3, pdf_3, pdf_min)]
v_perp_max_array = []
v_perp_min_array = []
v_para_max_array = []
v_para_min_array = []
for arg in args:
    v_perp_max, v_perp_min, v_para_max, v_para_min = plot_distribution_function(arg)
    v_perp_max_array.append(v_perp_max)
    v_perp_min_array.append(v_perp_min)
    v_para_max_array.append(v_para_max)
    v_para_min_array.append(v_para_min)
v_perp_max = max(v_perp_max_array)
v_perp_min = min(v_perp_min_array)
v_para_max = max(v_para_max_array)
v_para_min = min(v_para_min_array)

ax.minorticks_on()
ax.grid(which='both', alpha=0.3)

# 400 eVごとに等高線を引く
v_perp_contour = np.linspace(v_perp_min, v_perp_max, 10000) * 1E3   # [m/s]
v_para_contour = np.linspace(v_para_min, v_para_max, 10000) * 1E3   # [m/s]
v_perp_mesh, v_para_mesh = np.meshgrid(v_perp_contour, v_para_contour)
v_para_max, v_perp_max = v_para_max * 1E3, v_perp_max * 1E3

Lorentz_factor_mesh = 1E0 / np.sqrt(1E0 - (v_perp_mesh**2 + v_para_mesh**2) / speed_of_light**2)
kinetic_energy_mesh_keV = mass_series[2] * speed_of_light**2 * (Lorentz_factor_mesh - 1E0) / elementary_charge * 1E-3

Lorentz_factor_max = 1E0 / np.sqrt(1E0 - (v_perp_max**2) / speed_of_light**2)
kinetic_energy_max_keV = mass_series[2] * speed_of_light**2 * (Lorentz_factor_max - 1E0) / elementary_charge * 1E-3

print(v_perp_max, Lorentz_factor_max, kinetic_energy_max_keV)

label_cont = np.arange(0, kinetic_energy_max_keV, 0.4)
cont = ax.contour(v_perp_mesh*1E-3, v_para_mesh*1E-3, kinetic_energy_mesh_keV, label_cont, colors='black', alpha=0.5)
cont.clabel(fmt='%1.1f keV', inline=True)
plt.tight_layout()
plt.savefig(path_fig_name)
plt.close()