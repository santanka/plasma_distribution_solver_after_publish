import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import stats
import os

grid_ionosphere_middle = 14
grid_middle_magnetosphere = 109
grid_fix = 175

BC_number = 11
min_number = 155

grid_number = 277

series_1 = 10
series_1_name = r'$\mathrm{e^{-}}$ (Magnetosphere) ($v_{\parallel i} > 0$)'
series_2 = 12
series_2_name = r'$\mathrm{e^{-}}$ (Trapped)'

channel = 1

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

file_name_1 = f'velocity_distribution_function/min_{str(min_number).zfill(3)}_grid_{str(grid_number).zfill(3)}_series_{str(series_1).zfill(2)}.csv'
path_name_1 = f'{dir_name}{dir_BC_name}{file_name_1}'

file_name_2 = f'velocity_distribution_function/min_{str(min_number).zfill(3)}_grid_{str(grid_number).zfill(3)}_series_{str(series_2).zfill(2)}.csv'
path_name_2 = f'{dir_name}{dir_BC_name}{file_name_2}'

print(path_name_1)
print(path_name_2)

data_1 = np.genfromtxt(path_name_1, delimiter=',', unpack=True)
data_2 = np.genfromtxt(path_name_2, delimiter=',', unpack=True)

mlat_deg_1_array            = data_1[0, :]
v_perp_i_1                  = data_1[1, :]
v_para_i_1                  = data_1[2, :]
v_perp_b_1                  = data_1[3, :]
v_para_b_1                  = data_1[4, :]
distribution_function_1     = data_1[5, :]
differential_flux_1         = data_1[6, :]
mo, _ = stats.mode(mlat_deg_1_array)
mlat_deg_1 = mo[0]
length2planet_1 = planet_radius * l_shell * np.cos(np.deg2rad(mlat_deg_1))**2E0

mlat_deg_2_array            = data_2[0, :]
v_perp_i_2                  = data_2[1, :]
v_para_i_2                  = data_2[2, :]
v_perp_b_2                  = data_2[3, :]
v_para_b_2                  = data_2[4, :]
distribution_function_2     = data_2[5, :]
differential_flux_2         = data_2[6, :]
mo, _ = stats.mode(mlat_deg_2_array)
mlat_deg_2 = mo[0]
length2planet_2 = planet_radius * l_shell * np.cos(np.deg2rad(mlat_deg_2))**2E0

def make_nan(v_perp, v_para, distribution_function):
    length = len(distribution_function)
    max_distribution_function = np.nanmax(distribution_function)
    min_distribution_function = np.nanmin(distribution_function)
    if (np.log10(min_distribution_function) < np.log10(max_distribution_function) - 20.):
        for count_i in range(length):
            if(np.log10(distribution_function[count_i]) < np.log10(max_distribution_function)-20.):
                v_perp[count_i] = np.nan
                v_para[count_i] = np.nan
                distribution_function[count_i] = np.nan
    return v_perp, v_para, distribution_function

def make_nan_vpara(v_perp, v_para, distribution_function):
    length = len(distribution_function)
    max_distribution_function = np.nanmax(distribution_function)
    min_distribution_function = np.nanmin(distribution_function)
    if (np.log10(min_distribution_function) < np.log10(max_distribution_function) - 20.):
        for count_i in range(length):
            if(v_para[count_i] < 0.):
                v_perp[count_i] = np.nan
                v_para[count_i] = np.nan
                distribution_function[count_i] = np.nan
    return v_perp, v_para, distribution_function

mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = r'\usepackage{bm}'
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = ['Computer Modern Roman']
mpl.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams["font.size"] = 20

fig = plt.figure(figsize=(14, 18), dpi=100, tight_layout=True)
fig.suptitle(r"Probability Density Function (scale=log10)")

def plot_distribution_function_nonvpara(place, label_x, label_y, axtitle, cbar_title, v_perp, v_para, distribution_function):
    v_perp, v_para, distribution_function = make_nan(v_perp, v_para, distribution_function)
    ax = fig.add_subplot(place, xlabel=label_x, ylabel=label_y, title=axtitle)
    mappable = ax.scatter(v_para, v_perp, c=np.log10(distribution_function), vmin=np.floor(np.nanmin(np.log10(distribution_function))), vmax=np.trunc(np.nanmax(np.log10(distribution_function))), cmap='turbo', s=700, alpha=0.7)
    fig.colorbar(mappable=mappable, ax=ax, label=cbar_title)
    ax.minorticks_on()
    ax.grid(which='both', alpha=0.3)
    ax.set_axisbelow(True)

def plot_distribution_function_vpara(place, label_x, label_y, axtitle, cbar_title, v_perp, v_para, distribution_function):
    v_perp, v_para, distribution_function = make_nan_vpara(v_perp, v_para, distribution_function)
    v_perp, v_para, distribution_function = make_nan(v_perp, v_para, distribution_function)
    ax = fig.add_subplot(place, xlabel=label_x, ylabel=label_y, title=axtitle)
    mappable = ax.scatter(v_para, v_perp, c=np.log10(distribution_function), vmin=np.floor(np.nanmin(np.log10(distribution_function))), vmax=np.trunc(np.nanmax(np.log10(distribution_function))), cmap='turbo', s=700, alpha=0.7)
    fig.colorbar(mappable=mappable, ax=ax, label=cbar_title)
    ax.minorticks_on()
    ax.grid(which='both', alpha=0.3)
    ax.set_axisbelow(True)


if (channel == 1):
    axlabel_x = r'$v_{\parallel i}$ [$\mathrm{m/s}$] ($+$ : South $\rightarrow$ North, $-$ : North $\rightarrow$ South)'
    axlabel_y = r'$v_{\perp i}$ [$\mathrm{m/s}$]'

    axtitle_1 = r'(a) ' + series_1_name + r' at Altitude $=$ ' + str(round(length2planet_1 / planet_radius, 2)) + r'$\mathrm{R_E}$'
    axtitle_2 = r'(b) ' + series_2_name + r' at Altitude $=$ ' + str(round(length2planet_2 / planet_radius, 2)) + r'$\mathrm{R_E}$'

    cbar_title = r'$\log_{10} (f(r_{\parallel i}, \bm{v}_{i}) / n(r_{\parallel i}))$'

    plot_distribution_function_vpara(211, axlabel_x, axlabel_y, axtitle_1, cbar_title, v_perp_i_1, v_para_i_1, distribution_function_1)
    plot_distribution_function_nonvpara(212, axlabel_x, axlabel_y, axtitle_2, cbar_title, v_perp_i_2, v_para_i_2, distribution_function_2)

    file_name_3 = f'plot/probably_density_function_comparison_min_{str(min_number).zfill(3)}_grid_{str(grid_number).zfill(3)}_series_{str(series_1).zfill(2)}_and_{str(series_2).zfill(2)}'
    plt.savefig(f'{dir_name}{dir_BC_name}{file_name_3}.png')
    #plt.show()