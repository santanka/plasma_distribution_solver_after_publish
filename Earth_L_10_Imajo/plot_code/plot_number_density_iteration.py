import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os

grid_ionosphere_middle = 14
grid_middle_magnetosphere = 109
grid_fix = 175

BC_number = 18
min_number = 174

l_shell = 10E0
series_number = 12

if (series_number == 12):
    mass_series = np.array([2.677950266103E-26, 1.67262192369E-27, 9.1093837015E-31, 9.1093837015E-31, 2.677950266103E-26, 1.67262192369E-27, 9.1093837015E-31, 9.1093837015E-31, 1.67262192369E-27, 9.1093837015E-31, 9.1093837015E-31, 9.1093837015E-31])
    charge_series = np.array([1E0, 1E0, -1E0, -1E0, 1E0, 1E0, -1E0, -1E0, 1E0, -1E0, -1E0, -1E0])
elif (series_number == 10):
    mass_series = np.array([2.677950266103E-26, 1.67262192369E-27, 9.1093837015E-31, 2.677950266103E-26, 1.67262192369E-27, 9.1093837015E-31, 1.67262192369E-27, 9.1093837015E-31, 9.1093837015E-31, 9.1093837015E-31])
    charge_series = np.array([1E0, 1E0, -1E0, 1E0, 1E0, -1E0, 1E0, -1E0, -1E0, -1E0])

dir_name = f'/mnt/j/plasma_distribution_solver_after_publish/Earth_L_10_Imajo/alpha_perp_12_parallel_12/grid_{str(grid_ionosphere_middle).zfill(3)}_{str(grid_middle_magnetosphere).zfill(3)}_{str(grid_fix).zfill(3)}/'
dir_BC_name = f'boundary_condition_{str(BC_number)}/'
file_name = f'number_density_iteration/min_{str(min_number).zfill(3)}.csv'

path_file = f'{dir_name}{dir_BC_name}{file_name}'
path_figure = f'{path_file[:-4]}.png'

data = np.loadtxt(path_file, delimiter=',', skiprows=1, max_rows=2*grid_fix-1)

data_length = data.shape[0]
data_length_half = int(data_length / 2)

coordinate_FA_half                      = data[data_length_half:, 0]
length2planet_half                      = data[data_length_half:, 1]
mlat_rad_half                           = data[data_length_half:, 2]
mlat_degree_half                        = data[data_length_half:, 3]
magnetic_flux_density_half              = data[data_length_half:, 4]
initial_electrostatic_potential_half    = data[data_length_half:, 5]
electrostatic_potential_half            = data[data_length_half:, 6]
number_density_half                     = data[data_length_half:, 7:series_number + 7]
charge_density_half                     = data[data_length_half:, series_number + 7]
charge_density_Poisson_half             = data[data_length_half:, series_number + 8]
convergence_number_half                 = data[data_length_half:, series_number + 9]

#constant parameter
elementary_charge = 1.602176634E-19 #[s A]
speed_of_light = 299792458E0    #[m s-1]
magnetic_constant = 1.25663706212E-6    #[kg m s-2 A-2]
electric_constant = 1E0 / magnetic_constant / speed_of_light**2E0   #[kg-1 m-3 s4 A2]
constant_of_gravitation = 6.67430E-11   #[kg-1 m3 s-2]
electron_mass = 9.1093837015E-31    #[kg]

planet_mass = 5.97243E24     #[kg]
planet_radius = 6.3781E6    #[m]

length2planet_half_planet_radius = length2planet_half / planet_radius

mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = ['Computer Modern Roman']
mpl.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams["font.size"] = 20

fig = plt.figure(figsize=(14, 28), dpi=100, tight_layout=True)
gs = fig.add_gridspec(16, 7)

ax1 = fig.add_subplot(gs[0:5, 0:6], xlabel=r'Geocentric distance [$\mathrm{R}_{\mathrm{E}}$]', ylabel=r'Number density $n$ [$\mathrm{cm^{-3}}$]', yscale='log', ylim=(1E-4, 1E6))
ax1.set_title(r'(a)', x=-0.05, y=1.05)
if (series_number == 10):
    ax1.plot(length2planet_half_planet_radius, number_density_half[:, 0]*1E-6 + number_density_half[:, 3]*1E-6, c='orange', label=r'$\mathrm{O^+}$(I)', linestyle='dotted', linewidth='4')
    ax1.plot(length2planet_half_planet_radius, number_density_half[:, 1]*1E-6 + number_density_half[:, 4]*1E-6, c='dimgrey', label=r'$\mathrm{H^+}$(I)', linestyle='dotted', linewidth='4')
    ax1.plot(length2planet_half_planet_radius, number_density_half[:, 2]*1E-6 + number_density_half[:, 5]*1E-6, c='blue', label=r'$\mathrm{e^-}$(I)', linestyle='dotted', linewidth='4')
    ax1.plot(length2planet_half_planet_radius, number_density_half[:, 6]*1E-6, c='purple', label=r'$\mathrm{H^+}$(M)', linestyle='dotted', linewidth='4')
    ax1.plot(length2planet_half_planet_radius, number_density_half[:, 7]*1E-6, c='green', label=r'$\mathrm{e^-}$(M)', linestyle='dotted', linewidth='4')
    ax1.plot(length2planet_half_planet_radius, number_density_half[:, 8]*1E-6 + number_density_half[:, 9]*1E-6, c='deepskyblue', label=r'$\mathrm{e^-}$(T)', linestyle='dotted', linewidth='4')
elif (series_number == 12):
    ax1.plot(length2planet_half_planet_radius, number_density_half[:, 0]*1E-6 + number_density_half[:, 4]*1E-6, c='orange', label=r'$\mathrm{O^+}$(I)', linestyle='dotted', linewidth='4')
    ax1.plot(length2planet_half_planet_radius, number_density_half[:, 1]*1E-6 + number_density_half[:, 5]*1E-6, c='dimgrey', label=r'$\mathrm{H^+}$(I)', linestyle='dotted', linewidth='4')
    ax1.plot(length2planet_half_planet_radius, number_density_half[:, 2]*1E-6 + number_density_half[:, 6]*1E-6, c='blue', label=r'cold $\mathrm{e^-}$(I)', linestyle='dotted', linewidth='4')
    ax1.plot(length2planet_half_planet_radius, number_density_half[:, 3]*1E-6 + number_density_half[:, 7]*1E-6, c='hotpink', label=r'hot $\mathrm{e^-}$(I)', linestyle='dotted', linewidth='4')
    ax1.plot(length2planet_half_planet_radius, number_density_half[:, 8]*1E-6, c='purple', label=r'$\mathrm{H^+}$(M)', linestyle='dotted', linewidth='4')
    ax1.plot(length2planet_half_planet_radius, number_density_half[:, 9]*1E-6, c='green', label=r'$\mathrm{e^-}$(M)', linestyle='dotted', linewidth='4')
    ax1.plot(length2planet_half_planet_radius, number_density_half[:, 10]*1E-6 + number_density_half[:, 11]*1E-6, c='deepskyblue', label=r'$\mathrm{e^-}$(T)', linestyle='dotted', linewidth='4')
ax1.minorticks_on()
ax1.grid(which='both', alpha=0.3)
ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)

ax2 = fig.add_subplot(gs[5:8, 0:6], xlabel=r'Geocentric distance [$\mathrm{R}_{\mathrm{E}}$]', ylabel=r'Electrostatic potential $\phi$ [$\mathrm{kV}$]')
ax2.set_title(r'(b)', x=-0.1, y=0.95)
ax2.plot(length2planet_half_planet_radius, electrostatic_potential_half*1E-3, linewidth='4', c='blue', linestyle='solid', label=r'$\phi \, [\mathrm{kV}]$', alpha=0.7)
ax2.minorticks_on()
ax2.grid(axis='both', which='both', alpha=0.3)

ax3 = fig.add_subplot(gs[8:10, 0:6], xlabel=r'Geocentric distance [$\mathrm{R}_{\mathrm{E}}$]', ylabel=r'$\phi$ [$\mathrm{V}$]')
ax3.set_title(r'(c)', x=-0.1, y=0.95)
ax3.plot(length2planet_half_planet_radius, electrostatic_potential_half, linewidth='4', c='blue')
electrostatic_potential_half_min = np.min(electrostatic_potential_half)
ax3.set_ylim(electrostatic_potential_half_min-1, -electrostatic_potential_half_min+1)
ax3.minorticks_on()
ax3.grid(which="both", alpha=0.3)

ax4 = fig.add_subplot(gs[10:, 0:6], xlabel=r'Geocentric distance [$\mathrm{R}_{\mathrm{E}}$]', ylabel=r'Convergence number', yscale='log')
ax4.set_title(r'(d)', x=-0.1, y=0.95)
ax4.plot(length2planet_half_planet_radius, convergence_number_half, linewidth='4', c='blue')
ax4.minorticks_on()
ax4.grid(which="both", alpha=0.3)

plt.tight_layout()
plt.subplots_adjust(wspace=1, hspace=1)
plt.savefig(path_figure)