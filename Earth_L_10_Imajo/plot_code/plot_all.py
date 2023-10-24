import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os

grid_ionosphere_middle = 14
grid_middle_magnetosphere = 109
grid_fix = 175

BC_number = 3
min_number = 163

channel = 1

l_shell = 10E0
series_number = 8

#mass_series = np.array([2.677950266103E-26, 1.67262192369E-27, 9.1093837015E-31, 9.1093837015E-31, 2.677950266103E-26, 1.67262192369E-27, 9.1093837015E-31, 9.1093837015E-31, 1.67262192369E-27, 9.1093837015E-31])
mass_series = np.array([2.677950266103E-26, 1.67262192369E-27, 9.1093837015E-31, 2.677950266103E-26, 1.67262192369E-27, 9.1093837015E-31, 1.67262192369E-27, 9.1093837015E-31])
#charge_series = np.array([1E0, 1E0, -1E0, -1E0, 1E0, 1E0, -1E0, -1E0, 1E0, -1E0])
charge_series = np.array([1E0, 1E0, -1E0, 1E0, 1E0, -1E0, 1E0, -1E0])

dir_name = f'/mnt/j/plasma_distribution_solver_after_publish/Earth_L_10_Imajo/alpha_perp_12_parallel_12/grid_{str(grid_ionosphere_middle).zfill(3)}_{str(grid_middle_magnetosphere).zfill(3)}_{str(grid_fix).zfill(3)}/'
dir_BC_name = f'boundary_condition_{str(BC_number)}/'
file_name = f'all/all_{str(min_number).zfill(3)}.csv'

path_dir = f'{dir_name}{dir_BC_name}plot/'

path_name = f'{dir_name}{dir_BC_name}{file_name}'

print(path_name)

data = np.genfromtxt(path_name, delimiter=',', unpack=True)

coordinate_FA = data[0, :]
length2planet = data[1, :]
mlat_rad = data[2, :]
mlat_degree = data[3, :]
magnetic_flux_density = data[4, :]
initial_electrostatic_potential = data[5, :]
electrostatic_potential = data[6, :]
number_density = data[7:series_number + 7, :]
charge_density = data[series_number + 7, :]
charge_density_Poisson = data[series_number + 8, :]
convergence_number = data[series_number + 9, :]
particle_flux_density = data[series_number + 10 : 2 * series_number + 10, :]
parallel_mean_velocity = data[2 * series_number + 10 : 3 * series_number + 10, :]
pressure_perp = data[3 * series_number + 10 : 4 * series_number + 10, :]
pressure_para = data[4 * series_number + 10 : 5 * series_number + 10, :]
pressure_dynamic = data[5 * series_number + 10 : 6 * series_number + 10, :]
temperature_perp = data[6 * series_number + 10 : 7 * series_number + 10, :]
temperature_para = data[7 * series_number + 10 : 8 * series_number + 10, :]
Alfven_speed = data[8 * series_number + 10, :]
Alfven_speed_per_lightspeed = data[8 * series_number + 11, :]
ion_inertial_length = data[8 * series_number + 12, :]
electron_inertial_length = data[8 * series_number + 13, :]
ion_Larmor_radius = data[8 * series_number + 14, :]
ion_acoustic_gyroradius = data[8 * series_number + 15, :]
electron_Larmor_radius = data[8 * series_number + 16, :]
current_density = data[8 * series_number + 17, :]

data_length = len(mlat_degree)
data_length_half = int((1+data_length)/2)

#data half
coordinate_FA_half                      = data[0, data_length_half-1:]
length2planet_half                      = data[1, data_length_half-1:]
mlat_rad_half                           = data[2, data_length_half-1:]
mlat_degree_half                        = data[3, data_length_half-1:]
magnetic_flux_density_half              = data[4, data_length_half-1:]
initial_electrostatic_potential_half    = data[5, data_length_half-1:]
electrostatic_potential_half            = data[6, data_length_half-1:]
number_density_half                     = data[7:series_number + 7, data_length_half-1:]
charge_density_half                     = data[series_number + 7, data_length_half-1:]
charge_density_Poisson_half             = data[series_number + 8, data_length_half-1:]
convergence_number_half                 = data[series_number + 9, data_length_half-1:]
particle_flux_density_half              = data[series_number + 10 : 2 * series_number + 10, data_length_half-1:]
parallel_mean_velocity_half             = data[2 * series_number + 10 : 3 * series_number + 10, data_length_half-1:]
pressure_perp_half                      = data[3 * series_number + 10 : 4 * series_number + 10, data_length_half-1:]
pressure_para_half                      = data[4 * series_number + 10 : 5 * series_number + 10, data_length_half-1:]
pressure_dynamic_half                   = data[5 * series_number + 10 : 6 * series_number + 10, data_length_half-1:]
temperature_perp_half                   = data[6 * series_number + 10 : 7 * series_number + 10, data_length_half-1:]
temperature_para_half                   = data[7 * series_number + 10 : 8 * series_number + 10, data_length_half-1:]
Alfven_speed_half                       = data[8 * series_number + 10, data_length_half-1:]
Alfven_speed_per_lightspeed_half        = data[8 * series_number + 11, data_length_half-1:]
ion_inertial_length_half                = data[8 * series_number + 12, data_length_half-1:]
electron_inertial_length_half           = data[8 * series_number + 13, data_length_half-1:]
ion_Larmor_radius_half                  = data[8 * series_number + 14, data_length_half-1:]
ion_acoustic_gyroradius_half            = data[8 * series_number + 15, data_length_half-1:]
electron_Larmor_radius_half             = data[8 * series_number + 16, data_length_half-1:]
current_density_half                    = data[8 * series_number + 17, data_length_half-1:]


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

if (channel == 1):

    electron_perp_pressure_half = np.zeros(data_length_half)
    ion_perp_pressure_half      = np.zeros(data_length_half)

    electron_para_pressure_half = np.zeros(data_length_half)
    ion_para_pressure_half      = np.zeros(data_length_half)

    for count_i in range(series_number):
        if charge_series[count_i] < 0:
            electron_perp_pressure_half += pressure_perp_half[count_i, :]
            electron_para_pressure_half += pressure_para_half[count_i, :]
        else:
            ion_perp_pressure_half += pressure_perp_half[count_i, :]
            ion_para_pressure_half += pressure_para_half[count_i, :]
    
    all_perp_pressure_half      = electron_perp_pressure_half + ion_perp_pressure_half

    all_para_pressure_half      = electron_para_pressure_half + ion_para_pressure_half

    electron_total_pressure_half = (2E0 * electron_perp_pressure_half + electron_para_pressure_half) / 3E0
    ion_total_pressure_half      = (2E0 * ion_perp_pressure_half + ion_para_pressure_half) / 3E0
    all_total_pressure_half      = electron_total_pressure_half + ion_total_pressure_half

    electron_cyclotron_freq_MHz = elementary_charge * magnetic_flux_density_half / electron_mass / 2E0 / np.pi * 1E-6

    fig = plt.figure(figsize=(14, 28), dpi=100, tight_layout=True)
    gs = fig.add_gridspec(16, 7)

    ax1 = fig.add_subplot(gs[0:5, 0:6], xlabel=r'Geocentric distance [$\mathrm{R}_{\mathrm{E}}$]', ylabel=r'Number density $n$ [$\mathrm{cm^{-3}}$]', yscale='log', ylim=(1E-4, 1E6))
    ax1.set_title(r'(a)', x=-0.05, y=1.05)
    ax1.plot(length2planet_half_planet_radius, number_density_half[0, :]*1E-6 + number_density_half[3, :]*1E-6, c='orange', label=r'$\mathrm{O^+}$(I)', linestyle='dotted', linewidth='4')
    ax1.plot(length2planet_half_planet_radius, number_density_half[1, :]*1E-6 + number_density_half[4, :]*1E-6, c='dimgrey', label=r'$\mathrm{H^+}$(I)', linestyle='dotted', linewidth='4')
    ax1.plot(length2planet_half_planet_radius, number_density_half[2, :]*1E-6 + number_density_half[5, :]*1E-6, c='blue', label=r'$\mathrm{e^-}$(I)', linestyle='dotted', linewidth='4')
    ax1.plot(length2planet_half_planet_radius, number_density_half[6, :]*1E-6, c='purple', label=r'$\mathrm{H^+}$(M)', linestyle='dotted', linewidth='4')
    ax1.plot(length2planet_half_planet_radius, number_density_half[7, :]*1E-6, c='green', label=r'$\mathrm{e^-}$(M)', linestyle='dotted', linewidth='4')
    #ax1.plot(length2planet_half_planet_radius, number_density_half[0, :]*1E-6 + number_density_half[4, :]*1E-6, c='orange', label=r'$\mathrm{O^+}$(I)', linestyle='dotted', linewidth='4')
    #ax1.plot(length2planet_half_planet_radius, number_density_half[1, :]*1E-6 + number_density_half[5, :]*1E-6, c='dimgrey', label=r'$\mathrm{H^+}$(I)', linestyle='dotted', linewidth='4')
    #ax1.plot(length2planet_half_planet_radius, number_density_half[2, :]*1E-6 + number_density_half[6, :]*1E-6, c='blue', label=r'cold $\mathrm{e^-}$(I)', linestyle='dotted', linewidth='4')
    #ax1.plot(length2planet_half_planet_radius, number_density_half[3, :]*1E-6 + number_density_half[7, :]*1E-6, c='hotpink', label=r'hot $\mathrm{e^-}$(I)', linestyle='dotted', linewidth='4')
    #ax1.plot(length2planet_half_planet_radius, number_density_half[8, :]*1E-6, c='purple', label=r'$\mathrm{H^+}$(M)', linestyle='dotted', linewidth='4')
    #ax1.plot(length2planet_half_planet_radius, number_density_half[9, :]*1E-6, c='green', label=r'$\mathrm{e^-}$(M)', linestyle='dotted', linewidth='4')
    ax1.minorticks_on()
    ax1.grid(which='both', alpha=0.3)
    ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)


    ax2 = fig.add_subplot(gs[5:8, 0:], xlabel=r'Geocentric distance [$\mathrm{R}_{\mathrm{E}}$]', ylabel=r'Electrostatic potential $\phi$ [$\mathrm{kV}$]')
    ax2.set_title(r'(b)', x=-0.1, y=0.95)
    ax2_twinx = ax2.twinx()
    ax2_twinx.set_ylabel(r'Electron cyclotron frequency $\Omega_{e}$ [$\mathrm{MHz}$]')
    ax2.plot(length2planet_half_planet_radius, electrostatic_potential_half*1E-3, linewidth='4', c='blue', linestyle='solid', label=r'$\phi \, [\mathrm{kV}]$', alpha=0.7)
    ax2_twinx.plot(length2planet_half_planet_radius, electron_cyclotron_freq_MHz, linewidth='4', c=r'orange', linestyle='dotted', label=r'$\Omega_{\mathrm{e}} \, [\mathrm{MHz}]$')
    ax2.minorticks_on()
    ax2.grid(axis='both', which='both', alpha=0.3)
    ax2_twinx.minorticks_on()
    #ax2_twinx.grid(which='both', axis='y', alpha=0.3)
    h1, l1 = ax2.get_legend_handles_labels()
    h2, l2 = ax2_twinx.get_legend_handles_labels()
    ax2_twinx.legend(h1+h2, l1+l2)

    ax3 = fig.add_subplot(gs[8:10, 0:], xlabel=r'Geocentric distance [$\mathrm{R}_{\mathrm{E}}$]', ylabel=r'$\phi$ [$\mathrm{V}$]')
    ax3.set_title(r'(c)', x=-0.1, y=0.95)
    ax3.plot(length2planet_half_planet_radius, electrostatic_potential_half, linewidth='4', c='blue')
    ax3.set_ylim(-1, 27)
    ax3.minorticks_on()
    ax3.grid(which="both", alpha=0.3)

    ax4 = fig.add_subplot(gs[11:, 0:2], yscale='log', xlabel=r'Geocentric distance [$\mathrm{R}_{\mathrm{E}}$]', ylabel=r'Pressure [$\mathrm{nPa}$]')
    ax4.set_title(r'(d) Perpendicular $P_{\perp}$')
    ax4.plot(length2planet_half_planet_radius, all_perp_pressure_half*1E9, c='purple', label='all', linewidth='4', alpha=0.7)
    ax4.plot(length2planet_half_planet_radius, ion_perp_pressure_half*1E9, c='orange', label='ion', linewidth='4', alpha=0.7)
    ax4.plot(length2planet_half_planet_radius, electron_perp_pressure_half*1E9, c='blue', label='electron', linewidth='4', alpha=0.7)
    ax4.minorticks_on()
    ax4.grid(which="both", alpha=0.3)

    ax5 = fig.add_subplot(gs[11:, 2:4], yscale='log', xlabel=r'Geocentric distance [$\mathrm{R}_{\mathrm{E}}$]')
    ax5.set_title(r'(e) Parallel $P_{\parallel}$')
    ax5.plot(length2planet_half_planet_radius, all_para_pressure_half*1E9, c='purple', label='all', linewidth='4', alpha=0.7)
    ax5.plot(length2planet_half_planet_radius, ion_para_pressure_half*1E9, c='orange', label='ion', linewidth='4', alpha=0.7)
    ax5.plot(length2planet_half_planet_radius, electron_para_pressure_half*1E9, c='blue', label='electron', linewidth='4', alpha=0.7)
    ax5.minorticks_on()
    ax5.grid(which="both", alpha=0.3)

    ax6 = fig.add_subplot(gs[11:, 4:6], yscale='log', xlabel=r'Geocentric distance [$\mathrm{R}_{\mathrm{E}}$]')
    ax6.set_title(r'(f) Total $P$')
    ax6.plot(length2planet_half_planet_radius, all_total_pressure_half*1E9, c='purple', label='all', linewidth='4', alpha=0.7)
    ax6.plot(length2planet_half_planet_radius, ion_total_pressure_half*1E9, c='orange', label='ion', linewidth='4', alpha=0.7)
    ax6.plot(length2planet_half_planet_radius, electron_total_pressure_half*1E9, c='blue', label='electron', linewidth='4', alpha=0.7)
    ax6.minorticks_on()
    ax6.grid(which="both", alpha=0.3)

    ax6.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0, fontsize=20)

    plt.tight_layout()
    plt.subplots_adjust(wspace=1, hspace=1)

    
    #フォルダがなければ作成
    if not os.path.isdir(path_dir):
        os.makedirs(path_dir)

    plt.savefig(f'{path_dir}numberdensity_electrostaticpotential_pressure_BC{str(BC_number)}_min{str(min_number).zfill(3)}.png')
    plt.savefig(f'{path_dir}numberdensity_electrostaticpotential_pressure_BC{str(BC_number)}_min{str(min_number).zfill(3)}.pdf')