import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os

grid_ionosphere_middle = 14
grid_middle_magnetosphere = 109
grid_fix = 175

BC_number = 1
min_number = 175

channel = 1

l_shell = 10E0
series_number = 8

mass_series = np.array([2.677950266103E-26, 1.67262192369E-27, 9.1093837015E-31, 9.1093837015E-31, 2.677950266103E-26, 1.67262192369E-27, 9.1093837015E-31, 9.1093837015E-31, 1.67262192369E-27, 9.1093837015E-31])
#mass_series = np.array([2.677950266103E-26, 1.67262192369E-27, 9.1093837015E-31, 2.677950266103E-26, 1.67262192369E-27, 9.1093837015E-31, 1.67262192369E-27, 9.1093837015E-31])
charge_series = np.array([1E0, 1E0, -1E0, -1E0, 1E0, 1E0, -1E0, -1E0, 1E0, -1E0])
#charge_series = np.array([1E0, 1E0, -1E0, 1E0, 1E0, -1E0, 1E0, -1E0])

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
planet_rotation = 2E0 * np.pi / 86164.09053083288E0   #[rad s-1]

length2planet_half_planet_radius = length2planet_half / planet_radius

mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = ['Computer Modern Roman']
mpl.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams["font.size"] = 20

B_ratio = magnetic_flux_density_half / magnetic_flux_density_half[-1]
electron_Fluid_number_density_boundary = 2E11   #[m-3]
electron_Fluid_temperature_boundary = 1E-1 * elementary_charge #[J]
electron_Fluid_number_density = electron_Fluid_number_density_boundary * np.exp(elementary_charge * (electrostatic_potential_half - electrostatic_potential_half[-1]) / electron_Fluid_temperature_boundary) \
                                * B_ratio

electron_Fluid_temperature_boundary_2 = 5E-1 * elementary_charge #[J]
electron_Fluid_number_density_2 = electron_Fluid_number_density_boundary * np.exp(elementary_charge * (electrostatic_potential_half - electrostatic_potential_half[-1]) / electron_Fluid_temperature_boundary_2) \
                                * B_ratio

Oxygen_Fluid_number_density_boundary = 2E11   #[m-3]
Oxygen_Fluid_temperature_boundary = 1E-1 * elementary_charge #[J]
Oxygen_Fluid_number_density = Oxygen_Fluid_number_density_boundary * np.exp(2E0 * constant_of_gravitation * planet_mass * mass_series[0] / planet_radius / l_shell / Oxygen_Fluid_temperature_boundary \
                                                                            * ((1E0 / np.cos(mlat_rad_half) - 1E0 / np.cos(mlat_rad_half[-1])) \
                                                                               #+ 2.5E-1 * planet_rotation**2E0 * (planet_radius * l_shell)**3E0 / constant_of_gravitation / planet_mass * (np.cos(mlat_rad_half)**6E0 - np.cos(mlat_rad_half[-1])**6E0)) \
                                                                                - elementary_charge * (electrostatic_potential_half - electrostatic_potential_half[-1]) / Oxygen_Fluid_temperature_boundary)) * B_ratio

Oxygen_Fluid_temperature_boundary_2 = 5E-1 * elementary_charge #[J]
Oxygen_Fluid_number_density_2 = Oxygen_Fluid_number_density_boundary * np.exp(2E0 * constant_of_gravitation * planet_mass * mass_series[0] / planet_radius / l_shell / Oxygen_Fluid_temperature_boundary_2 \
                                                                            * ((1E0 / np.cos(mlat_rad_half) - 1E0 / np.cos(mlat_rad_half[-1])) \
                                                                                 #+ 2.5E-1 * planet_rotation**2E0 * (planet_radius * l_shell)**3E0 / constant_of_gravitation / planet_mass * (np.cos(mlat_rad_half)**6E0 - np.cos(mlat_rad_half[-1])**6E0)) \
                                                                                    - elementary_charge * (electrostatic_potential_half - electrostatic_potential_half[-1]) / Oxygen_Fluid_temperature_boundary_2)) * B_ratio


proton_Fluid_number_density_boundary = 2E8   #[m-3]
proton_Fluid_temperature_boundary = 1E-1 * elementary_charge #[J]
proton_Fluid_number_density = proton_Fluid_number_density_boundary * np.exp(2E0 * constant_of_gravitation * planet_mass * mass_series[1] / planet_radius / l_shell / proton_Fluid_temperature_boundary \
                                                                            * ((1E0 / np.cos(mlat_rad_half) - 1E0 / np.cos(mlat_rad_half[-1])) \
                                                                                #+ 2.5E-1 * planet_rotation**2E0 * (planet_radius * l_shell)**3E0 / constant_of_gravitation / planet_mass * (np.cos(mlat_rad_half)**6E0 - np.cos(mlat_rad_half[-1])**6E0)) \
                                                                                - elementary_charge * (electrostatic_potential_half - electrostatic_potential_half[-1]) / proton_Fluid_temperature_boundary)) * B_ratio

proton_Fluid_temperature_boundary_2 = 1E0 * elementary_charge #[J]
proton_Fluid_number_density_2 = proton_Fluid_number_density_boundary * np.exp(2E0 * constant_of_gravitation * planet_mass * mass_series[1] / planet_radius / l_shell / proton_Fluid_temperature_boundary_2 \
                                                                            * ((1E0 / np.cos(mlat_rad_half) - 1E0 / np.cos(mlat_rad_half[-1])) \
                                                                                #+ 2.5E-1 * planet_rotation**2E0 * (planet_radius * l_shell)**3E0 / constant_of_gravitation / planet_mass * (np.cos(mlat_rad_half)**6E0 - np.cos(mlat_rad_half[-1])**6E0)) \
                                                                                - elementary_charge * (electrostatic_potential_half - electrostatic_potential_half[-1]) / proton_Fluid_temperature_boundary_2)) * B_ratio


#electron_Fluid_number_densityとnumber_density_halfの比較

fig = plt.figure(figsize=(14, 14), dpi=100, tight_layout=True)
#ax = fig.add_subplot(111, xlabel=r'Geocentric distance [$\mathrm{R}_{\mathrm{E}}$]', ylabel=r'Number density $n$ [$\mathrm{cm^{-3}}$]', yscale='log', ylim=(1E-4, 1E6), xlim=(500, 10000))
ax = fig.add_subplot(111, xlabel=r'Altitude [$\mathrm{km}$]', ylabel=r'Number density $n$ [$\mathrm{cm^{-3}}$]', yscale='log', ylim=(1E-4, 1E6), xlim=(500, 10000))

ax.plot((length2planet_half - planet_radius) * 1E-3, electron_Fluid_number_density*1E-6, label='e- Boltzman Fluid (Ergun, 0.1 eV)', color='blue', linestyle='solid', linewidth=4, alpha=0.5)
ax.plot((length2planet_half - planet_radius) * 1E-3, electron_Fluid_number_density_2*1E-6, label='e- Boltzman Fluid (Ergun, 0.5 eV)', color='blue', linestyle='dashed', linewidth=4, alpha=0.5)
ax.plot((length2planet_half - planet_radius) * 1E-3, number_density_half[2, :]*1E-6 + number_density_half[5, :]*1E-6, label='e- (PDS, 0.5 eV)', color='blue', linestyle='dotted', linewidth=4, alpha=0.5)

ax.plot((length2planet_half - planet_radius) * 1E-3, Oxygen_Fluid_number_density*1E-6, label='O+ Fluid (Ergun, 0.1 eV)', color='orange', linestyle='solid', linewidth=4, alpha=0.5)
ax.plot((length2planet_half - planet_radius) * 1E-3, Oxygen_Fluid_number_density_2*1E-6, label='O+ Fluid (Ergun, 0.5 eV)', color='orange', linestyle='dashed', linewidth=4, alpha=0.5)
ax.plot((length2planet_half - planet_radius) * 1E-3, number_density_half[0, :]*1E-6 + number_density_half[3, :]*1E-6, label='O+ (PDS, 0.5 eV)', color='orange', linestyle='dotted', linewidth=4, alpha=0.5)

ax.plot((length2planet_half - planet_radius) * 1E-3, proton_Fluid_number_density*1E-6, label='H+ Fluid (Ergun, 0.1 eV)', color='green', linestyle='solid', linewidth=4, alpha=0.5)
ax.plot((length2planet_half - planet_radius) * 1E-3, proton_Fluid_number_density_2*1E-6, label='H+ Fluid (Ergun, 1.0 eV)', color='green', linestyle='dashed', linewidth=4, alpha=0.5)
ax.plot((length2planet_half - planet_radius) * 1E-3, number_density_half[1, :]*1E-6 + number_density_half[4, :]*1E-6, label='H+ (PDS, 1.0 eV)', color='green', linestyle='dotted', linewidth=4, alpha=0.5)


ax.minorticks_on()
ax.grid(which='both', alpha=0.3)
ax.legend(loc='best')

plt.show()