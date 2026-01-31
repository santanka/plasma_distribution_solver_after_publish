import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = ['Computer Modern Roman']
mpl.rcParams['mathtext.fontset'] = 'cm'

plt.rcParams["font.size"] = 50

fig = plt.figure(figsize=(22, 14), dpi=100, tight_layout=True)
ax = fig.add_subplot(111, xlabel=r'[$\mathrm{R_{E}}$]', ylabel=r'[$\mathrm{R_{E}}$]', xlim=(-1.5, 9.5), ylim=(0, 3.5))

#Earth plot
Earth_radius = 6.3781E6 #[m]
Earth_equatorial_radius = 6.3781E6 / Earth_radius #[RE]
Earth_polar_radius = 6.3568E6 / Earth_radius #[RE]

mlat_rad_2pi = np.linspace(0E0, 1E0 * np.pi, 10000)

Earth_x = Earth_equatorial_radius * np.cos(mlat_rad_2pi)
Earth_y = Earth_polar_radius * np.sin(mlat_rad_2pi)

#ax.plot(Earth_x, Earth_y, color='darkgoldenrod', linewidth='4')


#altitude 500 km plot
Earth_start_altitude = 500E3 / Earth_radius #[RE]

Earth_x_altitude = (Earth_equatorial_radius + Earth_start_altitude) * np.cos(mlat_rad_2pi)
Earth_y_altitude = (Earth_polar_radius + Earth_start_altitude) * np.sin(mlat_rad_2pi)

#ax.plot(Earth_x_altitude, Earth_y_altitude, color='darkgoldenrod', linewidth='4')
ax.plot(Earth_x_altitude, Earth_y_altitude, color='k', linewidth='4')


#magnetic field line plot
Earth_l_shell = 9E0

a_req_b = Earth_equatorial_radius**2E0 + 2E0 * Earth_l_shell * Earth_start_altitude - np.abs(Earth_polar_radius)**2E0
Earth_mlat_1 = (a_req_b + np.sqrt(a_req_b**2E0 + 4E0 * Earth_l_shell**2E0 * (Earth_polar_radius**2E0 - np.abs(Earth_start_altitude)**2E0))) / 2E0 / Earth_l_shell**2E0
Earth_mlat_1 = - np.arccos(np.sqrt(Earth_mlat_1))
Earth_mlat_2 = - Earth_mlat_1

Earth_mlat_1 = 0

print(Earth_mlat_1*180E0/np.pi)

mlat_rad = np.linspace(Earth_mlat_1, Earth_mlat_2, 10000)

field_line_x = Earth_l_shell * np.cos(mlat_rad)**2E0 * np.cos(mlat_rad)
field_line_y = Earth_l_shell * np.cos(mlat_rad)**2E0 * np.sin(mlat_rad)

#ax.plot(field_line_x, field_line_y, color='purple', linewidth='4')
ax.plot(field_line_x, field_line_y, color='k', linewidth='4')

#potential drop point
#first_drop_R_from_center = 3E0 #[RE]
#second_drop_R_from_center = 8E0 #[RE]
#
#first_drop_point_mlat = np.arccos(np.sqrt(first_drop_R_from_center / Earth_l_shell))
#second_drop_point_mlat = np.arccos(np.sqrt(second_drop_R_from_center / Earth_l_shell))
#
#first_drop_point_x = first_drop_R_from_center * np.cos(first_drop_point_mlat)
#first_drop_point_y = first_drop_R_from_center * np.sin(first_drop_point_mlat)
#
#second_drop_point_x = second_drop_R_from_center * np.cos(second_drop_point_mlat)
#second_drop_point_y = second_drop_R_from_center * np.sin(second_drop_point_mlat)
#
#ax.scatter(first_drop_point_x, first_drop_point_y, marker='D', s=200, c='orangered', zorder=5)
#ax.scatter(first_drop_point_x, -first_drop_point_y, marker='D', s=200, c='orangered', zorder=5)
#ax.scatter(second_drop_point_x, second_drop_point_y, marker='D', s=200, c='orangered', zorder=5)
#ax.scatter(second_drop_point_x, -second_drop_point_y, marker='D', s=200, c='orangered', zorder=5)
#
#
##boundary dot plot
#ax.scatter(field_line_x[0], field_line_y[0], marker='s', s=200, c='orangered', zorder=5)
#ax.scatter(field_line_x[-1], field_line_y[-1], marker='s', s=200, c='orangered', zorder=5)
#ax.scatter(Earth_l_shell, 0, marker='s', s=200, c='orangered', zorder=5)
#
#
##text
#ax.text(0, 0, r"Earth", ha='center', va='center', fontsize=80)
#
#ax.text(1.6, field_line_y[-1] + 0.3, r"North", ha='center', va='center', fontsize=50)
#ax.text(1.6, field_line_y[-1] + 0.0, r"ionospheric", ha='center', va='center', fontsize=50)
#ax.text(1.6, field_line_y[-1] - 0.3, r"end", ha='center', va='center', fontsize=50)
#
#ax.text(1.6, field_line_y[0] + 0.3, r"South", ha='center', va='center', fontsize=50)
#ax.text(1.6, field_line_y[0] + 0.0, r"ionospheric", ha='center', va='center', fontsize=50)
#ax.text(1.6, field_line_y[0] - 0.3, r"end", ha='center', va='center', fontsize=50)
#
#ax.text(10.8, 0.4, r"magnetic", ha='center', va='center', fontsize=50)
#ax.text(10.8, 0.1, r"equator", ha='center', va='center', fontsize=50)
#
#ax.text(0,  1.4, r"$3.5 \, \mathrm{kV}$", ha='center', va='center', fontsize=50)
#ax.text(0, -1.4, r"$3.5 \, \mathrm{kV}$", ha='center', va='center', fontsize=50)
#ax.text(10.5, -0.3, r"$0 \, \mathrm{V}$", ha='center', va='center', fontsize=50)
#
#ax.text(first_drop_point_x-0.6, first_drop_point_y, r"$3 \, \mathrm{R_{E}}$", ha='center', va='center', fontsize=50)
#ax.text(second_drop_point_x+0.6, second_drop_point_y, r"$8 \, \mathrm{R_{E}}$", ha='center', va='center', fontsize=50)
#ax.text(first_drop_point_x-0.6, -first_drop_point_y, r"$3 \, \mathrm{R_{E}}$", ha='center', va='center', fontsize=50)
#ax.text(second_drop_point_x+0.6, -second_drop_point_y, r"$8 \, \mathrm{R_{E}}$", ha='center', va='center', fontsize=50)

ax.minorticks_on()
ax.grid(which='both', alpha=0.3)

plt.axis('scaled')
plt.show()


quit()

#Jupiter plot
Jupiter_radius              = 7.1492E7                      #[m]
Jupiter_equatorial_radius   = 7.1492E7 / Jupiter_radius     #[RJ]
Jupiter_polar_radius        = 6.6854E7 / Jupiter_radius     #[RJ]

mlat_rad_2pi = np.linspace(0E0, 2E0 * np.pi, 10000)

Jupiter_x = Jupiter_equatorial_radius * np.cos(mlat_rad_2pi)
Jupiter_y = Jupiter_polar_radius * np.sin(mlat_rad_2pi)

ax.plot(Jupiter_x, Jupiter_y, color='darkgoldenrod', linewidth='4')


#altitude 2500 km plot
Jupiter_start_altitude = 2.5E6 / Jupiter_radius        #[RJ]

Jupiter_x_altitude = (Jupiter_equatorial_radius + Jupiter_start_altitude) * np.cos(mlat_rad_2pi)
Jupiter_y_altitude = (Jupiter_polar_radius + Jupiter_start_altitude) * np.sin(mlat_rad_2pi)

ax.plot(Jupiter_x_altitude, Jupiter_y_altitude, color='darkgoldenrod', linewidth='4')


#magnetic field line plot
Jupiter_l_shell = 5.91E0

a_req_b = Jupiter_equatorial_radius**2E0 + 2E0 * Jupiter_l_shell * Jupiter_start_altitude - np.abs(Jupiter_polar_radius)**2E0
Jupiter_mlat_1 = (a_req_b + np.sqrt(a_req_b**2E0 + 4E0 * Jupiter_l_shell**2E0 * (Jupiter_polar_radius**2E0 - np.abs(Jupiter_start_altitude)**2E0))) / 2E0 / Jupiter_l_shell**2E0
Jupiter_mlat_1 = - np.arccos(np.sqrt(Jupiter_mlat_1))
Jupiter_mlat_2 = - Jupiter_mlat_1

mlat_rad = np.linspace(Jupiter_mlat_1, Jupiter_mlat_2, 10000)

field_line_x = Jupiter_l_shell * np.cos(mlat_rad)**2E0 * np.cos(mlat_rad)
field_line_y = Jupiter_l_shell * np.cos(mlat_rad)**2E0 * np.sin(mlat_rad)

ax.plot(field_line_x, field_line_y, color='purple', linewidth='4')


#Io plot
Io_rotation_radius  = 4.217E8  / Jupiter_radius     #[RJ]
Io_radius           = 1.8216E6 / Jupiter_radius     #[RJ]

Io_x = Io_rotation_radius + Io_radius * np.cos(mlat_rad_2pi)
Io_y =                  0 + Io_radius * np.sin(mlat_rad_2pi)

ax.plot(Io_x, Io_y, color='darkgoldenrod', linewidth='4')


#boundary dot plot
ax.scatter(field_line_x[0], field_line_y[0], marker='s', s=100, c='orangered', zorder=5)
ax.scatter(field_line_x[-1], field_line_y[-1], marker='s', s=100, c='orangered', zorder=5)
ax.scatter(Jupiter_l_shell, 0, marker='s', s=100, c='orangered', zorder=5)


#text
ax.text(0, 0, r"Jupiter", ha='center', va='center', fontsize=100)

ax.text(1.2, field_line_y[-1] + 0.2, r"North", ha='center', va='center', fontsize=50)
ax.text(1.2, field_line_y[-1] + 0.0, r"ionospheric", ha='center', va='center', fontsize=50)
ax.text(1.2, field_line_y[-1] - 0.2, r"end", ha='center', va='center', fontsize=50)

ax.text(1.2, field_line_y[0] + 0.2, r"South", ha='center', va='center', fontsize=50)
ax.text(1.2, field_line_y[0] + 0.0, r"ionospheric", ha='center', va='center', fontsize=50)
ax.text(1.2, field_line_y[0] - 0.2, r"end", ha='center', va='center', fontsize=50)

ax.text(5.65, 0, r"Io", ha='center', va='center', fontsize=50)

ax.text(6.5, 0.3, r"magnetic", ha='center', va='center', fontsize=50)
ax.text(6.5, 0.1, r"equator", ha='center', va='center', fontsize=50)

ax.text(0,  1.2, r"$30 \, \mathrm{kV}$", ha='center', va='center', fontsize=70)
ax.text(0, -1.2, r"$30 \, \mathrm{kV}$", ha='center', va='center', fontsize=70)
ax.text(6.5, -0.2, r"$0 \, \mathrm{V}$", ha='center', va='center', fontsize=70)

ax.minorticks_on()
ax.grid(which='both', alpha=0.3)

plt.axis('equal')
plt.show()