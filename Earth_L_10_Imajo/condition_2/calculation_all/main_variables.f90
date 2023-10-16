module main_variables
    use constant_parameter
    use constant_in_the_simulation
    use reference_results_setting, only: boundary_series_number

    implicit none
    
    !-----------------
    ! reference result
    !-----------------

    double precision, dimension(real_grid_number) :: mlat_rad, mlat_degree, length2planet, coordinate_FA
    double precision, dimension(real_grid_number) :: magnetic_flux_density, initial_electrostatic_potential
    double precision, dimension(real_grid_number) :: electrostatic_potential
    double precision, dimension(boundary_series_number, real_grid_number) :: number_density
    double precision, dimension(real_grid_number) :: charge_density, charge_density_poisson
    double precision, dimension(real_grid_number) :: convergence_number
    double precision :: convergence_number_sum

    double precision, dimension(boundary_series_number) :: boundary_number_density, boundary_temperature_perp
    double precision, dimension(boundary_series_number) :: boundary_temperature_para, charge_number, particle_mass
    integer, dimension(boundary_series_number) :: injection_grid_number


    !----------
    ! satellite
    !----------

    double precision, dimension(real_grid_number) :: length2satellite

    !-----------------------------
    ! Gaussian-Legendre quadrature
    !-----------------------------

    double precision, dimension(adiabatic_invariant_grid_number) :: zero_points_mu, weights_mu
    double precision, dimension(parallel_grid_number) :: zero_points_parallel, weights_parallel


    !----------
    ! variables
    !----------

    double precision, dimension(boundary_series_number, adiabatic_invariant_grid_number) :: adiabatic_invariant
    double precision, dimension(boundary_series_number, real_grid_number) :: potential_energy
    double precision, dimension(boundary_series_number, real_grid_number, adiabatic_invariant_grid_number) :: amin, alim, amax
    double precision, dimension(boundary_series_number, real_grid_number, adiabatic_invariant_grid_number) :: potential_plus_Bmu
    double precision, dimension(boundary_series_number, real_grid_number) :: particle_flux_density, pressure_perp, pressure_para
    double precision, dimension(boundary_series_number, real_grid_number) :: temperature_perp, temperature_para
    double precision, dimension(boundary_series_number, real_grid_number) :: parallel_mean_velocity, pressure_dynamic
    double precision, dimension(real_grid_number) :: Alfven_speed, Alfven_speed_per_lightspeed
    double precision, dimension(real_grid_number) :: electron_inertial_length, ion_inertial_length
    double precision, dimension(real_grid_number) :: ion_Larmor_radius, ion_acoustic_radius, electron_Larmor_radius
    double precision, dimension(real_grid_number) :: current_density


    !--------
    ! counter
    !--------

    integer :: count_s, count_i


    !------------
    ! file format
    !------------

    character(len = 180) :: format_character

    
end module main_variables