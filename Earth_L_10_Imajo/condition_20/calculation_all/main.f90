program main

    use constant_parameter
    use constant_in_the_simulation
    use reference_results_setting
    use main_variables
    !$use omp_lib

    implicit none

    call set_omp_num_threads(.true.)

    !-------------
    ! data loading
    !-------------

    open(20, file = reference_file, action = 'read', status = 'old')
    read(20, *) convergence_number_sum
    do count_i = 1, real_grid_number
        
        read(20, *) coordinate_FA(count_i), length2planet(count_i), mlat_rad(count_i), mlat_degree(count_i), &
            & magnetic_flux_density(count_i), initial_electrostatic_potential(count_i), electrostatic_potential(count_i), &
            & number_density(:, count_i), charge_density(count_i), charge_density_poisson(count_i), convergence_number(count_i)

    end do  !count_i

    do count_s = 1, boundary_series_number
        
        read(20, *) boundary_number_density(count_s), boundary_temperature_perp(count_s), boundary_temperature_para(count_s), &
            & charge_number(count_s), particle_mass(count_s), injection_grid_number(count_s)

    end do  !count_s
    close(20)

    charge_number = charge_number * elementary_charge
    boundary_temperature_perp = boundary_temperature_perp * elementary_charge
    boundary_temperature_para = boundary_temperature_para * elementary_charge


    !----------
    ! satellite
    !----------

    if ( satellite_mass /= 0d0 ) then
        length2satellite =  sqrt(length2planet**2d0 + (planet_radius * satellite_l_shell)**2d0 &
            & - 2d0 * length2planet * planet_radius * satellite_l_shell * cos(mlat_rad))
    end if

    !-----------------------------
    ! Gaussian-Legendre quadrature
    !-----------------------------

    call calculate_zeropoints_and_weights_for_Gaussian(adiabatic_invariant_grid_number, zero_points_mu, weights_mu)
    call calculate_zeropoints_and_weights_for_Gaussian(parallel_grid_number, zero_points_parallel, weights_parallel)


    !------------------------
    ! 1st adiabatic invariant
    !------------------------

    call make_adiabatic_invariant(boundary_temperature_perp, magnetic_flux_density, injection_grid_number, zero_points_mu, &
        & adiabatic_invariant)


    !----------
    ! potential
    !----------

    call make_potential_energy(mlat_rad, length2planet, length2satellite, charge_number, particle_mass, electrostatic_potential, &
        & potential_energy)

    call make_potential_plus_Bmu(potential_energy, adiabatic_invariant, magnetic_flux_density, potential_plus_Bmu)


    !--------------
    ! accessibility
    !--------------

    call make_amin(potential_plus_Bmu, injection_grid_number, amin)
    call make_alim(potential_plus_Bmu, injection_grid_number, amin, alim, boundary_temperature_para)
    call make_amax(potential_plus_Bmu, injection_grid_number, boundary_temperature_para, amin, alim, amax)

    print *, 'make accessibility'


    !------------
    ! calculation
    !------------

    call make_number_density(boundary_number_density, boundary_temperature_perp, boundary_temperature_para, &
        & potential_energy, magnetic_flux_density, adiabatic_invariant, injection_grid_number, amin, alim, amax, &
        & zero_points_parallel, weights_mu, weights_parallel, number_density)

    call cannot_reach_check(number_density, injection_grid_number)

    print *, 'make number density'

    call make_particle_flux_density(boundary_number_density, magnetic_flux_density, adiabatic_invariant, injection_grid_number, &
        & boundary_temperature_perp, boundary_temperature_para, particle_mass, potential_energy, alim, amax, &
        & zero_points_parallel, weights_mu, weights_parallel, number_density, particle_flux_density, parallel_mean_velocity)

    print *, 'make particle flux density'

    call make_pressure_perpendicular(boundary_number_density, boundary_temperature_perp, boundary_temperature_para, &
        & potential_energy, injection_grid_number, magnetic_flux_density, adiabatic_invariant, number_density, amin, alim, amax, &
        & zero_points_parallel, weights_mu, weights_parallel, pressure_perp, temperature_perp)

    print *, 'make pressure perpendicular'
    
    call make_pressure_parallel(boundary_number_density, boundary_temperature_perp, boundary_temperature_para, &
        & injection_grid_number, magnetic_flux_density, adiabatic_invariant, potential_energy, particle_mass, zero_points_parallel,&
        & weights_mu, weights_parallel, number_density, parallel_mean_velocity, amin, alim, amax, pressure_para, temperature_para)

    print *, 'make pressure parallel'

    call make_pressure_dynamic(number_density, parallel_mean_velocity, particle_mass, pressure_dynamic)

    call make_Alfven_speed(magnetic_flux_density, particle_mass, number_density, Alfven_speed, Alfven_speed_per_lightspeed)

    call make_inertial_length(charge_number, particle_mass, number_density, ion_inertial_length, electron_inertial_length)

    call make_Larmor_radius(magnetic_flux_density, number_density, pressure_perp, charge_number, particle_mass, &
        & ion_Larmor_radius, ion_acoustic_radius, electron_Larmor_radius)

    call make_current_density(charge_number, particle_flux_density, current_density)
    


    !------------
    ! data output
    !------------

    call makedirs(result_dir)
    call make_result_file_format(format_character)
    print *, format_character
    open(30, file = result_file)
    do count_i = 1, real_grid_number
        
        write(30, format_character) coordinate_FA(count_i), length2planet(count_i), mlat_rad(count_i), mlat_degree(count_i), &
            & magnetic_flux_density(count_i), initial_electrostatic_potential(count_i), electrostatic_potential(count_i), &
            & number_density(:, count_i), charge_density(count_i), charge_density_poisson(count_i), convergence_number(count_i), &
            & particle_flux_density(:, count_i), parallel_mean_velocity(:, count_i), pressure_perp(:, count_i), &
            & pressure_para(:, count_i), pressure_dynamic(:, count_i), temperature_perp(:, count_i), temperature_para(:, count_i), &
            & Alfven_speed(count_i), Alfven_speed_per_lightspeed(count_i), ion_inertial_length(count_i), &
            & electron_inertial_length(count_i), ion_Larmor_radius(count_i), ion_acoustic_radius(count_i), &
            & electron_Larmor_radius(count_i), current_density(count_i)

    end do  !count_i
    close(30)
    
end program main