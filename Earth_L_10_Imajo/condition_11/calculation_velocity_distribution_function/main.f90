program main

    use constant_parameter
    use constant_in_the_simulation
    use reference_results_setting
    use main_variables

    implicit none

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

    !call calculate_zeropoints_and_weights_for_Gaussian(adiabatic_invariant_grid_number, zero_points_mu, weights_mu)
    !call calculate_zeropoints_and_weights_for_Gaussian(parallel_grid_number, zero_points_parallel, weights_parallel)

    !------------------------
    ! 1st adiabatic invariant
    !------------------------

    call make_adiabatic_invariant(particle_mass, magnetic_flux_density, injection_grid_number, adiabatic_invariant)



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
    call make_alim(particle_mass, potential_plus_Bmu, injection_grid_number, amin, alim, boundary_temperature_para)
    call make_amax(particle_mass, potential_plus_Bmu, injection_grid_number, boundary_temperature_para, amin, alim, amax)


    !----------
    ! make file
    !----------
    
    print *, result_dir
    call makedirs(result_dir)
    
    do count_s = 1, boundary_series_number

        injection_grid = injection_grid_number(count_s)

        call make_result_file_name(count_s, result_file_name)
        print *, result_file_name
        open(30, file = result_file_name)
        print *, number_density(count_s, plot_grid_number)

        if (number_density(count_s, plot_grid_number) /= 0d0) then

        do count_mu = 1, adiabatic_invariant_grid_number
            
            v_perp_i = sqrt(2d0 * magnetic_flux_density(plot_grid_number) * adiabatic_invariant(count_s, count_mu) &
                & / particle_mass(count_s))
            v_perp_b = sqrt(2d0 * magnetic_flux_density(injection_grid) * adiabatic_invariant(count_s, count_mu) &
                & / particle_mass(count_s))
            
            !if (count_s == 5) then
            !    print *, amin(count_s, plot_grid_number, count_mu)* sqrt(2d0 / particle_mass(count_s)), &
            !        & alim(count_s, plot_grid_number, count_mu)* sqrt(2d0 / particle_mass(count_s)), &
            !        & amax(count_s, plot_grid_number, count_mu)* sqrt(2d0 / particle_mass(count_s))
            !end if
                
            if ( alim(count_s, plot_grid_number, count_mu) > amin(count_s, plot_grid_number, count_mu) ) then
                if ( amin(count_s, plot_grid_number, count_mu) == 0d0 ) then

                    do count_a = 2, parallel_grid_number
                        
                        a_i = (alim(count_s, plot_grid_number, count_mu) * 1d-8) &
                            & * 1d8**(dble(count_a - 2) / dble(parallel_grid_number - 2))
                        
                        if ( plot_grid_number >= injection_grid ) then
                            sign_a = 1d0

                        else if ( plot_grid_number < injection_grid ) then
                            sign_a = -1d0

                        end if

                        v_para_i = sign_a * sqrt(2d0 / particle_mass(count_s)) * a_i
                        v_para_b = sqrt(a_i**2d0 + potential_plus_Bmu(count_s, plot_grid_number, count_mu) &
                            & - potential_plus_Bmu(count_s, injection_grid, count_mu)) * sqrt(2d0 / particle_mass(count_s)) &
                            & * sign_a

                        distribution_function_i = boundary_number_density(count_s) / number_density(count_s, plot_grid_number) &
                            & / (2d0 * pi * boundary_temperature_perp(count_s) / particle_mass(count_s)) &
                            & / sqrt(2d0 * pi * boundary_temperature_para(count_s) / particle_mass(count_s)) &
                            & * exp(- particle_mass(count_s) / 2d0 &
                            & * (v_para_b**2d0 / boundary_temperature_para(count_s) &
                            & + v_perp_b**2d0 / boundary_temperature_perp(count_s)))

                        differential_flux_i = elementary_charge / particle_mass(count_s) * abs(v_para_i) &
                            & * number_density(count_s, plot_grid_number) &
                            & * sqrt(v_perp_i**2d0 + v_para_i**2d0) * distribution_function_i * 1d-1
    
                        if (distribution_function_i > 1d-50) then
                            write(30, '(1PE25.15E3, 6(",", 1PE25.15E3))') mlat_degree(plot_grid_number), v_perp_i, v_para_i, &
                                & v_perp_b, v_para_b, distribution_function_i, differential_flux_i
                        end if

                    end do  !count_a


                else if ( amin(count_s, plot_grid_number, count_mu) /= 0d0 ) then

                    do count_a = 1, parallel_grid_number

                        a_i = amin(count_s, plot_grid_number, count_mu) &
                            & * (alim(count_s, plot_grid_number, count_mu) / amin(count_s, plot_grid_number, count_mu)) &
                            & **(dble(count_a - 1) / dble(parallel_grid_number - 1))
                        
                        if ( plot_grid_number >= injection_grid ) then
                            sign_a = 1d0

                        else if ( plot_grid_number < injection_grid ) then
                            sign_a = -1d0

                        end if

                        v_para_i = sign_a * sqrt(2d0 / particle_mass(count_s)) * a_i
                        v_para_b = sqrt(a_i**2d0 + potential_plus_Bmu(count_s, plot_grid_number, count_mu) &
                            & - potential_plus_Bmu(count_s, injection_grid, count_mu)) * sqrt(2d0 / particle_mass(count_s)) &
                            & * sign_a

                        distribution_function_i = boundary_number_density(count_s) / number_density(count_s, plot_grid_number) &
                            & / (2d0 * pi * boundary_temperature_perp(count_s) / particle_mass(count_s)) &
                            & / sqrt(2d0 * pi * boundary_temperature_para(count_s) / particle_mass(count_s)) &
                            & * exp(- particle_mass(count_s) / 2d0 &
                            & * (v_para_b**2d0 / boundary_temperature_para(count_s) &
                            & + v_perp_b**2d0 / boundary_temperature_perp(count_s)))

                        differential_flux_i = elementary_charge / particle_mass(count_s) * abs(v_para_i) &
                            & * number_density(count_s, plot_grid_number) &
                            & * sqrt(v_perp_i**2d0 + v_para_i**2d0) * distribution_function_i * 1d-1
    
                        if (distribution_function_i > 1d-50) then
                            write(30, '(1PE25.15E3, 6(",", 1PE25.15E3))') mlat_degree(plot_grid_number), v_perp_i, v_para_i, &
                                & v_perp_b, v_para_b, distribution_function_i, differential_flux_i
                        end if

                    end do  !count_a

                end if

            end if
        
            if ( amax(count_s, plot_grid_number, count_mu) > amin(count_s, plot_grid_number, count_mu) ) then
                if ( amin(count_s, plot_grid_number, count_mu) == 0d0 ) then

                    do count_a = 2, parallel_grid_number

                        a_i = (amax(count_s, plot_grid_number, count_mu)*1d-8) &
                            & * 1d8**(dble(count_a - 2) / dble(adiabatic_invariant_grid_number - 2))
                        
                        if ( plot_grid_number >= injection_grid ) then
                            sign_a = -1d0

                        else if ( plot_grid_number < injection_grid ) then
                            sign_a = 1d0

                        end if

                        v_para_i = sign_a * sqrt(2d0 / particle_mass(count_s)) * a_i
                        v_para_b = sqrt(a_i**2d0 + potential_plus_Bmu(count_s, plot_grid_number, count_mu) &
                            & - potential_plus_Bmu(count_s, injection_grid, count_mu)) * sqrt(2d0 / particle_mass(count_s)) &
                            & * sign_a

                        distribution_function_i = boundary_number_density(count_s) / number_density(count_s, plot_grid_number) &
                            & / (2d0 * pi * boundary_temperature_perp(count_s) / particle_mass(count_s)) &
                            & / sqrt(2d0 * pi * boundary_temperature_para(count_s) / particle_mass(count_s)) &
                            & * exp(- particle_mass(count_s) / 2d0 &
                            & * (v_para_b**2d0 / boundary_temperature_para(count_s) &
                            & + v_perp_b**2d0 / boundary_temperature_perp(count_s)))

                        differential_flux_i = elementary_charge / particle_mass(count_s) * abs(v_para_i) &
                            & * number_density(count_s, plot_grid_number) &
                            & * sqrt(v_perp_i**2d0 + v_para_i**2d0) * distribution_function_i * 1d-1
    
                        if (distribution_function_i > 1d-50) then
                            write(30, '(1PE25.15E3, 6(",", 1PE25.15E3))') mlat_degree(plot_grid_number), v_perp_i, v_para_i, &
                                & v_perp_b, v_para_b, distribution_function_i, differential_flux_i
                        end if

                    end do  !count_a

                else if ( amin(count_s, plot_grid_number, count_mu) /= 0d0 ) then

                    do count_a = 1, parallel_grid_number

                        a_i = amin(count_s, plot_grid_number, count_mu) &
                            & * (amax(count_s, plot_grid_number, count_mu) / amin(count_s, plot_grid_number, count_mu)) &
                            & **(dble(count_a - 1) / dble(parallel_grid_number - 1))
                        
                        if ( plot_grid_number >= injection_grid ) then
                            sign_a = -1d0

                        else if ( plot_grid_number < injection_grid ) then
                            sign_a = 1d0

                        end if

                        v_para_i = sign_a * sqrt(2d0 / particle_mass(count_s)) * a_i
                        v_para_b = sqrt(a_i**2d0 + potential_plus_Bmu(count_s, plot_grid_number, count_mu) &
                            & - potential_plus_Bmu(count_s, injection_grid, count_mu)) * sqrt(2d0 / particle_mass(count_s)) &
                            & * sign_a

                        distribution_function_i = boundary_number_density(count_s) / number_density(count_s, plot_grid_number) &
                            & / (2d0 * pi * boundary_temperature_perp(count_s) / particle_mass(count_s)) &
                            & / sqrt(2d0 * pi * boundary_temperature_para(count_s) / particle_mass(count_s)) &
                            & * exp(- particle_mass(count_s) / 2d0 &
                            & * (v_para_b**2d0 / boundary_temperature_para(count_s) &
                            & + v_perp_b**2d0 / boundary_temperature_perp(count_s)))

                        differential_flux_i = elementary_charge / particle_mass(count_s) * abs(v_para_i) &
                            & * number_density(count_s, plot_grid_number) &
                            & * sqrt(v_perp_i**2d0 + v_para_i**2d0) * distribution_function_i * 1d-1
    
                        if (distribution_function_i > 1d-50) then
                            write(30, '(1PE25.15E3, 6(",", 1PE25.15E3))') mlat_degree(plot_grid_number), v_perp_i, v_para_i, &
                                & v_perp_b, v_para_b, distribution_function_i, differential_flux_i
                        end if

                    end do  !count_a

                end if

            end if

        end do  !count_mu

        end if

        close(30)

    end do  !count_s
    
end program main