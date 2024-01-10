program main

    use constant_parameter
    use constant_in_the_simulation
    use main_variables
    !$use omp_lib

    implicit none

    call date_and_time(dummy_1, dummy_2, dummy_3, date_time_start)
    print *, "start time: ", date_time_start(1), "/", date_time_start(2), "/", date_time_start(3), "  ", &
        & date_time_start(5), ":", date_time_start(6), ":", date_time_start(7), ".", date_time_start(8)

    print_omp = .true.
    call set_omp_num_threads(print_omp, omp_num)
    print_omp = .false.

    !-------------------------
    ! initial setting of field
    !-------------------------

    call calculate_planet_mlat(planet_mlat_1, planet_mlat_2)


    !mlat, length2planet, length2satellite, coordinate_FA
    call make_spatial_grid(planet_mlat_1, planet_mlat_2, mlat, length2planet, length2satellite, coordinate_FA)
    
    !diff_coordinate_FA
    do count_i = 1, real_grid_number - 1
        diff_coordinate_FA(count_i) = coordinate_FA(count_i + 1) - coordinate_FA(count_i)
    end do  !count_i

    !magnetic_flux_density
    magnetic_flux_density = magnetic_constant * planet_dipole_moment / 4d0 / pi / (planet_l_shell * planet_radius)**3d0 &
        & * sqrt(1d0 + 3d0 * sin(mlat)**2d0) / cos(mlat)**6d0

    print *, "initial setting of field is finished"

    
    !------------------
    ! initial condition
    !------------------
    
    do count_i = 1, real_grid_number

        if ( (count_i <= initial_grid_ionophere_middle_1 - 1) &
            & .or. (initial_grid_ionophere_middle_2 <= count_i .and. count_i <= real_grid_number) ) then
            initial_electrostatic_potential(count_i) = initial_electrostatic_potential_ionosphere

        else if ( (initial_grid_ionophere_middle_1 <= count_i .and. count_i <= initial_grid_middle_magnetosphere_1 - 1) &
            & .or. (initial_grid_middle_magnetosphere_2 <= count_i .and. count_i <= initial_grid_ionophere_middle_2 - 1) ) then
            
            !if (initial_grid_ionophere_middle_1 <= count_i .and. count_i <= initial_grid_middle_magnetosphere_1 - 1) then
            !    initial_electrostatic_potential(count_i) = initial_electrostatic_potential_middle &
            !        & - dble((count_i - initial_grid_ionophere_middle_1) * 0.1d0)
            !else if (initial_grid_middle_magnetosphere_2 <= count_i .and. count_i <= initial_grid_ionophere_middle_2 - 1) then
            !    initial_electrostatic_potential(count_i) = initial_electrostatic_potential_middle &
            !        & - dble((initial_grid_ionophere_middle_2 - 1 - count_i) * 0.1d0)
            !end if
            initial_electrostatic_potential(count_i) = initial_electrostatic_potential_middle
        
        else if ( initial_grid_middle_magnetosphere_1 <= count_i .and. count_i <= initial_grid_middle_magnetosphere_2 - 1 ) then
            initial_electrostatic_potential(count_i) = initial_electrostatic_potential_magnetosphere
        
        end if

    end do  !count_i

    electrostatic_potential = initial_electrostatic_potential

    print *, "initial condition is finished"


    !--------------------
    ! boundary conditions
    !--------------------

    call make_input_file_name(boundary_file)
    open(30, file = boundary_file, action = 'read', status = 'old')
    read(30, *)
    do count_s = 1, boundary_series_number

        read(30, *) dummy, initial_boundary_number_density(count_s), boundary_temperature_perp(count_s), &
            & boundary_temperature_para(count_s), charge_number(count_s), particle_mass(count_s), injection_grid_number(count_s)
    
    end do  !count_s
    close(30)

    charge_number = charge_number * elementary_charge
    boundary_temperature_perp = boundary_temperature_perp * elementary_charge
    boundary_temperature_para = boundary_temperature_para * elementary_charge
    boundary_number_density = initial_boundary_number_density

    print *, "boundary conditions are finished"


    !-----------------------------
    ! Gaussian-Legendre quadrature
    !-----------------------------

    call calculate_zeropoints_and_weights_for_Gaussian(adiabatic_invariant_grid_number, zero_points_mu, weights_mu)
    call calculate_zeropoints_and_weights_for_Gaussian(parallel_grid_number, zero_points_parallel, weights_parallel)

    print *, "Gaussian-Legendre quadrature is finished"


    !------------------------
    ! 1st adiabatic invariant
    !------------------------

    call make_adiabatic_invariant(boundary_temperature_perp, magnetic_flux_density, injection_grid_number, zero_points_mu, &
        & adiabatic_invariant)

    print *, "1st adiabatic invariant is finished"


    !----------------
    ! iteration start
    !----------------

    print *, "iteration start"
    

    do count_min_position = initial_min_grid_start, initial_min_grid_end

        initial_min_grid_1 = count_min_position
        initial_min_grid_2 = real_grid_number + 1 - initial_min_grid_1

        count_iteration = 0
        count_iteration_timer = 0
        convergence_number_sum_min = 1d4
        electrostatic_potential = initial_electrostatic_potential
        boundary_number_density = initial_boundary_number_density

        iteration : do  !count_iteration
            count_iteration = count_iteration + 1

            !-----------
            ! make _diff
            !-----------

            call make_electrostatic_potential_diff(electrostatic_potential, initial_min_grid_1, initial_min_grid_2, &
                & electrostatic_potential_diff)

            call make_boundary_number_density_diff(boundary_number_density, boundary_number_density_diff)

            call make_potential_energy_diff(mlat, length2planet, length2satellite, charge_number, particle_mass, &
                & electrostatic_potential_diff, potential_energy_diff)

            call make_potential_plus_Bmu_diff(potential_energy_diff, adiabatic_invariant, magnetic_flux_density, &
                & potential_plus_Bmu_diff)


            !--------------
            ! accessibility
            !--------------

            call make_amin(potential_plus_Bmu_diff, injection_grid_number, amin)

            call make_alim(potential_plus_Bmu_diff, injection_grid_number, amin, alim, boundary_temperature_para)

            call make_amax(potential_plus_Bmu_diff, injection_grid_number, boundary_temperature_para, amin, alim, amax)


            !--------
            ! density
            !--------

            call make_number_density(boundary_number_density_diff, boundary_temperature_perp, boundary_temperature_para, &
            & potential_energy_diff, magnetic_flux_density, adiabatic_invariant, injection_grid_number, amin, alim, amax, &
            & zero_points_parallel, weights_mu, weights_parallel, number_density_diff)

            call cannot_reach_check(number_density_diff, injection_grid_number)

            call make_charge_density_from_number_density(number_density_diff, charge_number, charge_density_diff, &
                & charge_density_plus_diff, charge_density_minus_diff)

            call make_charge_density_from_Poisson_eq(electrostatic_potential_diff, diff_coordinate_FA, charge_density_poisson_diff)


            !------------------
            ! convergence check
            !------------------

            call make_convergence_number(charge_density_diff, charge_density_plus_diff, charge_density_minus_diff, &
                & charge_density_poisson_diff, convergence_number_diff, convergence_number_sum)

            if ( convergence_number_sum_min > convergence_number_sum .or. convergence_number_sum /= convergence_number_sum ) then

                call make_result_file_name(initial_min_grid_1, result_file)
                call make_result_file_format(format_character)
                open(50, file = result_file)
                write(50, "(1PE25.15E3)") convergence_number_sum
                do count_i = 1, real_grid_number

                    write(50, format_character) coordinate_FA(count_i), length2planet(count_i), mlat(count_i), &
                        & mlat(count_i) * deg_per_rad, magnetic_flux_density(count_i), initial_electrostatic_potential(count_i), &
                        & electrostatic_potential_diff(1, count_i), number_density_diff(1, :, count_i), &
                        & charge_density_diff(1, count_i), charge_density_poisson_diff(1, count_i), &
                        & convergence_number_diff(1, count_i)

                end do
                do count_s = 1, boundary_series_number

                    write(50, "(1PE25.15E3, 4(',', 1PE25.15E3), ',', I25)") boundary_number_density_diff(1, count_s), &
                        & boundary_temperature_perp(count_s) / elementary_charge, &
                        & boundary_temperature_para(count_s) / elementary_charge, charge_number(count_s) / elementary_charge, &
                        & particle_mass(count_s), injection_grid_number(count_s)

                end do
                close(50)
            end if

            if ( convergence_number_sum_min > convergence_number_sum &
                & .and. (convergence_number_sum_min - convergence_number_sum) > convergence_number_sum_min * 1d-3) then
                convergence_number_sum_min = convergence_number_sum
                count_iteration_timer = count_iteration
            end if

            if ( count_iteration - count_iteration_timer == 5E3 ) then
                print *, "finish(time up)"
                exit iteration
            end if

            if ( convergence_number_sum /= convergence_number_sum ) then
                print *, "finish(error: NaN)"
                exit iteration
            end if


            !------------------------------
            ! update by using Newton method
            !------------------------------

            call Newton_method_for_electrostatic_potential(electrostatic_potential_diff, convergence_number_diff, &
                & initial_min_grid_1, initial_min_grid_2, electrostatic_potential)

            call Newton_method_for_boundary_number_density(boundary_number_density_diff, convergence_number_diff, &
                & injection_grid_number, boundary_number_density)


            !------
            ! print
            !------

            if ( mod(count_iteration, 50) == 1 ) then
                call set_omp_num_threads(print_omp, omp_num)
                call date_and_time(dummy_1, dummy_2, dummy_3, date_time)
                print *, count_min_position, count_iteration, count_iteration_timer, convergence_number_sum_min, &
                & convergence_number_sum, omp_num, &
                & electrostatic_potential(initial_min_grid_1), minval(electrostatic_potential), minloc(electrostatic_potential), &
                & maxval(electrostatic_potential), maxloc(electrostatic_potential), &
                & date_time(1), '/', date_time(2), '/', date_time(3), "  ", date_time(5), ":", date_time(6), &
                & ":", date_time(7), ".", date_time(8)
            end if

        end do iteration !count_iteration

    end do  !count_min_position

    !各ファイルのconvergence_number_sum_minを一つのファイルにまとめる
    call make_convergence_number_file_name(convergence_number_file)
    open(60, file = convergence_number_file, action = 'write')
    do count_min_position = initial_min_grid_start, initial_min_grid_end

        initial_min_grid_1 = count_min_position
        initial_min_grid_2 = real_grid_number + 1 - initial_min_grid_1

        call make_result_file_name(initial_min_grid_1, result_file)

        inquire(file = result_file, exist = file_exist)
        if ( file_exist .eqv. .false. ) then
            cycle
        end if

        open(50, file = result_file, action = 'read', status = 'old')
        read(50, *) convergence_number_sum
        close(50)

        write(60, "(I3, 5(1PE25.15E3), ',')") count_min_position, coordinate_FA(initial_min_grid_2), &
            & length2planet(initial_min_grid_2), mlat(initial_min_grid_2), mlat(initial_min_grid_2) * deg_per_rad, &
            & convergence_number_sum

    end do !count_min_position
    close(60)

    call date_and_time(dummy_1, dummy_2, dummy_3, date_time_end)
    print *, "finish"
    print *, "start time: ", date_time_start(1), "/", date_time_start(2), "/", date_time_start(3), "  ", &
        & date_time_start(5), ":", date_time_start(6), ":", date_time_start(7), ".", date_time_start(8)
    print *, "end time: ", date_time_end(1), "/", date_time_end(2), "/", date_time_end(3), "  ", &
        & date_time_end(5), ":", date_time_end(6), ":", date_time_end(7), ".", date_time_end(8)

end program main