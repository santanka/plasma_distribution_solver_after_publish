subroutine set_omp_num_threads(print_omp, omp_num)
    use omp_lib

    implicit none
    
    logical, intent(in) :: print_omp
    integer, intent(out) :: omp_num

    ! Get the number of available processors
    omp_num = min(omp_get_num_procs(), omp_get_max_threads())
    if ( print_omp .eqv. .true. ) print *, "Maximum threads: ", omp_num

    ! 最大スレッド数で並列処理を実行する
    call omp_set_num_threads(omp_num)

end subroutine set_omp_num_threads
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine calculate_planet_mlat(planet_mlat_1, planet_mlat_2)
    use boundary_and_initial_conditions

    implicit none
    
    double precision, intent(out) :: planet_mlat_1, planet_mlat_2

    double precision :: a_req_b

    a_req_b = planet_equatorial_radius**2d0 + 2d0 * planet_l_shell * planet_radius * planet_start_altitude &
        & - abs(planet_polar_radius)**2d0
    planet_mlat_1 = (a_req_b + sqrt(a_req_b**2d0 + 4d0 * (planet_l_shell * planet_radius)**2d0 &
        & * (planet_polar_radius**2d0 - abs(planet_start_altitude)**2d0))) / 2d0 / (planet_l_shell * planet_radius)**2d0
    planet_mlat_1 = - acos(sqrt(planet_mlat_1))
    planet_mlat_2 = - planet_mlat_1
    
end subroutine calculate_planet_mlat
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_spatial_grid(planet_mlat_1, planet_mlat_2, mlat, length2planet, length2satellite, coordinate_FA)
    use constant_parameter
    use constant_in_the_simulation

    implicit none
    
    double precision, intent(in) :: planet_mlat_1, planet_mlat_2
    double precision, dimension(real_grid_number), intent(out) :: mlat, length2planet, length2satellite, coordinate_FA

    double precision :: planet_coordinate_FA_1, planet_coordinate_FA_2, r_eq, function1, function2, mlat_before, mlat_after
    integer :: count_i

    r_eq = planet_equatorial_radius * planet_l_shell
    planet_coordinate_FA_1 = r_eq * (asinh(sqrt(3d0) * sin(planet_mlat_1)) / 2d0 / sqrt(3d0) &
        & + sin(planet_mlat_1) * sqrt(5d0 - 3d0 * cos(2d0 * planet_mlat_1)) / 2d0 / sqrt(2d0))
    planet_coordinate_FA_2 = r_eq * (asinh(sqrt(3d0) * sin(planet_mlat_2)) / 2d0 / sqrt(3d0) &
        & + sin(planet_mlat_2) * sqrt(5d0 - 3d0 * cos(2d0 * planet_mlat_2)) / 2d0 / sqrt(2d0))

    do count_i = 1, real_grid_number
        coordinate_FA(count_i) = planet_coordinate_FA_1 + (planet_coordinate_FA_2 - planet_coordinate_FA_1) &
            & * dble(count_i - 1) / dble(real_grid_number - 1)
        
        ! Using Newton method, calculate the magnetic latitude corresponding to the coordinate_FA
        if (count_i == 1) then
            mlat(count_i) = planet_mlat_1
        else if (count_i == real_grid_number) then
            mlat(count_i) = planet_mlat_2
        else
            mlat_after = mlat(count_i - 1)
            function1 = r_eq * (asinh(sqrt(3d0) * sin(mlat_after)) / 2d0 / sqrt(3d0) &
                & + sin(mlat_after) * sqrt(5d0 - 3d0 * cos(2d0 * mlat_after)) / 2d0 / sqrt(2d0)) - coordinate_FA(count_i)
            do while ( abs(function1) > 1d-7 )
                mlat_before = mlat_after
                function2 = r_eq * cos(mlat_before) * sqrt(1d0 + 3d0 * sin(mlat_before)**2d0)
                mlat_after = mlat_before - function1 / function2
                function1 = r_eq * (asinh(sqrt(3d0) * sin(mlat_after)) / 2d0 / sqrt(3d0) &
                    & + sin(mlat_after) * sqrt(5d0 - 3d0 * cos(2d0 * mlat_after)) / 2d0 / sqrt(2d0)) - coordinate_FA(count_i)
            end do
            mlat(count_i) = mlat_after
        end if
    end do  !count_i

    length2planet = r_eq * cos(mlat)**2d0

    if ( satellite_mass /= 0d0 ) then
        length2satellite = sqrt(length2planet**2d0 + r_eq**2d0 - 2d0 * length2planet * r_eq * cos(mlat))
    end if
    
end subroutine make_spatial_grid
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_input_file_name(boundary_file)
    use boundary_and_initial_conditions

    implicit none
    
    character(len = 25), intent(out) ::  boundary_file

    boundary_file = boundary_file_front
    
end subroutine make_input_file_name
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine calculate_zeropoints_and_weights_for_Gaussian(order, zero_points, weights)
    ! based on GSLGIT (https://doi.org/10.1585/jspf1958.64.397)
    use constant_parameter
    !$use omp_lib

    implicit none
    
    integer, intent(in) :: order
    double precision, dimension(order), intent(out) :: zero_points
    double precision, dimension(order), intent(out) :: weights

    double precision :: order_double_precision
    double precision :: Legendre_0, Legendre_1, Legendre_2, Legendre_derivation_2, Legendre_derivation_1
    double precision :: zero_point_update, update
    integer :: count_i, count_j, count_k
    integer :: order_half

    order_double_precision = dble(order)
    order_half = order / 2

    !$omp parallel private(zero_point_update, Legendre_0, Legendre_1, Legendre_2, Legendre_derivation_2, Legendre_derivation_1, &
    !$omp & update, count_j, count_k)
    !$omp do
    do count_i = 1, order_half
        
        zero_point_update = cos(pi / 2d0 / order_double_precision * dble(2 * count_i - 1))

        iteration_end : do count_k = 1, 10000000
        
            Legendre_0 = 1d0
            Legendre_1 = zero_point_update
            
            do count_j = 2, order
            
                Legendre_2 = (dble(2 * count_j - 1) * zero_point_update * Legendre_1 - dble(count_j - 1) * Legendre_0)/dble(count_j)
                Legendre_0 = Legendre_1
                Legendre_1 = Legendre_2
    
            end do  !count_j

            Legendre_derivation_1 = (Legendre_0 - zero_point_update * Legendre_1) &
                & / (1d0 + zero_point_update) / (1d0 - zero_point_update) * order_double_precision

            Legendre_derivation_2 = (2d0 * zero_point_update * Legendre_derivation_1 &
                & - order_double_precision * (order_double_precision + 1d0) * Legendre_1) &
                & / (1d0 + zero_point_update) / (1d0 - zero_point_update)

            update = (1d0 + Legendre_1 * Legendre_derivation_2 / 2d0 / Legendre_derivation_1**2d0) *Legendre_1/Legendre_derivation_1

            zero_point_update = zero_point_update - update

            if ( abs(update) < 1d-16 ) then

                exit iteration_end
                
            end if

            if ( count_k == 10000000 ) then
                print *, "Error!: solution is not found. count_i = ", count_i
                stop
            end if
            
        end do iteration_end  !count_k

        Legendre_0 = 1d0
        Legendre_1 = zero_point_update
            
        do count_j = 2, order
    
            Legendre_2 = (dble(2 * count_j - 1) * zero_point_update * Legendre_1 - dble(count_j - 1) * Legendre_0) / dble(count_j)
            Legendre_0 = Legendre_1
            Legendre_1 = Legendre_2

        end do  !count_j

        zero_points(count_i) = zero_point_update

        weights(count_i) = 2d0 * (1d0 + zero_point_update) * (1d0 - zero_point_update) / (order_double_precision * Legendre_0)**2d0

        !print *, count_i, zero_points(count_i), weights(count_i), Legendre_1
        
        zero_points(order + 1 - count_i) = zero_points(count_i)
        zero_points(count_i) = - zero_points(count_i)

        weights(order + 1 - count_i) = weights(count_i)

    end do  !count_i
    !$omp end do
    !$omp end parallel

    if ( order_half == (order-1)/2 ) then
        
        zero_points((order+1)/2) = 0d0

        Legendre_0 = 1d0
        Legendre_1 = zero_point_update
            
        do count_j = 2, order
    
            Legendre_2 = (dble(2 * count_j - 1) * 0d0 * Legendre_1 - dble(count_j - 1) * Legendre_0) / dble(count_j)
            Legendre_0 = Legendre_1
            Legendre_1 = Legendre_2

        end do  !count_j

        weights((order+1)/2) = 2d0 * (1d0 + 0d0) * (1d0 - 0d0) / (order_double_precision * Legendre_0)**2d0

        !print *, (order+1)/2, zero_points((order+1)/2), weights((order+1)/2), Legendre_1

    end if

    
end subroutine calculate_zeropoints_and_weights_for_Gaussian
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_adiabatic_invariant(boundary_temperature_perp, magnetic_flux_density, injection_grid_number, zero_points_mu, &
    & adiabatic_invariant)
    use constant_parameter
    use constant_in_the_simulation
    use boundary_and_initial_conditions
    !$use omp_lib

    implicit none

    double precision, dimension(boundary_series_number), intent(in) :: boundary_temperature_perp
    double precision, dimension(real_grid_number), intent(in) :: magnetic_flux_density
    integer, dimension(boundary_series_number), intent(in) :: injection_grid_number
    double precision, dimension(adiabatic_invariant_grid_number), intent(in) :: zero_points_mu
    double precision, dimension(boundary_series_number, adiabatic_invariant_grid_number), intent(out) :: adiabatic_invariant

    integer :: count_s, count_mu
    double precision :: max_mu

    do count_s = 1, boundary_series_number

        max_mu = alpha_perp * boundary_temperature_perp(count_s) / magnetic_flux_density(injection_grid_number(count_s))

        !$omp parallel
        !$omp do
        do count_mu = 1, adiabatic_invariant_grid_number

            adiabatic_invariant(count_s, count_mu) = max_mu / 2d0 * (zero_points_mu(count_mu) + 1d0)

        end do  !count_mu
        !$omp end do
        !$omp end parallel

    end do  !count_s
    
end subroutine make_adiabatic_invariant
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_electrostatic_potential_diff(electrostatic_potential, initial_min_grid_1, initial_min_grid_2, &
    & electrostatic_potential_diff)
    use constant_parameter
    use constant_in_the_simulation
    use boundary_and_initial_conditions

    implicit none
    
    double precision, dimension(real_grid_number), intent(in) :: electrostatic_potential
    integer, intent(in) :: initial_min_grid_1, initial_min_grid_2
    double precision, dimension(3, real_grid_number), intent(out) ::  electrostatic_potential_diff

    electrostatic_potential_diff(1, :) = electrostatic_potential
    electrostatic_potential_diff(2, :) = electrostatic_potential + 1d-6
    electrostatic_potential_diff(3, :) = electrostatic_potential - 1d-6
    electrostatic_potential_diff(2, 1) = electrostatic_potential(1)
    electrostatic_potential_diff(2, real_grid_number) = electrostatic_potential(real_grid_number)
    electrostatic_potential_diff(3, 1) = electrostatic_potential(1)
    electrostatic_potential_diff(3, real_grid_number) = electrostatic_potential(real_grid_number)
    if ( initial_fix_grid > initial_min_grid_1 .and. initial_fix_grid < initial_min_grid_2 ) then
        electrostatic_potential_diff(2, initial_fix_grid) = electrostatic_potential(initial_fix_grid)
        electrostatic_potential_diff(3, initial_fix_grid) = electrostatic_potential(initial_fix_grid)
    end if
        
end subroutine make_electrostatic_potential_diff
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_boundary_number_density_diff(boundary_number_density,  boundary_number_density_diff)
    use constant_parameter
    use constant_in_the_simulation
    use boundary_and_initial_conditions

    implicit none
    
    double precision, dimension(boundary_series_number), intent(in) :: boundary_number_density
    double precision, dimension(3, boundary_series_number), intent(out) ::  boundary_number_density_diff

    boundary_number_density_diff(1, :) = boundary_number_density
    boundary_number_density_diff(2, :) = boundary_number_density
    boundary_number_density_diff(3, :) = boundary_number_density

    boundary_number_density_diff(2, boundary_ionosphere_1_variable_species) &
        & = boundary_number_density(boundary_ionosphere_1_variable_species) * (1d0 + 1d-7)
    boundary_number_density_diff(3, boundary_ionosphere_1_variable_species) &
        & = boundary_number_density(boundary_ionosphere_1_variable_species) * (1d0 - 1d-7)
    
    boundary_number_density_diff(2, boundary_ionosphere_2_variable_species) &
        & = boundary_number_density(boundary_ionosphere_2_variable_species) * (1d0 + 1d-7)
    boundary_number_density_diff(3, boundary_ionosphere_2_variable_species) &
        & = boundary_number_density(boundary_ionosphere_2_variable_species) * (1d0 - 1d-7)

    if ( 1 <= boundary_magnetosphere_variable_species .and. boundary_magnetosphere_variable_species <= boundary_series_number ) then
        boundary_number_density_diff(2, boundary_magnetosphere_variable_species) &
            & = boundary_number_density(boundary_magnetosphere_variable_species) * (1d0 + 1d-7)
        boundary_number_density_diff(3, boundary_magnetosphere_variable_species) &
            & = boundary_number_density(boundary_magnetosphere_variable_species) * (1d0 - 1d-7)
    end if

end subroutine make_boundary_number_density_diff
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_potential_energy_diff(mlat, length2planet, length2satellite, charge_number, particle_mass, &
    & electrostatic_potential_diff, potential_energy_diff)
    use constant_parameter
    use constant_in_the_simulation
    use boundary_and_initial_conditions
    !$use omp_lib

    implicit none
    
    double precision, dimension(real_grid_number), intent(in) :: mlat, length2planet, length2satellite
    double precision, dimension(boundary_series_number), intent(in) :: charge_number, particle_mass
    double precision, dimension(3, real_grid_number), intent(in) :: electrostatic_potential_diff
    double precision, dimension(3, boundary_series_number, real_grid_number), intent(out) ::  potential_energy_diff

    integer :: count_h, count_i

    do count_h = 1, 3
        
        !$omp parallel
        !$omp do
        do count_i = 1, real_grid_number

            !gravity of planet
            potential_energy_diff(count_h, :, count_i) = - constant_of_gravitation * planet_mass * particle_mass &
                & / length2planet(count_i)
            
            !centrifugal force of planet
            potential_energy_diff(count_h, :, count_i) = potential_energy_diff(count_h, :, count_i) &
                & - particle_mass * (planet_rotation * length2planet(count_i) * cos(mlat(count_i)))**2d0 / 2d0
            
            !gravity of satellite
            if ( satellite_mass /= 0d0 ) then
                potential_energy_diff(count_h, :, count_i) = potential_energy_diff(count_h, :, count_i) &
                    & - constant_of_gravitation * satellite_mass * particle_mass / length2satellite(count_i)
            end if

            !Coulomb force
            potential_energy_diff(count_h, :, count_i) = potential_energy_diff(count_h, :, count_i) &
                & + charge_number * electrostatic_potential_diff(count_h, count_i)
        
        end do  !count_i
        !$omp end do
        !$omp end parallel
    
    end do  !count_h
    
end subroutine make_potential_energy_diff
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_potential_plus_Bmu_diff(potential_energy_diff, adiabatic_invariant, magnetic_flux_density, &
    & potential_plus_Bmu_diff)
    use constant_parameter
    use constant_in_the_simulation
    use boundary_and_initial_conditions
    !$use omp_lib

    implicit none
    
    double precision, dimension(3, boundary_series_number, real_grid_number), intent(in) :: potential_energy_diff
    double precision, dimension(boundary_series_number, adiabatic_invariant_grid_number), intent(in) :: adiabatic_invariant
    double precision, dimension(real_grid_number), intent(in) ::  magnetic_flux_density
    double precision, dimension(3, boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(out) :: &
        & potential_plus_Bmu_diff

    integer :: count_h, count_s, count_i

    do count_h = 1, 3
        
        do count_s = 1, boundary_series_number

            !$omp parallel
            !$omp do
            do count_i = 1, real_grid_number
                
                potential_plus_Bmu_diff(count_h, count_s, count_i, :) = potential_energy_diff(count_h, count_s, count_i) &
                    & + magnetic_flux_density(count_i) * adiabatic_invariant(count_s, :)

            end do  !count_i
            !$omp end do
            !$omp end parallel
            
        end do  !count_s

    end do  !count_h

end subroutine make_potential_plus_Bmu_diff
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_amin(potential_plus_Bmu_diff, injection_grid_number, amin)
    use constant_parameter
    use constant_in_the_simulation
    use boundary_and_initial_conditions
    !$use omp_lib

    implicit none
    
    double precision, dimension(3, boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(in) :: &
        & potential_plus_Bmu_diff
    integer, dimension(boundary_series_number), intent(in) :: injection_grid_number
    double precision, dimension(3, boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(out) :: amin

    integer :: count_h, count_s, count_i, count_mu, Emax_grid, count4max

    do count_s = 1, boundary_series_number

        !$omp parallel private(Emax_grid, count_h, count_mu, count4max)
        !$omp do
        do count_i = 1, real_grid_number

            if ( count_i /= injection_grid_number(count_s) ) then
                do count_h = 1, 3
                    
                    do count_mu = 1, adiabatic_invariant_grid_number
                        
                        if ( count_i < injection_grid_number(count_s) ) then
                            Emax_grid = count_i + 1
                            do count4max = count_i + 1, injection_grid_number(count_s)

                                if ( potential_plus_Bmu_diff(1, count_s, count4max, count_mu) &
                                    & > potential_plus_Bmu_diff(1, count_s, Emax_grid, count_mu) ) then
                                    Emax_grid = count4max
                                end if
                            
                            end do  !count4max
                        
                        else if ( count_i > injection_grid_number(count_s) ) then
                            Emax_grid = injection_grid_number(count_s)
                            do count4max = injection_grid_number(count_s), count_i - 1

                                if ( potential_plus_Bmu_diff(1, count_s, count4max, count_mu) &
                                    & > potential_plus_Bmu_diff(1, count_s, Emax_grid, count_mu) ) then
                                Emax_grid = count4max
                                end if
                            
                            end do  !count4max
                        end if

                        if ( potential_plus_Bmu_diff(count_h, count_s, count_i, count_mu) &
                            & > potential_plus_Bmu_diff(1, count_s, Emax_grid, count_mu) ) then
                            amin(count_h, count_s, count_i, count_mu) = 0d0

                        else
                            amin(count_h, count_s, count_i, count_mu) &
                                & = sqrt(potential_plus_Bmu_diff(1, count_s, Emax_grid, count_mu) &
                                & - potential_plus_Bmu_diff(count_h, count_s, count_i, count_mu))

                        end if

                    end do  !count_mu

                end do  !count_h

            else if ( count_i == injection_grid_number(count_s) ) then
                amin(:, count_s, count_i, :) = 0d0

            end if
        end do  !count_i
        !$omp end do
        !$omp end parallel
        
    end do  !count_s

    
end subroutine make_amin
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_alim(potential_plus_Bmu_diff, injection_grid_number, amin, alim, boundary_temperature_para)
    use constant_parameter
    use constant_in_the_simulation
    use boundary_and_initial_conditions
    !$use omp_lib

    implicit none
    
    double precision, dimension(3, boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(in) :: &
        & potential_plus_Bmu_diff
    integer, dimension(boundary_series_number), intent(in) :: injection_grid_number
    double precision, dimension(boundary_series_number), intent(in) :: boundary_temperature_para
    double precision, dimension(3, boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(in) :: amin
    double precision, dimension(3, boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(out) :: alim

    integer :: count_h, count_s, count_mu, count_i

    do count_h = 1, 3
        
        do count_s = 1, boundary_series_number

            !$omp parallel private(count_i)
            !$omp do
            do count_mu = 1, adiabatic_invariant_grid_number
                
                alim(count_h, count_s, :, count_mu) = sqrt( boundary_temperature_para(count_s) * alpha_parallel &
                    & + potential_plus_Bmu_diff(1, count_s, injection_grid_number(count_s), count_mu) &
                    & - potential_plus_Bmu_diff(count_h, count_s, :, count_mu))

                do count_i = 1, real_grid_number

                    if ( alim(count_h, count_s, count_i, count_mu) <= amin(count_h, count_s, count_i, count_mu) ) then
                        alim(count_h, count_s, count_i, count_mu) = amin(count_h, count_s, count_i, count_mu)
                    end if
                
                end do  !count_i

            end do  !count_mu
            !$omp end do
            !$omp end parallel

        end do  !count_s

    end do  !count_h
    
end subroutine make_alim
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_amax(potential_plus_Bmu_diff, injection_grid_number, boundary_temperature_para, amin, alim, amax)
    use constant_parameter
    use constant_in_the_simulation
    use boundary_and_initial_conditions
    !$use omp_lib

    implicit none
    
    double precision, dimension(3, boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(in) :: &
        & potential_plus_Bmu_diff
    integer, dimension(boundary_series_number), intent(in) :: injection_grid_number
    double precision, dimension(boundary_series_number), intent(in) :: boundary_temperature_para
    double precision, dimension(3, boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(in) :: amin
    double precision, dimension(3, boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(in) :: alim
    double precision, dimension(3, boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(out) :: amax

    integer :: count_h, count_s, count_i, count_mu, Emax_grid, count4max
    double precision, dimension(3) :: energy_difference, energy_difference_boundary

    do count_s = 1, boundary_series_number

        !$omp parallel private(count4max, Emax_grid, energy_difference, energy_difference_boundary, count_h, count_mu)
        !$omp do
        do count_i = 1, real_grid_number

            if ( (count_i == 1 .and. injection_grid_number(count_s) /= 1) &
                & .or. (count_i == real_grid_number .and. injection_grid_number(count_s) /= real_grid_number) ) then
                amax(:, count_s, count_i, :) = amin(:, count_s, count_i, :)

            else if ( count_i == injection_grid_number(count_s) .and. injection_grid_number(count_s) /= 1 &
                & .and. injection_grid_number(count_s) /= real_grid_number ) then
                    amax(:, count_s, count_i, :) = alim(:, count_s, count_i, :)
            
            else
                do count_mu = 1, adiabatic_invariant_grid_number

                    if ( count_i < injection_grid_number(count_s) .or. injection_grid_number(count_s) == real_grid_number ) then
                        Emax_grid = 1
                        do count4max = 1, count_i - 1

                            if ( potential_plus_Bmu_diff(1, count_s, count4max, count_mu) &
                                & > potential_plus_Bmu_diff(1, count_s, Emax_grid, count_mu) ) then
                                Emax_grid = count4max
                            end if

                        end do  !count4max
                    
                    else if (count_i > injection_grid_number(count_s) .or. injection_grid_number(count_s) == 1) then
                        Emax_grid = count_i + 1
                        do count4max = count_i + 1, real_grid_number
                            
                            if ( potential_plus_Bmu_diff(1, count_s, count4max, count_mu) &
                            & > potential_plus_Bmu_diff(1, count_s, Emax_grid, count_mu) ) then
                                Emax_grid = count4max
                            end if

                        end do  !count4max
                    end if

                    energy_difference = potential_plus_Bmu_diff(1, count_s, Emax_grid, count_mu) &
                        & - potential_plus_Bmu_diff(:, count_s, count_i, count_mu)
                    
                    energy_difference_boundary = potential_plus_Bmu_diff(1, count_s, injection_grid_number(count_s), count_mu) &
                    & - potential_plus_Bmu_diff(:, count_s, count_i, count_mu) + boundary_temperature_para(count_s) * alpha_parallel

                    do count_h = 1, 3
                        
                        if ( energy_difference(count_h) <= amin(count_h, count_s, count_i, count_mu)**2d0 .or. &
                            & energy_difference_boundary(count_h) <= amin(count_h, count_s, count_i, count_mu)**2d0 ) then
                            amax(count_h, count_s, count_i, count_mu) = amin(count_h, count_s, count_i, count_mu)

                        else if ( amin(count_h, count_s, count_i, count_mu)**2d0 < energy_difference(count_h) .and. &
                            & energy_difference(count_h) <= energy_difference_boundary(count_h) ) then
                            amax(count_h, count_s, count_i, count_mu) = sqrt(energy_difference(count_h))
                        
                        else if ( amin(count_h, count_s, count_i, count_mu)**2d0 < energy_difference(count_h) .and. &
                            & energy_difference(count_h) > energy_difference_boundary(count_h) ) then
                            amax(count_h, count_s, count_i, count_mu) = sqrt(energy_difference_boundary(count_h))

                        end if

                        !if ( energy_difference(count_h) <= amin(count_h, count_s, count_i, count_mu)**2d0 ) then
                        !    amax(count_h, count_s, count_i, count_mu) = amin(count_h, count_s, count_i, count_mu)
!
                        !else if ( amin(count_h, count_s, count_i, count_mu)**2d0 < energy_difference(count_h) ) then
                        !    amax(count_h, count_s, count_i, count_mu) = sqrt(energy_difference(count_h))
!
                        !end if

                    end do  !count_h
                    
                end do  !count_mu

            end if
            
        end do  !count_i
        !$omp end do
        !$omp end parallel

    end do  !count_s

    
end subroutine make_amax
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_number_density(boundary_number_density_diff, boundary_temperature_perp, boundary_temperature_para, &
    & potential_energy_diff, magnetic_flux_density, adiabatic_invariant, injection_grid_number, amin, alim, amax, &
    & zero_points_parallel, weights_mu, weights_parallel, number_density_diff)
    use constant_parameter
    use constant_in_the_simulation
    use boundary_and_initial_conditions

    implicit none
    
    double precision, dimension(3, boundary_series_number), intent(in) :: boundary_number_density_diff
    double precision, dimension(boundary_series_number), intent(in) :: boundary_temperature_perp, boundary_temperature_para
    double precision, dimension(3, boundary_series_number, real_grid_number), intent(in) :: potential_energy_diff
    double precision, dimension(real_grid_number), intent(in) :: magnetic_flux_density
    double precision, dimension(boundary_series_number, adiabatic_invariant_grid_number), intent(in) :: adiabatic_invariant
    integer, dimension(boundary_series_number), intent(in) :: injection_grid_number
    double precision, dimension(3, boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(in) :: &
        & amin, alim, amax
    double precision, dimension(parallel_grid_number), intent(in) :: zero_points_parallel, weights_parallel
    double precision, dimension(adiabatic_invariant_grid_number), intent(in) :: weights_mu
    double precision, dimension(3, boundary_series_number, real_grid_number), intent(out) :: number_density_diff

    integer :: count_h, count_s, count_i
    double precision :: integral, sqrt_temperature_para, coefficient4integral
    double precision, dimension(adiabatic_invariant_grid_number) :: y_variable
    double precision, dimension(adiabatic_invariant_grid_number) :: xmin, xlim, xmax

    do count_h = 1, 3
        
        do count_s = 1, boundary_series_number

            do count_i = 1, real_grid_number
                
                integral = 0d0
                if (count_i /= injection_grid_number(count_s) .or. count_h == 1) then
                    coefficient4integral = - (potential_energy_diff(count_h, count_s, count_i) &
                        & - potential_energy_diff(1, count_s, injection_grid_number(count_s))) / boundary_temperature_para(count_s)
                else if (count_i == injection_grid_number(count_s) .and. count_h /= 1) then
                    coefficient4integral = 0
                end if
                
                y_variable = (magnetic_flux_density(injection_grid_number(count_s)) / boundary_temperature_perp(count_s) &
                    & + (magnetic_flux_density(count_i) - magnetic_flux_density(injection_grid_number(count_s))) &
                    & / boundary_temperature_para(count_s)) * adiabatic_invariant(count_s, :)

                sqrt_temperature_para = sqrt(boundary_temperature_para(count_s))
                
                xmin = amin(count_h, count_s, count_i, :) / sqrt_temperature_para
                xlim = alim(count_h, count_s, count_i, :) / sqrt_temperature_para
                xmax = amax(count_h, count_s, count_i, :) / sqrt_temperature_para

                call calculation_x_mu_integral(xmin, xlim, xmax, y_variable, coefficient4integral, zero_points_parallel, &
                    & weights_mu, weights_parallel, integral)

                !number_density_diff(count_h, count_s, count_i) = integral / 2d0 * alpha &
                !    & * (1d0 + boundary_temperature_perp(count_s) / boundary_temperature_para(count_s) &
                !    & * (magnetic_flux_density(count_i) / magnetic_flux_density(injection_grid_number(count_s)) - 1d0)) / sqrt(pi)

                if ( count_i == injection_grid_number(count_s) ) then
                    number_density_diff(count_h, count_s, count_i) = boundary_number_density_diff(count_h, count_s) / sqrt(pi) &
                        & * magnetic_flux_density(count_i) / magnetic_flux_density(injection_grid_number(count_s)) * alpha_perp &
                        & / 2d0 * integral
                
                else
                    number_density_diff(count_h, count_s, count_i) = boundary_number_density_diff(1, count_s) / sqrt(pi) &
                        & * magnetic_flux_density(count_i) / magnetic_flux_density(injection_grid_number(count_s)) * alpha_perp &
                        & / 2d0 * integral
                
                end if

            end do  !count_i

        end do  !count_s

    end do  !count_h
    
end subroutine make_number_density
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine calculation_x_mu_integral(xmin, xlim, xmax, y_variable, coefficient4integral, zero_points_parallel, weights_mu, &
    & weights_parallel, integral_result)
    use constant_parameter
    use constant_in_the_simulation
    !$use omp_lib

    implicit none
    
    double precision, dimension(adiabatic_invariant_grid_number), intent(in) :: xmin, xlim, xmax, y_variable
    double precision, intent(in) :: coefficient4integral
    double precision, dimension(parallel_grid_number), intent(in) :: zero_points_parallel, weights_parallel
    double precision, dimension(adiabatic_invariant_grid_number), intent(in) :: weights_mu
    double precision, intent(out) :: integral_result

    integer :: count_mu, count_x
    double precision, dimension(adiabatic_invariant_grid_number) :: result_halfway, result_halfway_plus, result_halfway_minus

    integral_result = 0d0
    result_halfway = 0d0
    result_halfway_plus = 0d0
    result_halfway_minus = 0d0

    !$omp parallel private(count_x)
    !$omp do
    do count_mu = 1, adiabatic_invariant_grid_number

        if ( xlim(count_mu) > xmin(count_mu) ) then

            do count_x = 1, parallel_grid_number

                result_halfway_plus(count_mu) = result_halfway_plus(count_mu) &
                    & + exp(- ((xlim(count_mu) - xmin(count_mu)) / 2d0 * zero_points_parallel(count_x) &
                    & + (xlim(count_mu) + xmin(count_mu)) / 2d0)**2d0 - y_variable(count_mu) + coefficient4integral) &
                    & * weights_parallel(count_x) * (xlim(count_mu) - xmin(count_mu)) / 2d0

                if ( xmax(count_mu) > xmin(count_mu) ) then

                    result_halfway_minus(count_mu) = result_halfway_minus(count_mu) &
                        & + exp(- ((xmax(count_mu) - xmin(count_mu)) / 2d0 * zero_points_parallel(count_x) &
                        & + (xmax(count_mu) + xmin(count_mu)) / 2d0)**2d0 - y_variable(count_mu) + coefficient4integral) &
                        & * weights_parallel(count_x) * (xmax(count_mu) - xmin(count_mu)) / 2d0

                end if

            end do  !count_x
            
        end if
        
    end do  !count_mu
    !$omp end do
    !$omp end parallel

    result_halfway = result_halfway_plus + result_halfway_minus

    do count_mu = 1, adiabatic_invariant_grid_number
        
        integral_result = integral_result + result_halfway(count_mu) * weights_mu(count_mu)

    end do  !count_mu

end subroutine calculation_x_mu_integral
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine cannot_reach_check(number_density_diff, injection_grid_number)
    use constant_parameter
    use constant_in_the_simulation
    use boundary_and_initial_conditions
    !$use omp_lib

    implicit none
    
    double precision, dimension(3, boundary_series_number, real_grid_number), intent(inout) :: number_density_diff
    integer, dimension(boundary_series_number), intent(in) :: injection_grid_number

    integer :: cannot_reach_point
    integer :: count_h, count_s, count_i

    do count_h = 1, 3

        !$omp parallel private(count_i, cannot_reach_point)
        !$omp do
        do count_s = 1, boundary_series_number
            
            if ( injection_grid_number(count_s) /= real_grid_number ) then
                cannot_reach_point = 0
                do count_i = injection_grid_number(count_s), real_grid_number

                    if( number_density_diff(count_h, count_s, count_i) < 1d-5 .and. cannot_reach_point == 0) then
                        number_density_diff(count_h, count_s, count_i) = 0d0
                        cannot_reach_point = count_i
                    end if

                    if ( cannot_reach_point /= 0 .and. count_i >= cannot_reach_point ) then
                        number_density_diff(count_h, count_s, count_i) = 0d0
                    end if
                    
                end do  !count_i
            end if

            if ( injection_grid_number(count_s) /= 1 ) then
                cannot_reach_point = 0
                do count_i = injection_grid_number(count_s), 1, -1

                    if( number_density_diff(count_h, count_s, count_i) < 1d-5 .and. cannot_reach_point == 0) then
                        number_density_diff(count_h, count_s, count_i) = 0d0
                        cannot_reach_point = count_i
                    end if

                    if ( cannot_reach_point /= 0 .and. count_i <= cannot_reach_point ) then
                        number_density_diff(count_h, count_s, count_i) = 0d0
                    end if
                    
                end do  !count_i
            end if

        end do  !count_s
        !$omp end do
        !$omp end parallel

    end do  !count_h

end subroutine cannot_reach_check
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_charge_density_from_number_density(number_density_diff, charge_number, charge_density_diff, &
    & charge_density_plus_diff, charge_density_minus_diff)
    use constant_parameter
    use constant_in_the_simulation
    use boundary_and_initial_conditions

    implicit none
    
    double precision, dimension(3, boundary_series_number, real_grid_number), intent(in) :: number_density_diff
    double precision, dimension(boundary_series_number), intent(in) :: charge_number
    double precision, dimension(3, real_grid_number), intent(out) :: charge_density_diff, charge_density_plus_diff
    double precision, dimension(3, real_grid_number), intent(out) :: charge_density_minus_diff

    integer :: count_s

    charge_density_diff = 0d0
    charge_density_plus_diff = 0d0
    charge_density_minus_diff = 0d0

    do count_s = 1, boundary_series_number
        
        charge_density_diff = charge_density_diff + charge_number(count_s) * number_density_diff(:, count_s, :)

        if ( charge_number(count_s) > 0d0 ) then
            charge_density_plus_diff = charge_density_plus_diff + charge_number(count_s) * number_density_diff(:, count_s, :)

        else if ( charge_number(count_s) < 0d0 ) then
            charge_density_minus_diff = charge_density_minus_diff + charge_number(count_s) * number_density_diff(:, count_s, :)
            
        end if

    end do  !count_s

end subroutine make_charge_density_from_number_density
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_charge_density_from_Poisson_eq(electrostatic_potential_diff, diff_coordinate_FA, charge_density_poisson_diff)
    use constant_parameter
    use constant_in_the_simulation
    !$use omp_lib

    implicit none
    
    double precision, dimension(3, real_grid_number), intent(in) :: electrostatic_potential_diff
    double precision, dimension(real_grid_number - 1), intent(in) :: diff_coordinate_FA
    double precision, dimension(3, real_grid_number), intent(out) :: charge_density_poisson_diff

    integer :: count_i

    !$omp parallel private(count_i)
    !$omp do
    do count_i = 1, real_grid_number
        
        if ( count_i == 1 .or. count_i == real_grid_number ) then
            charge_density_poisson_diff(:, count_i) = 0d0

        else
            charge_density_poisson_diff(:, count_i) = electrostatic_potential_diff(1, count_i - 1) &
                & / diff_coordinate_FA(count_i - 1) / (diff_coordinate_FA(count_i - 1) + diff_coordinate_FA(count_i))
            
            charge_density_poisson_diff(:, count_i) = charge_density_poisson_diff(:, count_i) &
                & + electrostatic_potential_diff(1, count_i + 1) / diff_coordinate_FA(count_i) &
                & / (diff_coordinate_FA(count_i - 1) + diff_coordinate_FA(count_i))
            
            charge_density_poisson_diff(:, count_i) = charge_density_poisson_diff(:, count_i) &
                & - electrostatic_potential_diff(:, count_i) / diff_coordinate_FA(count_i - 1) / diff_coordinate_FA(count_i)

            charge_density_poisson_diff(:, count_i) = - 2d0 * electric_constant * charge_density_poisson_diff(:, count_i)

        end if

    end do  !count_i
    !$omp end do
    !$omp end parallel

end subroutine make_charge_density_from_Poisson_eq
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_convergence_number(charge_density_diff, charge_density_plus_diff, charge_density_minus_diff, &
    & charge_density_poisson_diff, convergence_number_diff, convergence_number_sum)
    use constant_parameter
    use constant_in_the_simulation

    implicit none
    
    double precision, dimension(3, real_grid_number), intent(in) :: charge_density_diff, charge_density_plus_diff
    double precision, dimension(3, real_grid_number), intent(in) :: charge_density_minus_diff, charge_density_poisson_diff
    double precision, dimension(3, real_grid_number), intent(out) :: convergence_number_diff
    double precision, intent(out) :: convergence_number_sum

    convergence_number_diff = (charge_density_diff - charge_density_poisson_diff)**2d0
    convergence_number_diff = convergence_number_diff / (- charge_density_minus_diff) / charge_density_plus_diff
    
    convergence_number_sum = sqrt(sum(convergence_number_diff(1, :)**2d0) / real_grid_number)
   
end subroutine make_convergence_number
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
!subroutine make_result_min_log(result_file)
!    use constant_parameter
!    use constant_in_the_simulation
!    use boundary_and_initial_conditions
!
!    implicit none
!
!    character(len = 179), intent(out) ::  result_file
!    character(len = 3) :: grid_front, grid_back, file_number, grid_fix
!    character(len = 2) :: alpha_parallel_str, alpha_perp_str
!    character(len = 1) :: boundary_file_number_str
!
!    write(boundary_file_number_str, "(I1.1)") boundary_file_number
!    write(grid_front, "(I3.3)") initial_grid_ionophere_middle_1
!    write(grid_back, "(I3.3)") initial_grid_middle_magnetosphere_1
!    write(grid_fix, "(I3.3)") initial_fix_grid
!    write(file_number, "(I3.3)") boundary_file_number
!    write(alpha_parallel_str, "(I2.2)") int(alpha_parallel)
!    write(alpha_perp_str, "(I2.2)") int(alpha_perp)
!
!    result_file = result_file_front // alpha_perp_str // "_parallel_" // alpha_parallel_str // "/grid_" // grid_front // "_" // &
!        & grid_back // "_" // grid_fix // "/boundary_conditions_" // boundary_file_number_str // "/number_density_iteration/" &
!        & // "min_convergence_log.csv"
!    
!end subroutine make_result_min_log
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine makedirs(input_dir)

    implicit none
    
    character(len=*), intent(in) :: input_dir
    logical :: exist

    character(len=1024) :: command

    inquire(file = trim(input_dir), exist = exist)

    if (exist .eqv. .false.) then
        command = "mkdir -p " // trim(input_dir)
        command = trim(command)

        call system(command)
    end if

end subroutine makedirs
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_result_file_name(initial_min_grid_1, result_file)
    use constant_parameter
    use constant_in_the_simulation
    use boundary_and_initial_conditions

    implicit none
    
    integer, intent(in) :: initial_min_grid_1
    character(len = *), intent(out) ::  result_file
    character(len = 174) :: result_dir
    character(len = 3) :: grid_front, grid_back, file_number, min_grid, grid_fix
    character(len = 2) :: alpha_parallel_str, alpha_perp_str
    character(len = 2) :: boundary_file_number_str

    write(boundary_file_number_str, "(I2.2)") boundary_file_number
    write(grid_front, "(I3.3)") initial_grid_ionophere_middle_1
    write(grid_back, "(I3.3)") initial_grid_middle_magnetosphere_1
    write(grid_fix, "(I3.3)") initial_fix_grid
    write(file_number, "(I3.3)") boundary_file_number
    write(min_grid, "(I3.3)") initial_min_grid_1
    write(alpha_parallel_str, "(I2.2)") int(alpha_parallel)
    write(alpha_perp_str, "(I2.2)") int(alpha_perp)

    result_dir = result_file_front // alpha_perp_str // "_parallel_" // alpha_parallel_str // "/grid_" // grid_front // "_" // &
        & grid_back // "_" // grid_fix // "/boundary_condition_" // boundary_file_number_str // "/number_density_iteration"

    call makedirs(result_dir)

    result_file = result_file_front // alpha_perp_str // "_parallel_" // alpha_parallel_str // "/grid_" // grid_front // "_" // &
        & grid_back // "_" // grid_fix // "/boundary_condition_" // boundary_file_number_str // "/number_density_iteration/min_" &
        & // min_grid // ".csv"
    
end subroutine make_result_file_name
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_result_file_format(format_character)
    use constant_parameter
    use constant_in_the_simulation
    use boundary_and_initial_conditions

    implicit none
    
    character(len = 33), intent(out) ::  format_character
    character(len = 2) :: series_number

    write(series_number, "(I2.1)") boundary_series_number + 9

    format_character = "(1PE25.15E3, " // series_number // "(',', 1PE25.15E3))"

end subroutine make_result_file_format
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine Newton_method_for_electrostatic_potential(electrostatic_potential_diff, convergence_number_diff, initial_min_grid_1, &
    & initial_min_grid_2, electrostatic_potential)
    use constant_parameter
    use constant_in_the_simulation
    use boundary_and_initial_conditions
    !$use omp_lib

    implicit none
    
    double precision, dimension(3, real_grid_number), intent(in) :: electrostatic_potential_diff, convergence_number_diff
    integer, intent(in) :: initial_min_grid_1, initial_min_grid_2
    double precision, dimension(real_grid_number), intent(out) :: electrostatic_potential

    integer :: count_i, min_loc(1)
    double precision :: update
    double precision, dimension(initial_min_grid_1 - 1) :: sorting_potential_1, sorting_potential_1_reverse
    double precision, dimension(real_grid_number - initial_min_grid_2) :: sorting_potential_2
    double precision, dimension(initial_fix_grid - initial_min_grid_1) :: sorting_potential_3
    double precision, dimension(initial_min_grid_2 - initial_fix_grid) :: sorting_potential_4, sorting_potential_4_reverse

    !----------
    ! iteration
    !----------

    electrostatic_potential = electrostatic_potential_diff(1, :)

    !$omp parallel private(update)
    !$omp do
    do count_i = 2, real_grid_number - 1
        
        if ( count_i /= initial_fix_grid .and. (convergence_number_diff(1, count_i) > convergence_number_diff(2, count_i) &
            & .or. convergence_number_diff(1, count_i) > convergence_number_diff(3, count_i)) &
            & .and. convergence_number_diff(2, count_i) /= convergence_number_diff(3, count_i) ) then
            
            if ( convergence_number_diff(1, count_i) == convergence_number_diff(2, count_i) ) then
                update = (electrostatic_potential_diff(1, count_i) - electrostatic_potential_diff(3, count_i)) &
                    & / (convergence_number_diff(1, count_i) - convergence_number_diff(3, count_i)) &
                    & * convergence_number_diff(1, count_i)

            else if ( convergence_number_diff(1, count_i) == convergence_number_diff(3, count_i) ) then
                update = (electrostatic_potential_diff(2, count_i) - electrostatic_potential_diff(1, count_i)) &
                    & / (convergence_number_diff(2, count_i) - convergence_number_diff(1, count_i)) &
                    & * convergence_number_diff(1, count_i)

            else
                update = ((electrostatic_potential_diff(1, count_i) - electrostatic_potential_diff(3, count_i)) &
                    & / (convergence_number_diff(1, count_i) - convergence_number_diff(3, count_i)) &
                    & + (electrostatic_potential_diff(2, count_i) - electrostatic_potential_diff(1, count_i)) &
                    & / (convergence_number_diff(2, count_i) - convergence_number_diff(1, count_i))) / 2d0 &
                    & * convergence_number_diff(1, count_i)
                
            end if

            if ( abs(update) <= 1d0 .and. update == update) then
                electrostatic_potential(count_i) = electrostatic_potential_diff(1, count_i) - update
            
            else if ( abs(update) > 1d0 .and. sqrt(convergence_number_diff(1, count_i)) <= 1d1 .and. update == update ) then
                electrostatic_potential(count_i) = electrostatic_potential_diff(1, count_i) - update / abs(update) &
                    & * sqrt(convergence_number_diff(1, count_i))

            else if ( abs(update) > 1d0 .and. sqrt(convergence_number_diff(1, count_i)) > 1d1 .and. update == update ) then
                electrostatic_potential(count_i) = electrostatic_potential_diff(1, count_i) - update / abs(update) * sqrt(1d1)
                
            end if

        end if

    end do  !count_i
    !$omp end do
    !$omp end parallel


    !--------
    ! sorting
    !--------

    sorting_potential_1 = electrostatic_potential(2:initial_min_grid_1)
    do count_i = 1, initial_min_grid_1 - 1

        if ( sorting_potential_1(count_i) > electrostatic_potential(1) ) then
            sorting_potential_1(count_i) = 2d0 * electrostatic_potential(1) &
                & - sorting_potential_1(count_i) 
        end if

    end do  !count_i    
    call heapsort(initial_min_grid_1 - 1, sorting_potential_1)


    sorting_potential_2 = electrostatic_potential(initial_min_grid_2:real_grid_number-1)
    do count_i = 1, real_grid_number - initial_min_grid_2

        
        if ( sorting_potential_2(count_i) > electrostatic_potential(real_grid_number) ) then
            sorting_potential_2(count_i) = 2d0 * electrostatic_potential(real_grid_number) &
                & - sorting_potential_2(count_i) 
        end if

    end do

    call heapsort(real_grid_number - initial_min_grid_2, sorting_potential_2)

    do count_i = 1, initial_min_grid_1 - 1
        
        sorting_potential_1_reverse(count_i) = sorting_potential_1(initial_min_grid_1 - count_i)

    end do  !count_i

    if ( initial_fix_grid > initial_min_grid_1 .and. initial_fix_grid < initial_min_grid_2 ) then
        sorting_potential_3 = electrostatic_potential(initial_min_grid_1:initial_fix_grid - 1)
        sorting_potential_4 = electrostatic_potential(initial_fix_grid + 1:initial_min_grid_2)

        do count_i = 1, initial_fix_grid - initial_min_grid_1

            if ( sorting_potential_3(count_i) > electrostatic_potential(initial_fix_grid) ) then
                sorting_potential_3(count_i) = 2d0 * electrostatic_potential(initial_fix_grid) &
                    & - sorting_potential_3(count_i)
                !print *, count_i, sorting_potential_3(count_i)
            end if
            
        end do  !count_i

        do count_i = 1, initial_min_grid_2 - initial_fix_grid

            if ( sorting_potential_4(count_i) > electrostatic_potential(initial_fix_grid) ) then
                sorting_potential_4(count_i) = 2d0 * electrostatic_potential(initial_fix_grid) &
                    & - sorting_potential_4(count_i)
            end if
            
        end do  !count_i

        call heapsort(initial_fix_grid - initial_min_grid_1, sorting_potential_3)
        call heapsort(initial_min_grid_2 - initial_fix_grid, sorting_potential_4)
        
        do count_i = 1, initial_min_grid_2 - initial_fix_grid
        
            sorting_potential_4_reverse(count_i) = sorting_potential_4(initial_min_grid_2 - initial_fix_grid - count_i + 1)
    
        end do  !count_i

    end if

    !$omp parallel
    !$omp do
    do count_i = 1, real_grid_number
        
        if ( count_i > 1 .and. count_i < initial_min_grid_1) then
            electrostatic_potential(count_i) = sorting_potential_1_reverse(count_i - 1)
            !print *, count_i, count_i-1, sorting_potential_1_reverse(count_i - 1)
        
        else if ( initial_min_grid_1 <= count_i .and. count_i < initial_fix_grid .and. initial_fix_grid > initial_min_grid_1 &
            & .and. initial_fix_grid < initial_min_grid_2 ) then
            electrostatic_potential(count_i) = sorting_potential_3(count_i - initial_min_grid_1 + 1)
            !print *, count_i, count_i - initial_min_grid_1 + 1, sorting_potential_3(count_i - initial_min_grid_1 + 1)
        
        else if ( initial_fix_grid < count_i .and. count_i <= initial_min_grid_2 .and. initial_fix_grid > initial_min_grid_1 &
            & .and. initial_fix_grid < initial_min_grid_2 ) then
            electrostatic_potential(count_i) = sorting_potential_4_reverse(count_i - initial_fix_grid)
            !print *, count_i, count_i - initial_fix_grid, sorting_potential_4_reverse(count_i - initial_fix_grid),&
            !& electrostatic_potential(real_grid_number-count_i+1)
        
        else if ( count_i > initial_min_grid_2 .and. count_i < real_grid_number ) then
            electrostatic_potential(count_i) = sorting_potential_2(count_i - initial_min_grid_2 + 1)
            !print *, count_i, count_i-initial_min_grid_2, sorting_potential_2(count_i - initial_min_grid_2 + 1),&
            !& electrostatic_potential(real_grid_number-count_i+1)

        end if

    end do  !count_i
    !$omp end do
    !$omp end parallel

    !$omp parallel
    !$omp do
    do count_i = 2, (real_grid_number - 1)/2
        
        if ( electrostatic_potential(count_i) /= electrostatic_potential(count_i) ) then
            electrostatic_potential(count_i) = electrostatic_potential_diff(1, count_i)
            print *, count_i
        end if 

        electrostatic_potential(count_i) = electrostatic_potential(real_grid_number + 1 - count_i)

    end do  !count_i
    !$omp end do
    !$omp end parallel

    do
        min_loc = minloc(electrostatic_potential)
        if ( electrostatic_potential(initial_min_grid_1) == electrostatic_potential(min_loc(1)) ) then
            exit
        else
            electrostatic_potential(min_loc(1)) = electrostatic_potential(initial_min_grid_1)
        end if
    end do

end subroutine Newton_method_for_electrostatic_potential
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine heapsort(n,array)
    ! https://slpr.sakura.ne.jp/qp/sortf90/
    implicit none
    integer,intent(in) :: n
    double precision,intent(inout) :: array(1:n)
   
    integer ::i,k,j,l
    double precision :: t
   
    if(n.le.0)then
       write(6,*)"Error, at heapsort"; stop
    endif
    if(n.eq.1)return
  
    l=n/2+1
    k=n
    do while(k.ne.1)
       if(l.gt.1)then
          l=l-1
          t=array(L)
       else
          t=array(k)
          array(k)=array(1)
          k=k-1
          if(k.eq.1) then
             array(1)=t
             exit
          endif
       endif
       i=l
       j=l+l
       do while(j.le.k)
          if(j.lt.k)then
             if(array(j).lt.array(j+1))j=j+1
          endif
          if (t.lt.array(j))then
             array(i)=array(j)
             i=j
             j=j+j
          else
             j=k+1
          endif
       enddo
       array(i)=t
    enddo
  
    return
end subroutine heapsort
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine Newton_method_for_boundary_number_density(boundary_number_density_diff, convergence_number_diff, injection_grid_number, &
    & boundary_number_density)
    use constant_parameter
    use constant_in_the_simulation
    use boundary_and_initial_conditions

    implicit none
    
    double precision, dimension(3, boundary_series_number), intent(in) :: boundary_number_density_diff
    double precision, dimension(3, real_grid_number), intent(in) :: convergence_number_diff
    integer, dimension(boundary_series_number), intent(in) :: injection_grid_number
    double precision, dimension(boundary_series_number), intent(out) :: boundary_number_density

    integer :: count_s, check_s, check_i
    double precision :: update

    boundary_number_density = boundary_number_density_diff(1, :)

    do count_s = 1, boundary_series_number

        if ( count_s == boundary_ionosphere_1_variable_species .or. count_s == boundary_ionosphere_2_variable_species &
            & .or. count_s == boundary_magnetosphere_variable_species ) then
            check_s = count_s
            check_i = injection_grid_number(check_s)
            
            if ( convergence_number_diff(1, check_i) > convergence_number_diff(2, check_i) &
            & .or. convergence_number_diff(1, check_i) > convergence_number_diff(3, check_i) ) then

                if ( convergence_number_diff(1, check_i) == convergence_number_diff(2, check_i) ) then
                    update = (boundary_number_density_diff(1, check_s) - boundary_number_density_diff(3, check_s)) &
                        & / (convergence_number_diff(1, check_i) - convergence_number_diff(3, check_i)) &
                        & * convergence_number_diff(1, check_i)

                else if ( convergence_number_diff(1, check_i) == convergence_number_diff(3, check_i) ) then
                    update = (boundary_number_density_diff(2, check_s) - boundary_number_density_diff(1, check_s)) &
                        & / (convergence_number_diff(2, check_i) - convergence_number_diff(1, check_i)) &
                        & * convergence_number_diff(1, check_i)
                
                else
                    update = ((boundary_number_density_diff(1, check_s) - boundary_number_density_diff(3, check_s)) &
                        & / (convergence_number_diff(1, check_i) - convergence_number_diff(3, check_i)) &
                        & + (boundary_number_density_diff(2, check_s) - boundary_number_density_diff(1, check_s)) &
                        & / (convergence_number_diff(2, check_i) - convergence_number_diff(1, check_i))) / 2d0 &
                        & * convergence_number_diff(1, check_i)
                    
                end if

                if ( abs(update) < boundary_number_density_diff(1, check_s) * 1d-1 ) then
                    boundary_number_density(check_s) = boundary_number_density_diff(1, check_s) - update
                else if ( update == update ) then
                    boundary_number_density(check_s) = boundary_number_density_diff(1, check_s) &
                        & - update / abs(update) * boundary_number_density_diff(1, check_s) * 1d-1
                end if

                if ( boundary_number_density(check_s) < 0d0 ) then
                    print *, "boundary_number_density < 0d0", boundary_number_density_diff(1, check_s), update
                    stop
                end if
                
            end if

        end if

    end do  !count_i

    
end subroutine Newton_method_for_boundary_number_density
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_convergence_number_file_name(convergence_number_file)
    use constant_parameter
    use constant_in_the_simulation
    use boundary_and_initial_conditions

    implicit none

    character(len = 206), intent(out) ::  convergence_number_file

    character(len = 3) :: grid_front, grid_back, file_number, min_grid, grid_fix
    character(len = 2) :: alpha_parallel_str, alpha_perp_str
    character(len = 2) :: boundary_file_number_str

    write(boundary_file_number_str, "(I2.2)") boundary_file_number
    write(grid_front, "(I3.3)") initial_grid_ionophere_middle_1
    write(grid_back, "(I3.3)") initial_grid_middle_magnetosphere_1
    write(grid_fix, "(I3.3)") initial_fix_grid
    write(alpha_parallel_str, "(I2.2)") int(alpha_parallel)
    write(alpha_perp_str, "(I2.2)") int(alpha_perp)

    convergence_number_file = result_file_front // alpha_perp_str // "_parallel_" // alpha_parallel_str // "/grid_" // grid_front &
        & // "_" // grid_back // "_" // grid_fix // "/boundary_condition_" // boundary_file_number_str // &
        & "/number_density_iteration/convergence_number_sum_list.csv"
    
end subroutine make_convergence_number_file_name