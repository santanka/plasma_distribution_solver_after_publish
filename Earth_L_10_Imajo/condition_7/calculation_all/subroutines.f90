subroutine set_omp_num_threads(print_omp)
    use omp_lib

    implicit none
    
    logical, intent(in) :: print_omp
    integer :: max_threads

    ! Get the number of available processors
    max_threads = omp_get_max_threads()!min(omp_get_num_procs(), omp_get_max_threads())
    if ( print_omp .eqv. .true. ) print *, "Maximum threads: ", max_threads

    ! 最大スレッド数で並列処理を実行する
    call omp_set_num_threads(max_threads)

end subroutine set_omp_num_threads
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
    use reference_results_setting
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
subroutine make_potential_energy(mlat_rad, length2planet, length2satellite, charge_number, particle_mass, &
    & electrostatic_potential, potential_energy)
    use constant_parameter
    use constant_in_the_simulation
    use reference_results_setting, only: boundary_series_number
    !$use omp_lib

    implicit none
    
    double precision, dimension(real_grid_number), intent(in) :: mlat_rad, length2planet, length2satellite
    double precision, dimension(boundary_series_number), intent(in) :: charge_number, particle_mass
    double precision, dimension(real_grid_number), intent(in) :: electrostatic_potential
    double precision, dimension(boundary_series_number, real_grid_number), intent(out) ::  potential_energy

    integer :: count_i

    !$omp parallel
    !$omp do
    do count_i = 1, real_grid_number
        
        !gravity of planet
        potential_energy(:, count_i) = - constant_of_gravitation * planet_mass * particle_mass / length2planet(count_i)

        !centrifugal force of planet
        potential_energy(:, count_i) = potential_energy(:, count_i) &
            & - particle_mass * (planet_rotation * length2planet(count_i) * cos(mlat_rad(count_i)))**2d0 / 2d0

        !gravity of satellite
        if ( satellite_mass /= 0d0 ) then
            potential_energy(:, count_i) = potential_energy(:, count_i) &
                & - constant_of_gravitation * satellite_mass * particle_mass / length2satellite(count_i)
        end if
        
        !Coulomb force
        potential_energy(:, count_i) = potential_energy(:, count_i) &
            & + charge_number * electrostatic_potential(count_i)

    end do  !count_i
    !$omp end do
    !$omp end parallel

end subroutine make_potential_energy
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_potential_plus_Bmu(potential_energy, adiabatic_invariant, magnetic_flux_density, potential_plus_Bmu)
    use constant_in_the_simulation
    use reference_results_setting, only: boundary_series_number
    !$use omp_lib

    implicit none
    
    double precision, dimension(boundary_series_number, real_grid_number), intent(in) ::  potential_energy
    double precision, dimension(boundary_series_number, adiabatic_invariant_grid_number), intent(in) :: adiabatic_invariant
    double precision, dimension(real_grid_number), intent(in) :: magnetic_flux_density
    double precision, dimension(boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(out) :: &
        & potential_plus_Bmu

    integer :: count_s, count_i

    do count_s = 1, boundary_series_number

        !$omp parallel
        !$omp do
        do count_i = 1, real_grid_number
            
            potential_plus_Bmu(count_s, count_i, :) = potential_energy(count_s, count_i) &
                & + magnetic_flux_density(count_i) * adiabatic_invariant(count_s, :)

        end do  !count_i
        !$omp end do
        !$omp end parallel
        
    end do  !count_s

end subroutine make_potential_plus_Bmu
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_amin(potential_plus_Bmu, injection_grid_number, amin)
    use constant_in_the_simulation
    use reference_results_setting, only: boundary_series_number
    !$use omp_lib

    implicit none
    
    double precision, dimension(boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(in) :: &
        & potential_plus_Bmu
    integer, dimension(boundary_series_number), intent(in) :: injection_grid_number
    double precision, dimension(boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(out) :: amin

    integer :: count_s, count_i, count_mu, Emax_grid, count4max

    do count_s = 1, boundary_series_number
        
        !$omp parallel private(Emax_grid, count_mu, count4max)
        !$omp do
        do count_i = 1, real_grid_number

            if ( count_i /= injection_grid_number(count_s) ) then
            
                do count_mu = 1, adiabatic_invariant_grid_number

                    if ( count_i < injection_grid_number(count_s) ) then
                        Emax_grid = count_i + 1
                        do count4max = count_i + 1, injection_grid_number(count_s)
                            if ( potential_plus_Bmu(count_s, count4max, count_mu) &
                                & > potential_plus_Bmu(count_s, Emax_grid, count_mu) ) then
                                Emax_grid = count4max
                            end if
                        end do  !count4max
                
                    else if ( count_i > injection_grid_number(count_s) ) then
                        Emax_grid = injection_grid_number(count_s)
                        do count4max = injection_grid_number(count_s), count_i - 1
                            if ( potential_plus_Bmu(count_s, count4max, count_mu) &
                                & > potential_plus_Bmu(count_s, Emax_grid, count_mu) ) then
                            Emax_grid = count4max
                            end if
                        end do  !count4max

                    end if

                    amin(count_s, count_i, count_mu) = sqrt(potential_plus_Bmu(count_s, Emax_grid, count_mu) &
                        & - potential_plus_Bmu(count_s, injection_grid_number(count_s), count_mu))

                    if ( potential_plus_Bmu(count_s, count_i, count_mu) &
                        & > potential_plus_Bmu(count_s, Emax_grid, count_mu) ) then
                        amin(count_s, count_i, count_mu) = 0d0

                    else
                        amin(count_s, count_i, count_mu) &
                            & = sqrt(potential_plus_Bmu(count_s, Emax_grid, count_mu) &
                            & - potential_plus_Bmu(count_s, count_i, count_mu))

                    end if
                
                end do  !count_mu

            else if ( count_i == injection_grid_number(count_s) ) then
                amin(count_s, count_i, :) = 0d0

            end if

        end do  !count_i
        !$omp end do
        !$omp end parallel

    end do  !count_s

end subroutine make_amin
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_alim(potential_plus_Bmu, injection_grid_number, amin, alim, boundary_temperature_para)
    use constant_parameter
    use constant_in_the_simulation
    use reference_results_setting, only: boundary_series_number
    !$use omp_lib

    implicit none
    
    double precision, dimension(boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(in) :: &
        & potential_plus_Bmu
    integer, dimension(boundary_series_number), intent(in) :: injection_grid_number
    double precision, dimension(boundary_series_number), intent(in) :: boundary_temperature_para
    double precision, dimension(boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(in) :: amin
    double precision, dimension(boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(out) :: alim

    integer :: count_s, count_i, count_mu

    do count_s = 1, boundary_series_number

        !$omp parallel private(count_i)
        !$omp do
        do count_mu = 1, adiabatic_invariant_grid_number
            
            alim(count_s, :, count_mu) = boundary_temperature_para(count_s) * alpha_parallel &
                & + potential_plus_Bmu(count_s, injection_grid_number(count_s), count_mu) &
                & - potential_plus_Bmu(count_s, :, count_mu)

            do count_i = 1, real_grid_number

                if ( alim(count_s, count_i, count_mu) <= amin(count_s, count_i, count_mu)**2d0 ) then
                    alim(count_s, count_i, count_mu) = amin(count_s, count_i, count_mu)**2d0
                end if
            
            end do  !count_i

            alim(count_s, :, count_mu) = sqrt(alim(count_s, :, count_mu))

        end do  !count_mu
        !$omp end do
        !$omp end parallel

    end do  !count_s

end subroutine make_alim
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_amax(potential_plus_Bmu, injection_grid_number, boundary_temperature_para, amin, alim, amax)
    use constant_parameter
    use constant_in_the_simulation
    use reference_results_setting, only: boundary_series_number
    !$use omp_lib

    implicit none
    
    double precision, dimension(boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(in) :: &
        & potential_plus_Bmu
    integer, dimension(boundary_series_number), intent(in) :: injection_grid_number
    double precision, dimension(boundary_series_number), intent(in) :: boundary_temperature_para
    double precision, dimension(boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(in) :: amin
    double precision, dimension(boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(in) :: alim
    double precision, dimension(boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(out) :: amax

    integer :: count_s, count_i, count_mu, Emax_grid, count4max
    double precision :: energy_difference, energy_difference_boundary

    do count_s = 1, boundary_series_number
        
        !$omp parallel private(count_mu, count4max, Emax_grid, energy_difference, energy_difference_boundary)
        !$omp do
        do count_i = 1, real_grid_number
            
            if ( (count_i == 1 .and. injection_grid_number(count_s) /= 1) &
                & .or. (count_i == real_grid_number .and. injection_grid_number(count_s) /= real_grid_number) ) then
                amax(count_s, count_i, :) = amin(count_s, count_i, :)

            else if ( count_i == injection_grid_number(count_s) .and. injection_grid_number(count_s) /= 1 &
                & .and. injection_grid_number(count_s) /= real_grid_number ) then
                amax(count_s, count_i, :) = alim(count_s, count_i, :)
            
            else
                do count_mu = 1, adiabatic_invariant_grid_number
                    
                    if ( count_i < injection_grid_number(count_s) .or. injection_grid_number(count_s) == real_grid_number ) then
                        Emax_grid = 1
                        do count4max = 1, count_i - 1

                            if ( potential_plus_Bmu(count_s, count4max, count_mu) &
                                & > potential_plus_Bmu(count_s, Emax_grid, count_mu) ) then
                                Emax_grid = count4max
                            end if

                        end do  !count4max
                    
                    else if (count_i > injection_grid_number(count_s) .or. injection_grid_number(count_s) == 1) then
                        Emax_grid = count_i + 1
                        do count4max = count_i + 1, real_grid_number
                            
                            if ( potential_plus_Bmu(count_s, count4max, count_mu) &
                            & > potential_plus_Bmu(count_s, Emax_grid, count_mu) ) then
                                Emax_grid = count4max
                            end if

                        end do  !count4max

                    end if

                    energy_difference = potential_plus_Bmu(count_s, Emax_grid, count_mu) &
                        & - potential_plus_Bmu(count_s, count_i, count_mu)

                    energy_difference_boundary = potential_plus_Bmu(count_s, injection_grid_number(count_s), count_mu) &
                        & - potential_plus_Bmu(count_s, count_i, count_mu) + boundary_temperature_para(count_s) * alpha_parallel

                    if ( energy_difference <= amin(count_s, count_i, count_mu)**2d0 .or. &
                        & energy_difference_boundary <= amin(count_s, count_i, count_mu)**2d0 ) then
                        amax(count_s, count_i, count_mu) = amin(count_s, count_i, count_mu)

                    else if ( amin(count_s, count_i, count_mu)**2d0 < energy_difference .and. &
                        & energy_difference <= energy_difference_boundary ) then
                        amax(count_s, count_i, count_mu) = sqrt(energy_difference)
                    
                    else if ( amin(count_s, count_i, count_mu)**2d0 < energy_difference .and. &
                        & energy_difference > energy_difference_boundary ) then
                        amax(count_s, count_i, count_mu) = sqrt(energy_difference_boundary)

                    end if

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
subroutine make_number_density(boundary_number_density, boundary_temperature_perp, boundary_temperature_para, &
    & potential_energy, magnetic_flux_density, adiabatic_invariant, injection_grid_number, amin, alim, amax, &
    & zero_points_parallel, weights_mu, weights_parallel, number_density)
    use constant_parameter
    use constant_in_the_simulation
    use reference_results_setting, only: boundary_series_number
    !$use omp_lib

    implicit none
    
    double precision, dimension(boundary_series_number), intent(in) :: boundary_number_density
    double precision, dimension(boundary_series_number), intent(in) :: boundary_temperature_perp, boundary_temperature_para
    double precision, dimension(boundary_series_number, real_grid_number), intent(in) :: potential_energy
    double precision, dimension(real_grid_number), intent(in) :: magnetic_flux_density
    double precision, dimension(boundary_series_number, adiabatic_invariant_grid_number), intent(in) :: adiabatic_invariant
    integer, dimension(boundary_series_number), intent(in) :: injection_grid_number
    double precision, dimension(boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(in) :: &
        & amin, alim, amax
    double precision, dimension(parallel_grid_number), intent(in) :: zero_points_parallel, weights_parallel
    double precision, dimension(adiabatic_invariant_grid_number), intent(in) :: weights_mu
    double precision, dimension(boundary_series_number, real_grid_number), intent(out) :: number_density

    integer :: count_s, count_i
    double precision :: integral, sqrt_temperature_para, coefficient4integral
    double precision, dimension(adiabatic_invariant_grid_number) :: y_variable
    double precision, dimension(adiabatic_invariant_grid_number) :: xmin, xlim, xmax
        
    do count_s = 1, boundary_series_number

        do count_i = 1, real_grid_number
            
            integral = 0d0
            coefficient4integral = - (potential_energy(count_s, count_i) &
                & - potential_energy(count_s, injection_grid_number(count_s))) / boundary_temperature_para(count_s)
            
            y_variable = (magnetic_flux_density(injection_grid_number(count_s)) / boundary_temperature_perp(count_s) &
                & + (magnetic_flux_density(count_i) - magnetic_flux_density(injection_grid_number(count_s))) &
                & / boundary_temperature_para(count_s)) * adiabatic_invariant(count_s, :)

            sqrt_temperature_para = sqrt(boundary_temperature_para(count_s))
            
            xmin = amin(count_s, count_i, :) / sqrt_temperature_para
            xlim = alim(count_s, count_i, :) / sqrt_temperature_para
            xmax = amax(count_s, count_i, :) / sqrt_temperature_para

            call calculation_x_mu_integral(xmin, xlim, xmax, y_variable, coefficient4integral, zero_points_parallel, &
                & weights_mu, weights_parallel, integral)

            
            number_density(count_s, count_i) = boundary_number_density(count_s) / sqrt(pi) &
                & * magnetic_flux_density(count_i) / magnetic_flux_density(injection_grid_number(count_s)) * alpha_perp &
                & / 2d0 * integral

        end do  !count_i

    end do  !count_s
    
end subroutine make_number_density
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine calculation_x_mu_integral(xmin, xlim, xmax, y_variable, coefficient4integral, zero_points_parallel, weights_mu, &
    & weights_parallel, integral_result)
    use constant_parameter
    use constant_in_the_simulation
    use omp_lib

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
subroutine cannot_reach_check(number_density, injection_grid_number)
    use constant_parameter
    use constant_in_the_simulation
    use reference_results_setting, only: boundary_series_number
    !$use omp_lib

    implicit none
    
    double precision, dimension(boundary_series_number, real_grid_number), intent(inout) :: number_density
    integer, dimension(boundary_series_number), intent(in) :: injection_grid_number

    integer :: cannot_reach_point
    integer :: count_s, count_i

    !$omp parallel private(count_i, cannot_reach_point)
    !$omp do
    do count_s = 1, boundary_series_number
            
        if ( injection_grid_number(count_s) /= real_grid_number ) then
            cannot_reach_point = 0
            do count_i = injection_grid_number(count_s), real_grid_number

                if( number_density(count_s, count_i) < 1d-5 .and. cannot_reach_point == 0) then
                    number_density(count_s, count_i) = 0d0
                    cannot_reach_point = count_i
                end if

                if ( cannot_reach_point /= 0 .and. count_i >= cannot_reach_point ) then
                    number_density(count_s, count_i) = 0d0
                end if
                    
            end do  !count_i
        end if

        if ( injection_grid_number(count_s) /= 1 ) then
            cannot_reach_point = 0
            do count_i = injection_grid_number(count_s), 1, -1

                if( number_density(count_s, count_i) < 1d-5 .and. cannot_reach_point == 0) then
                    number_density(count_s, count_i) = 0d0
                    cannot_reach_point = count_i
                end if

                if ( cannot_reach_point /= 0 .and. count_i <= cannot_reach_point ) then
                    number_density(count_s, count_i) = 0d0
                end if
                    
            end do  !count_i
        end if

    end do  !count_s
    !$omp end do
    !$omp end parallel

end subroutine cannot_reach_check
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_particle_flux_density(boundary_number_density, magnetic_flux_density, adiabatic_invariant, injection_grid_number, &
    & boundary_temperature_perp, boundary_temperature_para, particle_mass, potential_energy, alim, amax, &
    & zero_points_parallel, weights_mu, weights_parallel, number_density, particle_flux_density, parallel_mean_velocity)
    use constant_parameter
    use constant_in_the_simulation
    use reference_results_setting, only: boundary_series_number

    implicit none
    
    double precision, dimension(boundary_series_number), intent(in) :: boundary_number_density, particle_mass
    double precision, dimension(boundary_series_number), intent(in) :: boundary_temperature_perp, boundary_temperature_para
    double precision, dimension(real_grid_number), intent(in) :: magnetic_flux_density
    double precision, dimension(boundary_series_number, adiabatic_invariant_grid_number), intent(in) :: adiabatic_invariant
    double precision, dimension(boundary_series_number, real_grid_number), intent(in) :: potential_energy
    double precision, dimension(boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(in) :: &
        & alim, amax
    integer, dimension(boundary_series_number), intent(in) :: injection_grid_number
    double precision, dimension(adiabatic_invariant_grid_number), intent(in) :: weights_mu
    double precision, dimension(parallel_grid_number), intent(in) :: zero_points_parallel, weights_parallel
    double precision, dimension(boundary_series_number, real_grid_number), intent(in) :: number_density
    double precision, dimension(boundary_series_number, real_grid_number), intent(out) :: particle_flux_density
    double precision, dimension(boundary_series_number, real_grid_number), intent(out) :: parallel_mean_velocity

    integer :: count_s, count_i
    double precision :: integral, sqrt_temperature_para, coefficient4integral, flux_direction
    double precision, dimension(adiabatic_invariant_grid_number) :: xlim, xmax, y_variable

    do count_s = 1, boundary_series_number
        
        do count_i = 1, real_grid_number

            if ( number_density(count_s, count_i) /= 0d0 ) then

                if ( count_i < injection_grid_number(count_s) .or. injection_grid_number(count_s) == real_grid_number ) then
                    flux_direction = - 1d0

                else
                    flux_direction = 1d0

                end if

                integral = 0d0
                coefficient4integral = - (potential_energy(count_s, count_i) &
                    & - potential_energy(count_s, injection_grid_number(count_s))) / boundary_temperature_para(count_s)
            
                sqrt_temperature_para = sqrt(boundary_temperature_para(count_s))

                xlim = alim(count_s, count_i, :) / sqrt_temperature_para
                xmax = amax(count_s, count_i, :) / sqrt_temperature_para

                y_variable = (magnetic_flux_density(injection_grid_number(count_s)) / boundary_temperature_perp(count_s) &
                    & + (magnetic_flux_density(count_i) - magnetic_flux_density(injection_grid_number(count_s))) &
                    & / boundary_temperature_para(count_s)) * adiabatic_invariant(count_s, :)

                call calculation_x_exp(xlim, xmax, y_variable, coefficient4integral, zero_points_parallel, weights_mu, &
                    & weights_parallel, integral)
                
                !print *, count_s, count_i, integral, number_density(count_s, count_i)
                if ( integral /= integral ) then
                    stop
                end if

                particle_flux_density(count_s, count_i) = boundary_number_density(count_s) / sqrt(pi) &
                    & * magnetic_flux_density(count_i) / magnetic_flux_density(injection_grid_number(count_s)) &
                    & * sqrt(2d0 * boundary_temperature_para(count_s) / particle_mass(count_s)) * alpha_perp / 2d0 * integral &
                    & * flux_direction
                
                parallel_mean_velocity(count_s, count_i) = particle_flux_density(count_s, count_i) /number_density(count_s, count_i)
                
            else if ( number_density(count_s, count_i) == 0d0 ) then
                particle_flux_density(count_s, count_i) = 0d0
                parallel_mean_velocity(count_s, count_i) = 0d0

            end if

        end do  !count_i

    end do  !count_s

end subroutine make_particle_flux_density
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine calculation_x_exp(xlim, xmax, y_variable, coefficient4integral, zero_points_parallel, weights_mu, &
    & weights_parallel, integral_result)
    use constant_in_the_simulation
    !$use omp_lib

    implicit none
    
    double precision, dimension(adiabatic_invariant_grid_number), intent(in) :: xlim, xmax, y_variable
    double precision, intent(in) :: coefficient4integral
    double precision, dimension(parallel_grid_number), intent(in) :: zero_points_parallel, weights_parallel
    double precision, dimension(adiabatic_invariant_grid_number), intent(in) :: weights_mu
    double precision, intent(out) :: integral_result

    integer :: count_mu, count_x
    double precision, dimension(adiabatic_invariant_grid_number) :: result_halfway

    integral_result = 0d0
    result_halfway = 0d0

    !$omp parallel
    !$omp do
    do count_mu = 1, adiabatic_invariant_grid_number
        
        if ( xlim(count_mu) > xmax(count_mu) ) then

            do count_x = 1, parallel_grid_number
                
                result_halfway(count_mu) = result_halfway(count_mu) + weights_parallel(count_x) &
                    & * ( (xlim(count_mu) - xmax(count_mu)) / 2d0 &
                    & * zero_points_parallel(count_x) + (xlim(count_mu) + xmax(count_mu)) / 2d0 ) &
                    & * exp( - ( (xlim(count_mu) - xmax(count_mu)) / 2d0 * zero_points_parallel(count_x) &
                    & + (xlim(count_mu) + xmax(count_mu)) / 2d0 )**2d0 - y_variable(count_mu) + coefficient4integral)

            end do !count_x

        end if

    end do  !count_y
    !$omp end do
    !$omp end parallel

    integral_result = sum(weights_mu * (xlim - xmax) / 2d0 * result_halfway)
    if ( integral_result /= integral_result ) then
        do count_mu = 1, adiabatic_invariant_grid_number
            print *, count_mu, xlim(count_mu), xmax(count_mu), result_halfway(count_mu)
        end do
    end if

end subroutine calculation_x_exp
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_pressure_perpendicular(boundary_number_density, boundary_temperature_perp, boundary_temperature_para, &
    & potential_energy, injection_grid_number, magnetic_flux_density, adiabatic_invariant, number_density, amin, alim, amax, &
    & zero_points_parallel, weights_mu, weights_parallel, pressure_perp, temperature_perp)
    use constant_parameter
    use constant_in_the_simulation
    use reference_results_setting, only: boundary_series_number
    !$use omp_lib

    implicit none
    
    double precision, dimension(boundary_series_number), intent(in) :: boundary_number_density
    double precision, dimension(boundary_series_number), intent(in) :: boundary_temperature_perp, boundary_temperature_para
    double precision, dimension(boundary_series_number, real_grid_number), intent(in) :: potential_energy
    double precision, dimension(real_grid_number), intent(in) :: magnetic_flux_density
    double precision, dimension(boundary_series_number, adiabatic_invariant_grid_number), intent(in) :: adiabatic_invariant
    integer, dimension(boundary_series_number), intent(in) :: injection_grid_number
    double precision, dimension(boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(in) :: &
        & amin, alim, amax
    double precision, dimension(parallel_grid_number), intent(in) :: zero_points_parallel, weights_parallel
    double precision, dimension(adiabatic_invariant_grid_number), intent(in) :: weights_mu
    double precision, dimension(boundary_series_number, real_grid_number), intent(in) :: number_density
    double precision, dimension(boundary_series_number, real_grid_number), intent(out) :: pressure_perp, temperature_perp

    integer :: count_s, count_i
    double precision :: sqrt_temperature_para, integral, coefficient4integral
    double precision, dimension(adiabatic_invariant_grid_number) :: xmin, xlim, xmax, y_variable

    do count_s = 1, boundary_series_number
        
        do count_i = 1, real_grid_number
            
            if ( number_density(count_s, count_i) /= 0d0 ) then

                integral = 0d0
                coefficient4integral = - (potential_energy(count_s, count_i) &
                    & - potential_energy(count_s, injection_grid_number(count_s))) / boundary_temperature_para(count_s)
            
                y_variable = (magnetic_flux_density(injection_grid_number(count_s)) / boundary_temperature_perp(count_s) &
                    & + (magnetic_flux_density(count_i) - magnetic_flux_density(injection_grid_number(count_s))) &
                    & / boundary_temperature_para(count_s)) * adiabatic_invariant(count_s, :)

                sqrt_temperature_para = sqrt(boundary_temperature_para(count_s))
            
                xmin = amin(count_s, count_i, :) / sqrt_temperature_para
                xlim = alim(count_s, count_i, :) / sqrt_temperature_para
                xmax = amax(count_s, count_i, :) / sqrt_temperature_para

                call calculation_mu_exp_integral(xmin, xlim, xmax, y_variable, coefficient4integral, zero_points_parallel, &
                    & weights_mu, weights_parallel, integral)

                pressure_perp(count_s, count_i) = boundary_number_density(count_s) * boundary_temperature_perp(count_s) / sqrt(pi) &
                    & * boundary_temperature_para(count_s) * magnetic_flux_density(count_i) &
                    & / ((boundary_temperature_para(count_s) - boundary_temperature_perp(count_s)) &
                    & * magnetic_flux_density(injection_grid_number(count_s)) + boundary_temperature_perp(count_s) &
                    & * magnetic_flux_density(count_i)) * magnetic_flux_density(count_i) &
                    & / magnetic_flux_density(injection_grid_number(count_s)) * alpha_perp / 2d0 * integral

                temperature_perp(count_s, count_i) = pressure_perp(count_s, count_i) / number_density(count_s, count_i)
                
            else
                pressure_perp(count_s, count_i) = 0d0
                temperature_perp(count_s, count_i) = 0d0

            end if

        end do  !count_i

    end do  !count_s

end subroutine make_pressure_perpendicular
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine calculation_mu_exp_integral(xmin, xlim, xmax, y_variable, coefficient4integral, zero_points_parallel, weights_mu, &
    & weights_parallel, integral_result)
    use constant_parameter
    use constant_in_the_simulation
    use omp_lib

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
        
        integral_result = integral_result + result_halfway(count_mu) * weights_mu(count_mu) * y_variable(count_mu)

    end do  !count_mu

end subroutine calculation_mu_exp_integral
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_pressure_parallel(boundary_number_density, boundary_temperature_perp, boundary_temperature_para, &
    & injection_grid_number, magnetic_flux_density, adiabatic_invariant, potential_energy, particle_mass, zero_points_parallel, &
    & weights_mu, weights_parallel, number_density, parallel_mean_velocity, amin, alim, amax, pressure_para, temperature_para)
    use constant_parameter
    use constant_in_the_simulation
    use reference_results_setting, only: boundary_series_number

    implicit none
    
    double precision, dimension(boundary_series_number), intent(in) :: boundary_number_density, particle_mass
    double precision, dimension(boundary_series_number), intent(in) :: boundary_temperature_perp, boundary_temperature_para
    integer, dimension(boundary_series_number), intent(in) :: injection_grid_number
    double precision, dimension(real_grid_number), intent(in) :: magnetic_flux_density
    double precision, dimension(boundary_series_number, adiabatic_invariant_grid_number), intent(in) :: adiabatic_invariant
    double precision, dimension(boundary_series_number, real_grid_number), intent(in) :: potential_energy
    double precision, dimension(boundary_series_number, real_grid_number), intent(in) :: number_density, parallel_mean_velocity
    double precision, dimension(boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(in) :: amin, alim
    double precision, dimension(boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(in) :: amax
    double precision, dimension(adiabatic_invariant_grid_number), intent(in) :: weights_mu
    double precision, dimension(adiabatic_invariant_grid_number), intent(in) :: zero_points_parallel, weights_parallel
    double precision, dimension(boundary_series_number, real_grid_number), intent(out) :: pressure_para, temperature_para

    integer :: count_s, count_i
    double precision :: sqrt_temperature_para, integral, coefficient4integral, normal_parallel_mean_velocity
    double precision, dimension(adiabatic_invariant_grid_number) :: y_variable, xlim, xmax, xmin

    do count_s = 1, boundary_series_number
        
        do count_i = 1, real_grid_number
            
            if ( number_density(count_s, count_i) /= 0d0 ) then

                integral = 0d0
                coefficient4integral = - (potential_energy(count_s, count_i) &
                    & - potential_energy(count_s, injection_grid_number(count_s))) / boundary_temperature_para(count_s)
            
                y_variable = (magnetic_flux_density(injection_grid_number(count_s)) / boundary_temperature_perp(count_s) &
                    & + (magnetic_flux_density(count_i) - magnetic_flux_density(injection_grid_number(count_s))) &
                    & / boundary_temperature_para(count_s)) * adiabatic_invariant(count_s, :)
                
                sqrt_temperature_para = sqrt(boundary_temperature_para(count_s))

                xmin = amin(count_s, count_i, :) / sqrt_temperature_para
                xlim = alim(count_s, count_i, :) / sqrt_temperature_para
                xmax = amax(count_s, count_i, :) / sqrt_temperature_para

                normal_parallel_mean_velocity = sqrt(particle_mass(count_s) / 2d0 / boundary_temperature_para(count_s)) &
                    & * abs(parallel_mean_velocity(count_s, count_i))

                call calculation_integral_for_parallel_pressure(xmin, xlim, xmax, y_variable, coefficient4integral, &
                    & normal_parallel_mean_velocity, zero_points_parallel, weights_mu, weights_parallel, integral)

                pressure_para(count_s, count_i) =  boundary_number_density(count_s) / sqrt(pi) * alpha_perp &
                & * boundary_temperature_para(count_s) * magnetic_flux_density(count_i) &
                & / magnetic_flux_density(injection_grid_number(count_s)) * integral

                temperature_para(count_s, count_i) = pressure_para(count_s, count_i) / number_density(count_s, count_i)
                
            else
                pressure_para(count_s, count_i) = 0d0
                temperature_para(count_s, count_i) = 0d0

            end if

        end do  !count_i
        
    end do  !count_s

end subroutine make_pressure_parallel
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine calculation_integral_for_parallel_pressure(xmin, xlim, xmax, y_variable, coefficient4integral, &
    & normal_parallel_mean_velocity, zero_points_parallel, weights_mu, weights_parallel, integral_result)
    use constant_parameter
    use constant_in_the_simulation
    !$use omp_lib

    implicit none
    
    double precision, dimension(adiabatic_invariant_grid_number), intent(in) :: xmin, xlim, xmax, y_variable
    double precision, dimension(adiabatic_invariant_grid_number), intent(in) :: weights_mu
    double precision, dimension(parallel_grid_number), intent(in) :: zero_points_parallel, weights_parallel
    double precision, intent(in) :: normal_parallel_mean_velocity, coefficient4integral
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
                    & + ((xlim(count_mu) - xmin(count_mu)) / 2d0 * zero_points_parallel(count_x) &
                    & + (xlim(count_mu) + xmin(count_mu)) / 2d0 - normal_parallel_mean_velocity)**2d0 &
                    & * exp(- ((xlim(count_mu) - xmin(count_mu)) / 2d0 * zero_points_parallel(count_x) &
                    & + (xlim(count_mu) + xmin(count_mu)) / 2d0)**2d0 - y_variable(count_mu) + coefficient4integral) &
                    & * weights_parallel(count_x) * (xlim(count_mu) - xmin(count_mu)) / 2d0

                if ( xmax(count_mu) > xmin(count_mu) ) then

                    result_halfway_minus(count_mu) = result_halfway_minus(count_mu) &
                        & + ((xmax(count_mu) - xmin(count_mu)) / 2d0 * zero_points_parallel(count_x) &
                        & + (xmax(count_mu) + xmin(count_mu)) / 2d0 - normal_parallel_mean_velocity)**2d0 &
                        & * exp(- ((xmax(count_mu) - xmin(count_mu)) / 2d0 * zero_points_parallel(count_x) &
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

end subroutine calculation_integral_for_parallel_pressure
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_pressure_dynamic(number_density, parallel_mean_velocity, particle_mass, pressure_dynamic)
    use constant_in_the_simulation
    use reference_results_setting, only: boundary_series_number

    implicit none
    
    double precision, dimension(boundary_series_number, real_grid_number), intent(in) :: number_density, parallel_mean_velocity
    double precision, dimension(boundary_series_number), intent(in) :: particle_mass
    double precision, dimension(boundary_series_number, real_grid_number), intent(out) :: pressure_dynamic

    integer :: count_s

    do count_s = 1, boundary_series_number
        
        pressure_dynamic(count_s, :) = 5d-1 * number_density(count_s, :) * particle_mass(count_s) &
            & * parallel_mean_velocity(count_s, :)**2d0

    end do  !count_s

    
end subroutine make_pressure_dynamic
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_Alfven_speed(magnetic_flux_density, particle_mass, number_density, Alfven_speed, Alfven_speed_per_lightspeed)
    use constant_parameter
    use constant_in_the_simulation
    use reference_results_setting, only: boundary_series_number

    implicit none
    
    double precision, dimension(real_grid_number), intent(in) :: magnetic_flux_density
    double precision, dimension(boundary_series_number), intent(in) :: particle_mass
    double precision, dimension(boundary_series_number, real_grid_number), intent(in) :: number_density
    double precision, dimension(real_grid_number), intent(out) :: Alfven_speed, Alfven_speed_per_lightspeed

    double precision, dimension(real_grid_number) :: mass_density
    integer :: count_i

    do count_i = 1, real_grid_number
        
        mass_density(count_i) = sum(particle_mass * number_density(:, count_i))

    end do  !count_i

    Alfven_speed_per_lightspeed = magnetic_flux_density / sqrt(magnetic_flux_density**2d0 + mass_density / electric_constant)
    Alfven_speed = Alfven_speed_per_lightspeed * speed_of_light
    
end subroutine make_Alfven_speed
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_inertial_length(charge_number, particle_mass, number_density, ion_inertial_length, electron_inertial_length)
    use constant_parameter
    use constant_in_the_simulation
    use reference_results_setting, only: boundary_series_number

    implicit none

    double precision, dimension(boundary_series_number), intent(in) :: charge_number, particle_mass
    double precision, dimension(boundary_series_number, real_grid_number), intent(in) :: number_density
    double precision, dimension(real_grid_number), intent(out) :: ion_inertial_length, electron_inertial_length

    integer :: count_s
    double precision, dimension(real_grid_number) :: charge_density_ion, mass_density_ion, number_density_electron

    charge_density_ion = 0d0
    mass_density_ion = 0d0
    number_density_electron = 0d0

    do count_s = 1, boundary_series_number
        
        if ( charge_number(count_s) > 0d0 ) then
            charge_density_ion = charge_density_ion + charge_number(count_s) * number_density(count_s, :)
            mass_density_ion = mass_density_ion + particle_mass(count_s) * number_density(count_s, :)

        else if ( charge_number(count_s) < 0d0 ) then
            number_density_electron = number_density_electron + number_density(count_s, :)

        end if

    end do  !count_s

    ion_inertial_length = sqrt(mass_density_ion / magnetic_constant) / charge_density_ion
    electron_inertial_length = sqrt(electron_mass / magnetic_constant / number_density_electron) / elementary_charge
    
end subroutine make_inertial_length
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_Larmor_radius(magnetic_flux_density, number_density, pressure_perp, charge_number, particle_mass, &
    & ion_Larmor_radius, ion_acoustic_radius, electron_Larmor_radius)
    use constant_parameter
    use constant_in_the_simulation
    use reference_results_setting, only: boundary_series_number

    implicit none

    double precision, dimension(real_grid_number), intent(in) :: magnetic_flux_density
    double precision, dimension(boundary_series_number), intent(in) :: charge_number, particle_mass
    double precision, dimension(boundary_series_number, real_grid_number), intent(in) :: number_density, pressure_perp
    double precision, dimension(real_grid_number), intent(out) :: ion_Larmor_radius, ion_acoustic_radius, electron_Larmor_radius

    integer :: count_s
    double precision, dimension(real_grid_number) :: pressure_perp_ion, pressure_perp_electron, mass_density_ion
    double precision, dimension(real_grid_number) :: number_density_electron, number_density_ion, charge_density_ion

    pressure_perp_ion = 0d0
    pressure_perp_electron = 0d0
    mass_density_ion = 0d0
    number_density_ion = 0d0
    number_density_electron = 0d0
    charge_density_ion = 0d0

    do count_s = 1, boundary_series_number
        
        if ( charge_number(count_s) > 0d0 ) then
            pressure_perp_ion = pressure_perp_ion + pressure_perp(count_s, :)
            mass_density_ion = mass_density_ion + particle_mass(count_s) * number_density(count_s, :)
            number_density_ion = number_density_ion + number_density(count_s, :)
            charge_density_ion = charge_density_ion + charge_number(count_s) * number_density (count_s, :)

        else if ( charge_number(count_s) < 0d0 ) then
            pressure_perp_electron = pressure_perp_electron + pressure_perp(count_s, :)
            number_density_electron = number_density_electron + number_density(count_s, :)

        end if

    end do  !count_s

    ion_Larmor_radius = sqrt(2d0 * pressure_perp_ion * mass_density_ion) / magnetic_flux_density / charge_density_ion
    ion_acoustic_radius = sqrt(2d0 * number_density_ion / number_density_electron * mass_density_ion * pressure_perp_electron) &
        & / magnetic_flux_density / charge_density_ion
    electron_Larmor_radius = sqrt(2d0 * pressure_perp_electron * number_density_electron * electron_mass) / magnetic_flux_density &
        & / number_density_electron / elementary_charge

end subroutine make_Larmor_radius
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_current_density(charge_number, particle_flux_density, current_density)
    use constant_in_the_simulation
    use reference_results_setting, only: boundary_series_number

    implicit none
    
    double precision, dimension(boundary_series_number), intent(in) :: charge_number
    double precision, dimension(boundary_series_number, real_grid_number), intent(in) :: particle_flux_density
    double precision, dimension(real_grid_number), intent(out) :: current_density
    
    integer :: count_s

    current_density = 0d0

    do count_s = 1, boundary_series_number
        
        current_density = current_density + charge_number(count_s) * particle_flux_density(count_s, :)

    end do  !count_s

end subroutine make_current_density
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_result_file_format(format_character)
    use reference_results_setting

    implicit none
    
    character(len = 34), intent(out) ::  format_character
    character(len = 3) :: series_number

    write(series_number, "(I3)") 8 * boundary_series_number + 17

    format_character = "(1PE25.15E3, " // series_number // "(',', 1PE25.15E3))"

end subroutine make_result_file_format
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