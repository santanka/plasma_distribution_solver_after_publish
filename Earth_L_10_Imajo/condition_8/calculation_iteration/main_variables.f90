module main_variables
    use constant_in_the_simulation
    use boundary_and_initial_conditions

    implicit none
    
    !------------------
    ! calculation field
    !------------------
    
    double precision, dimension(real_grid_number) :: mlat, length2planet, length2satellite, coordinate_FA
    double precision, dimension(real_grid_number - 1) :: diff_coordinate_FA
    double precision, dimension(real_grid_number) :: magnetic_flux_density
    double precision :: planet_mlat_1, planet_mlat_2


    !------------------
    ! initial condition
    !------------------

    double precision, dimension(real_grid_number) :: initial_electrostatic_potential


    !-------------------
    ! boundary condition
    !-------------------

    double precision, dimension(boundary_series_number) :: boundary_number_density, initial_boundary_number_density
    double precision, dimension(boundary_series_number) :: boundary_temperature_perp, boundary_temperature_para
    double precision, dimension(boundary_series_number) :: charge_number
    double precision, dimension(boundary_series_number) :: particle_mass
    integer, dimension(boundary_series_number) :: injection_grid_number


    !-----------------------------
    ! Gaussian-Legendre quadrature
    !-----------------------------

    double precision, dimension(adiabatic_invariant_grid_number) :: zero_points_mu, weights_mu
    double precision, dimension(parallel_grid_number) :: zero_points_parallel, weights_parallel

    !-----------------------------
    ! variables (count_h: 0, +, -)
    !-----------------------------

    double precision, dimension(boundary_series_number, adiabatic_invariant_grid_number) :: adiabatic_invariant
    double precision, dimension(3, real_grid_number) :: electrostatic_potential_diff
    double precision, dimension(3, boundary_series_number) :: boundary_number_density_diff
    double precision, dimension(3, boundary_series_number, real_grid_number) :: potential_energy_diff
    double precision, dimension(3, boundary_series_number, real_grid_number, adiabatic_invariant_grid_number) :: amin, amax, alim
    double precision, dimension(3, boundary_series_number, real_grid_number, adiabatic_invariant_grid_number) :: &
        & potential_plus_Bmu_diff
    double precision, dimension(3, boundary_series_number, real_grid_number) :: number_density_diff
    double precision, dimension(3, real_grid_number) :: charge_density_diff, charge_density_plus_diff, charge_density_minus_diff
    double precision, dimension(3, real_grid_number) :: charge_density_poisson_diff


    !--------------
    ! for iteration
    !--------------

    double precision, dimension(real_grid_number) :: electrostatic_potential
    double precision, dimension(3, real_grid_number) :: convergence_number_diff
    double precision :: convergence_number_sum
    double precision :: convergence_number_sum_min
    integer :: initial_min_grid_1, initial_min_grid_2


    !--------
    ! counter
    !--------

    integer :: count_iteration, count_iteration_timer, count_h, count_s, count_i, count_mu, count_min_position


    !----------
    ! file name
    !----------

    character(len = 188) :: result_file, format_character
    character(len = 205) :: convergence_number_file
    character(len = 25) :: boundary_file
    logical :: file_exist


    !------
    ! dummy
    !------

    character(len = 128) :: dummy
    logical :: print_omp
    integer :: omp_num

    !-----
    ! time
    !-----

    integer :: date_time(8), date_time_start(8), date_time_end(8)
    character(len = 10) :: dummy_1, dummy_2, dummy_3

end module main_variables