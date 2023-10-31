module constant_in_the_simulation
    use constant_parameter

    implicit none

    !------------
    ! grid number
    !------------

    integer, parameter :: real_grid_number = 349    !for real space
    integer, parameter :: adiabatic_invariant_grid_number = 60    !for adiabatic invariant
    integer, parameter :: parallel_grid_number = 30    !for a(parallel) (adiabatic_invariant_grid_number >= parallel_grid_number)

    
    !------------------------
    ! integral velocity limit
    !------------------------

    double precision, parameter :: alpha_parallel = 12d0
    double precision, parameter :: alpha_perp = 12d0


    !-------------------
    ! planet's constants
    !-------------------

    double precision, parameter :: planet_mass = 5.97243d24  ![kg]
    double precision, parameter :: planet_radius = 6.3781d6     ![m]
    double precision, parameter :: planet_equatorial_radius = 6.3781d6      ![m]
    double precision, parameter :: planet_polar_radius = 6.3568d6   ![m]
    double precision, parameter :: planet_start_altitude = 5.0d5    ![m]
    double precision, parameter :: planet_rotation = 2d0 * pi / 86164.09053083288d0   ![rad s-1]
    double precision, parameter :: planet_l_shell = 10d0
    double precision, parameter :: planet_dipole_moment = 4E0 * pi * (planet_radius/1d3)**3d0 / magnetic_constant &
        & * sqrt((29404.8d0)**2d0 + (1450.9d0)**2d0 + 4652.5d0**2d0)  ![A m2]
    ! IGRF-13, 2020, 7.73d22 [A m2]

    
    !--------------------------------------------------------------------------------------------
    ! satellite's constants (if not set the satellite at the equator, set these parameters to 0.)
    !--------------------------------------------------------------------------------------------

    double precision, parameter :: satellite_mass = 0d0     ![kg]
    double precision, parameter :: satellite_l_shell = 0d0

end module constant_in_the_simulation