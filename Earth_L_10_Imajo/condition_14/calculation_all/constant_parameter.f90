module constant_parameter
    implicit none

    !------------------------------------
    ! mathematical and physical constants
    !------------------------------------

    double precision, parameter :: pi = 4d0*atan(1d0)
    double precision, parameter :: elementary_charge = 1.602176634d-19  ![C] = [s A]
    double precision, parameter :: speed_of_light = 2.99792458d8    ![m s-1]
    double precision, parameter :: magnetic_constant = 1.25663706212d-6     ![kg m s-2 A-2]
    double precision, parameter :: electric_constant = 1d0 / magnetic_constant / speed_of_light**2d0    ![kg-1 m-3 s4 A2]
    double precision, parameter :: constant_of_gravitation = 6.67430d-11    ![kg-1 m3 s-2]
    double precision, parameter :: electron_mass = 9.1093837015d-31     ![kg]

    double precision, parameter :: deg_per_rad = 180d0 / pi
    double precision, parameter :: rad_per_deg = pi / 180d0

end module constant_parameter