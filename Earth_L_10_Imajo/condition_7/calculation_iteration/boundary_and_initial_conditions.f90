module boundary_and_initial_conditions
    use constant_in_the_simulation

    implicit none

    !------------------
    ! initial condition
    !------------------

    integer, parameter :: initial_fix_grid = 175
    !(initial_min_grid_1 < initial_fix_grid < initial_min_grid_2 => fix)
    !(initial_fix_grid <= initial_min_grid_1 or initial_min_grid_2 <= initial_fix_grid => not fix)

    integer, parameter :: initial_grid_ionophere_middle_1 = 14
    ! 1 ~ initial_grid_ionophere_middle_1 - 1   initial_electrostatic_potential_ionosphere
    
    integer, parameter :: initial_grid_middle_magnetosphere_1 = 109
    ! initial_grid_ionosphere_middle_1 ~ initial_grid_middle_magnetosphere_1 - 1    initial_electrostatic_potential_middle

    integer, parameter :: initial_grid_middle_magnetosphere_2 = real_grid_number + 2 - initial_grid_middle_magnetosphere_1
    ! initial_grid_middle_magnetosphere_1 ~ initial_grid_middle_magnetosphere_2 - 1    initial_electrostatic_potential_magnetosphere

    integer, parameter :: initial_grid_ionophere_middle_2 = real_grid_number + 2 - initial_grid_ionophere_middle_1
    ! initial_grid_middle_magnetosphere_2 ~ initial_grid_ionophere_middle_2 - 1     initial_electrostatic_potential_middle
    ! initial_grid_ionophere_middle_2 ~ real_grid_number    initial_electrostatic_potential_ionosphere

    integer, parameter :: initial_min_grid_start = initial_grid_middle_magnetosphere_1
    integer, parameter :: initial_min_grid_end = initial_fix_grid

    double precision, parameter :: initial_electrostatic_potential_ionosphere = 1d4     ![V] = [kg m2 s-3 A-1]
    double precision, parameter :: initial_electrostatic_potential_middle = 9d3       ![V] = [kg m2 s-3 A-1]
    double precision, parameter :: initial_electrostatic_potential_magnetosphere = 0d0  ![V] = [kg m2 s-3 A-1]


    !-------------------
    ! boundary condition
    !-------------------

    character(len=25), parameter :: boundary_file_front = '../boundary_condition.csv'
    integer, parameter :: boundary_file_number = 7  !boundary_file
    integer, parameter :: boundary_series_number = 8
    integer, parameter :: boundary_ionosphere_1_variable_species = 3
    integer, parameter :: boundary_ionosphere_2_variable_species = 6
    integer, parameter :: boundary_magnetosphere_variable_species = 7
    !(1 <= boundary_magnetosphere_variable_species <= boundary_series_number  =>  turn on,  else  =>  turn off)


    !------------
    ! result file
    !------------

    character(len = 96), parameter :: result_file_front = &
        & '../../../../../../../mnt/j/plasma_distribution_solver_after_publish/Earth_L_10_Imajo/alpha_perp_'

end module boundary_and_initial_conditions