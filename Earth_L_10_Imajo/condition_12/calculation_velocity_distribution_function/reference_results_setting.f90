module reference_results_setting

    implicit none
    
    !-----------------------
    ! reference results file
    !-----------------------

    character(len=186), parameter :: &
    & reference_file ='../../../../../../../mnt/j/plasma_distribution_solver_after_publish/Earth_L_10_Imajo/' // &
    &'alpha_perp_12_parallel_12/grid_014_109_175/boundary_condition_12/number_density_iteration/min_149.csv'
    integer, parameter :: boundary_series_number = 12
    integer, parameter :: plot_grid_number = 242


    !------------
    ! result file
    !------------

    character(len=150), parameter :: result_file_front = reference_file(1:150)
    character(len=44), parameter :: result_file_middle_1 = 'velocity_distribution_function/min_' // reference_file(180:182) &
        & // '_grid_'
    character(len=8), parameter :: result_file_middle_2 = '_series_'
    character(len=4), parameter :: result_file_back = '.csv'
    character(len=180) :: result_dir = result_file_front // result_file_middle_1(1:30)

    ! mlat_degree(1), v_perp_i(2), v_para_i(3), distribution_function_i(4), v_perp_b(5), v_para_b(6), distribution_function_b(7)


end module reference_results_setting