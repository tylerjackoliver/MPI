module problem_constants

    implicit none

    PUBLIC

    integer, parameter :: num_dims = 1
    integer, parameter :: nbrs_len = 2  ! Must be twice nbrs_len
    integer, parameter :: P = 2
    integer, parameter :: M = 192
    integer, parameter :: N = 128
    integer, parameter :: max_iters = 10000
    integer, parameter :: check_int = 100

end module problem_constants