# 1 "problem_constants.f90"
module problem_constants

    implicit none

    PUBLIC

    integer,            parameter :: num_dims           = 2                     ! Number of dimensions to compute - must match compilation flag
    integer,            parameter :: nbrs_len           = 4                     ! Must be twice nbrs_len
    integer,            parameter :: P                  = 24                    ! Number of processors to run on - must match TORQUE/SLURM
    integer,            parameter :: M                  = 768                   ! Width of image, pixels
    integer,            parameter :: N                  = 768                   ! Height of image, pixels
    integer,            parameter :: max_iters          = 10000                 ! Maximum number of iterations
    integer,            parameter :: check_int          = 100                   ! Iteration counts at which to check the value of the iteration array
    
    double precision,   parameter :: stopping_criterion = 0.1                   ! Stopping criterion

    character(*),       parameter :: fname              = "edge768x768.pgm"     ! Filename

end module problem_constants
