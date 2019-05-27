module problem_constants

    implicit none

    PUBLIC

    integer,            parameter :: num_dims           = 2                     ! Number of dimensions to compute - must match compilation flag
    integer,            parameter :: nbrs_len           = 4                     ! Must be twice nbrs_len
    integer                       :: P
    integer,            parameter :: M                  = 512                   ! Width of image, pixels
    integer,            parameter :: N                  = 384                   ! Height of image, pixels
    integer,            parameter :: max_iters          = 10000                 ! Maximum number of iterations
    integer,            parameter :: check_int          = 100                   ! Iteration counts at which to check the value of the iteration array
    
    double precision,   parameter :: stopping_criterion = 0.1                   ! Stopping criterion

    character(*),       parameter :: fname              = "edge512x384.pgm"     ! Filename

end module problem_constants
