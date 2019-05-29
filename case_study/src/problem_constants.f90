module problem_constants

    implicit none

    PUBLIC

    integer,            parameter :: M                  = 256                   ! Width of image, pixels
    integer,            parameter :: N                  = 192                   ! Height of image, pixels
    integer,            parameter :: max_iters          = 10000                 ! Maximum number of iterations
    integer,            parameter :: check_int          = 100                   ! Iteration counts at which to check the value of the iteration array
    
    double precision,   parameter :: stopping_criterion = 0.1                   ! Stopping criterion

    character(*),       parameter :: fname              = "edge256x192.pgm"     ! Filename

end module problem_constants
