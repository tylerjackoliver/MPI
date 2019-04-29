module problem_constants

    implicit none

    PUBLIC

    integer,            parameter :: num_dims           = 2
    integer,            parameter :: nbrs_len           = 4  ! Must be twice nbrs_len
    integer,            parameter :: P                  = 2
    integer,            parameter :: M                  = 768
    integer,            parameter :: N                  = 768
    integer,            parameter :: max_iters          = 10000
    integer,            parameter :: check_int          = 100
    
    double precision,   parameter :: stopping_criterion = 0.1

    character(*),       parameter :: fname              = "edge768x768.pgm"

end module problem_constants
