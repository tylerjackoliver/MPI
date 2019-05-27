module serial

    use problem_constants                                           ! Problem set-up
    use utility_functions                                           ! Non-MPI utility routines
    use pgmio                                                       ! PGM I/O
    use, intrinsic :: iso_fortran_env, only: std_err => error_unit  ! Get portable versions of stderr

    implicit none

    PUBLIC

    integer                             :: Mp                       ! Horizontal chunk per processor
    integer                             :: Np                       ! Vertical chunk per processor

    integer                             :: M_actual                 ! Actual image width: used for testing we have width = M
    integer                             :: N_actual                 ! Actual image height: used for testing we have height = N

    !
    ! Spoofed MPI variable initialisations
    !

    integer                             :: rank                     ! Processor rank
    integer                             :: pool_size                ! Size of the worker pool

    integer                             :: cart_comm                ! Spoofed cartesian communicator

    double precision                    :: average                  ! Average value of array

    integer, allocatable                :: nbrs(:)                  ! Spoofed nbrs array
    integer, allocatable                :: dims(:)                  ! Spoofed dims array

    integer                             :: ierr                     ! Spoofed ierr: we have to use MPI to get wtime()

    contains


        subroutine initialise(pool_size, rank, nbrs, dims, cart_comm)

            !
            ! "Initialises" the serial runtime; spoofs the values of MPI-related variables and prints
            ! welcome messages.
            !
            ! Unused, spoofed variables are given a value of -1.
            !

            integer,               intent(out)      :: pool_size
            integer,               intent(out)      :: rank

            integer, dimension(:), intent(inout)    :: nbrs
            integer, dimension(:), intent(inout)    :: dims
    
            integer,               intent(out)      :: cart_comm

            ! So that we can use mpi_wtime(), we still have to initialise
            ! MPI

            call MPI_INIT(ierr)

            ! Spoof the MPI variables
            
            pool_size = -1;
            rank = 0;

            ! Check that the image is the size we expect

            call pgmsize(fname, M_actual, N_actual)

            if (M_actual .ne. M .or. N_actual .ne. N) then

                ! Write error message to stderr

                write(std_err, *) "Incorrect image size specified."

                error stop "Incorrect image size specified."

            end if

            ! Set the image chunks to be just the image itself

            Mp = M
            Np = N

            ! Set cart_comm to -1
            
            cart_comm = -1

            call util_print_welcomemessage(rank)


        end subroutine initialise 

    
        subroutine send_data(to_send, size_to_send, to_recv, comm)

            !
            ! send_data spoofs the initial transfer of data from the master buffer
            ! to all other working processes.
            !
            ! Inputs
            ! ~~~~~~
            !
            ! to_send: Array of data to send. Type: 2-d double precision.
            ! counts:  Number of elements to send. Type: integer
            ! size_to_recv: Number of elements to receive. Type: integer
            ! comm: Communicator on which to perform the data send. Type: integer
            ! 
            ! Inputs & Outputs
            ! ~~~~~~~~~~~~~~~~
            ! 
            ! to_recv: Array of data to receive into. Type: 2-d double precision
            !
                
            double precision, dimension(:,:), intent(in)    :: to_send
            
            integer,                          intent(in)    :: size_to_send
            
            double precision, dimension(:,:), intent(inout) :: to_recv

            integer,                          intent(in)    :: comm

            to_recv = to_send 

        end subroutine send_data


        subroutine check_criterion(new, old, cart_comm, num_iters, average, delta_global)

            !
            ! check_criterion is a wrapping function that allows the determination of the global
            ! delta, local delta, average value of an array, and also prints a status message.
            !
            ! Here, local_delta is equal to global_delta.
            !
            ! Inputs
            ! ~~~~~~
            !
            ! new: Array of values at the current time-step. Type: 2-d double precision.
            ! old: Array of values ar the previous time-step. Type: 2-d double precision.
            ! cart_comm: Communicator used to retrieve data. Type: integer.
            ! num_iters: Number of iterations that have been performed. Type: integer.
            !
            ! Outputs
            ! ~~~~~~~
            ! 
            ! average: Average value of the new array. Type: 2-d double precision array.
            ! delta_global: Maximum difference between old and new across all processors. Type: double precision.
            !

            double precision, dimension(:,:), intent(in)    :: new
            double precision, dimension(:,:), intent(in)    :: old
    
            integer,                          intent(in)    :: cart_comm
            integer,                          intent(in)    :: num_iters
    
            double precision,                 intent(out)   :: average
            double precision,                 intent(out)   :: delta_global
    
            average = sum(new) / (M * N)

            call util_get_local_delta(new, old, delta_global)

            call util_print_average_max(average, delta_global, num_iters, rank)

        end subroutine check_criterion


        subroutine gather_data(to_send, to_recv, comm)

            !
            ! gather_data spoofs the data retrieval from all other processors at the end of
            ! the computation
            !
            ! Inputs
            ! ~~~~~~
            ! to_send: Array to send to the master process. Type: 2-d double precision array.
            ! send_size: Size of the array to send. Type: integer.
            ! recv_size: Size of the array to receive. Type: integer.
            ! comm: Communicator to gather from. Type: integer.
            !
            ! Outputs
            ! ~~~~~~~
            !
            ! to_recv: Array to receive data into. Type: 2-d double precision array.
            !
    
            double precision, dimension(:,:), intent(in)        :: to_send
    
            double precision, dimension(:,:), intent(out)       :: to_recv
    
            integer,                          intent(in)        :: comm
   
            ! Not using comm will throw a warning.

            to_recv = to_send
    
        end subroutine gather_data


        subroutine finalise()

            ! Placeholder subroutine to allow this function to be called in 
            ! all versions of the source code; nothing actually happens!

        end subroutine finalise


        subroutine util_printfinish(num_iters, time_start, time_finish, rank)

            !
            ! util_printfinish prints a message to stdout that gives the total time of the computation.
            !
            ! Inputs
            ! ~~~~~~
            ! num_iters: Number of iterations performed. Type: integer.
            ! time_start: CPU time at the start of the computation. Type: double precision.
            ! time_finish: CPU time at the end of the computation. Type: double precision.
            ! rank: Rank of the processor calling the function. Type: integer.
            !

            integer,            intent(in)  :: num_iters

            double precision,   intent(in)  :: time_start
            double precision,   intent(in)  :: time_finish

            integer,            intent(in)  :: rank

            if (num_iters .eq. max_iters) then

                ! Print a message to stderr

                write(std_err, *) "Maximum number of iterations reached."

                error stop "Maximum number of iterations reached. Stopping..."

            end if

            if (rank .eq. 0) then

                print '(A)',          "  ****************************************************************"
                print '(A, F8.4, A)', "             Finished computing in ", time_finish - time_start, " seconds."
                print '(A)',          "  ****************************************************************"

            end if

        end subroutine util_printfinish


end module serial
