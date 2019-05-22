module utility_functions

    use problem_constants

    implicit none

    contains

    subroutine util_array_copy(new, old)

        double precision, dimension(:,:), intent(inout) :: new
        double precision, dimension(:,:), intent(in)    :: old
    
        integer                                         :: dims_new(2)
        integer                                         :: dims_old(2)
    
        integer                                         :: i, j

        integer                                         :: ierr
    
        dims_new = shape(new)
        dims_old = shape(old)
    
        if ( (dims_new(1) .ne. dims_old(1)) .OR. (dims_new(2) .ne. dims_old(2)) ) then
    
            call MPI_FINALIZE(ierr)
            error stop "Error: tried to util_array_copy but the dimensions didn't match."
    
        end if
    
        ! Use array slicing for now
    
        do j = 1, dims_new(2)
    
            do i = 1,dims_new(1)
    
                new(i, j) = old(i, j)
    
            end do
    
        end do
    
    end subroutine util_array_copy


    subroutine util_array_init(array, init_val)

        double precision, dimension(:, :), intent(inout):: array
        double precision,                  intent(in)   :: init_val

        integer                                         :: dims(2)
        integer                                         :: i, j

        dims = shape(array)

        do j = 1, dims(2)

            do i = 1, dims(1)

                array(i, j) = init_val

            end do

        end do

    end subroutine util_array_init


    subroutine util_print_onrank0(msg, rank)

        character(*), intent(in)    :: msg
        integer,      intent(in)    :: rank

        if (rank .eq. 0) then

            print '(A)', msg

        end if

    end subroutine util_print_onrank0

    
    subroutine util_print_average_max(average, delta, num_iters, rank)

        double precision, intent(in)         :: average

        integer,      intent(in) :: num_iters
        integer,      intent(in) :: rank 

        double precision, intent(in) :: delta

        if (rank .eq. 0) then

            print '(A, I6, A, F6.1, A, F7.3)', "After ", num_iters, &
                                            " iterations the average is ", &
                                            average, " and the delta is ", delta

        end if

    end subroutine util_print_average_max


    subroutine util_get_local_delta(newval, oldval, loc_delta)

        double precision, dimension(:,:), intent(in)    :: newval
        double precision, dimension(:,:), intent(in)    :: oldval 

        double precision,                 intent(out)   :: loc_delta

        integer                                         :: i, j
        integer                                         :: dims(2)

        double precision                                :: temp

        dims = shape(newval)
        loc_delta = 0

        ! Altered loop variables based on zero-indexing

        do j = 2, dims(2)-1

            do i = 2, dims(1)-1

                temp = abs(newval(i, j) - oldval(i, j))

                if (temp > loc_delta) then

                    loc_delta = temp

                end if

            end do
            
        end do

    end subroutine util_get_local_delta

    subroutine util_print_welcomemessage(rank)

        integer, intent(in) :: rank

        if (rank .eq. 0) then

            write(*,*) " ***************************************************************** "
            write(*,*)
            write(*,*) " _____    _              ____       _            _   _"             
            write(*,*) "| ____|__| | __ _  ___  |  _ \  ___| |_ ___  ___| |_(_) ___  _ __  "
            write(*,*) "|  _| / _` |/ _` |/ _ \ | | | |/ _ \ __/ _ \/ __| __| |/ _ \| '_ \ " 
            write(*,*) "| |__| (_| | (_| |  __/ | |_| |  __/ ||  __/ (__| |_| | (_) | | | |"
            write(*,*) "|_____\__,_|\__, |\___| |____/ \___|\__\___|\___|\__|_|\___/|_| |_|"
            write(*,*) "            |___/                                                  "
            write(*,*)
            write(*,*) " ***************************************************************** "

        end if

    end subroutine util_print_welcomemessage


    subroutine util_write_results(num_iters, time_io, time_iterating, time_gather)

        !
        ! Appends the results of this run to file "run_times.dat"
        !

        integer,          intent(in) :: num_iters

        double precision, intent(in) :: time_io
        double precision, intent(in) :: time_iterating
        double precision, intent(in) :: time_gather

        logical :: exist

        inquire(file='run_times.dat', exist=exist)

        if (exist) then

            open(unit=55, file='run_times.dat', status='old', position='append')
            write(*,*) "P, Iterations, IO Time, Iteration time, Time per iteration, Gather time"

        else

            open(unit=55, file='run_times.dat', status='new', action='write')

        end if

        write(55, '(I3, I6, 4F10.7)') P, num_iters, time_io, time_iterating, time_iterating/num_iters, time_gather

        close(55)

    end subroutine util_write_results

end module utility_functions
