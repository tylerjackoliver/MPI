module utility_functions

    implicit none

    contains

    subroutine array_copy(new, old)

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
            error stop "Error: tried to array_copy but the dimensions didn't match."
    
        end if
    
        ! Use array slicing for now
    
        do j = 1, dims_new(2)
    
            do i = 1,dims_new(1)
    
                new(i, j) = old(i, j)
    
            end do
    
        end do
    
    end subroutine array_copy


    subroutine array_init(array, init_val)

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

    end subroutine array_init


    subroutine print_onrank0(msg, rank)

        character(*), intent(in)    :: msg
        integer,      intent(in)    :: rank

        if (rank .eq. 0) then

            print *, msg

        end if

    end subroutine print_onrank0

    
    subroutine print_average_max(average, delta, num_iters, rank)

        double precision, intent(in)         :: average

        integer,      intent(in) :: num_iters
        integer,      intent(in) :: rank 

        double precision, intent(in) :: delta

        if (rank .eq. 0) then

            print '(A, I6, A, F5.1, A, F6.3)', "After ", num_iters, &
                                            " iterations the average is ", &
                                            average, " and the delta is ", delta

        end if

    end subroutine print_average_max


    subroutine get_local_delta(newval, oldval, loc_delta)

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

    end subroutine get_local_delta


end module utility_functions
