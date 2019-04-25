!
! Jack Tyler: Image re/de-composition code; serial version.
!
! Modified 25/4 to compute the average value of the numbers every 10% of iterations.
! Modified 25/4 to cease computation whenever the difference slips below 0.1.
!

program main

    use neighbour_indexes
    use mpi                                                 ! Use MPI library
    use pgmio                                               ! Use .pgm reading library

    implicit none

    !
    ! Variable declarations for PGM file I/O
    !

    integer                             :: M                ! Number of pixels horizontally
    integer                             :: N                ! Number of pixels vertically

    integer                             :: Mp               ! Horizontal image slice
    integer                             :: Np               ! Vertical image slice

    integer                             :: i                ! Loop variables
    integer                             :: j                ! Loop variables
    
    integer                             :: num_iters

    integer, parameter                  :: max_iters = 10000! Number of iterations of jacobi
    integer, parameter                  :: P = 2            ! Number of cores

    integer, parameter                  :: check_int = 100

    !
    ! PGM data arrays
    !

    double precision, allocatable       :: old      (:,:)   ! "old" array: for iteration
    double precision, allocatable       :: buf      (:,:)   ! Buffer_array: initial read-in location
    double precision, allocatable       :: new(:,:)   ! "new" array: for iteration
    double precision, allocatable       :: edge     (:,:)   ! Array containing edge values
    double precision, allocatable       :: masterbuf(:,:)   ! Data array for rank 0

    double precision, parameter         :: frac = .25d0     ! Optimisation: compute fraction

    double precision                    :: delta_global = 10! Delta flag: set unnecessarily high

    double precision                    :: delta = 100

    !
    ! MPI variable initialisation
    !

    integer, parameter                  :: comm = MPI_COMM_WORLD
    integer                             :: rank
    integer                             :: pool_size
    integer                             :: ierr

    !
    ! Cartesian topology
    !

    integer                             :: cart_comm        ! Cartesian topology communicator
    integer, dimension(mpi_status_size) :: recv_status      ! Receive status, non-blocking send
    integer                             :: num_dims         ! Number of dimensions
    integer, allocatable                :: nbrs(:)          ! Array for addresses of neighbours

    integer, parameter                  :: x_dir = 0        ! What direction is +x?
    integer, parameter                  :: displacement = 1 ! Displacement for cart. top.

    double precision                    :: average

    !
    ! Pre-run checks
    !

    ! Initialise MPI, check the size

    call MPI_INIT(ierr)

    call MPI_COMM_SIZE(comm, pool_size, ierr)

    ! If size =/ P, then exit

    if (pool_size .ne. P) then

        call MPI_FINALIZE(ierr)
        error stop "Number of processors does not match the expected number of processes."

    end if


    !
    ! Execute program
    !

    ! Define x and y dimensions

    num_dims = 1

    M = 192
    N = 128

    Mp = M
    Np = ceiling(dble(N/P))                                ! Assumes N/P is perfect

    ! Now allocate data arrays

    allocate(masterbuf(M, N))
    allocate(old(0:Mp+1, 0:Np+1))
    allocate(new(0:Mp+1, 0:Np+1))
    allocate(edge(0:Mp+1, 0:Np+1))
    allocate(buf(Mp, Np))
    allocate(nbrs(2 * num_dims))

    ! Step 0: Initialise the cartesian topology

    call initialise_standard_topology_1d(pool_size, comm, cart_comm, rank, nbrs)

    ! Step 1: read the edges data file into the buffer

    if (rank .eq. 0) then

        call pgmread('edge192x128.pgm', masterbuf)

    end if

    ! Now we've read it on the master thread, time to scatter to all procs

    call mpi_send_data_1d(masterbuf, Mp*Np, buf, Mp * Np, cart_comm)

    ! Step 2: Copy arrays

    call print_onrank0("Initialising arrays...", rank)

    call array_copy(edge(1:Mp, 1:Np), buf(1:Mp, 1:Np))

    call array_init(old, 255.d0)

    ! Step 3: Iterate through our edges

    call print_onrank0("Iterating...", rank)

    num_iters = 0

    do while (delta_global > 0.1 .and. num_iters < max_iters)

        ! We first need to send the halos to/from the left and right processes

        call send_halos_1d(old, Np, M, nbrs, cart_comm)

        delta = 0

        do j = 1, Np

            do i = 1, Mp ! Column major

                new(i, j) = frac * ( old(i-1, j) + old(i+1, j) + &
                old(i, j+1) +  old(i, j-1) - edge(i, j) )

            end do

        end do

        if ((mod(num_iters, check_int) .eq. 0)) then

            call get_average(new, M, N, cart_comm, average)
            call get_local_delta(new, old, delta)
            call get_global_delta(delta, cart_comm, delta_global)
            call print_average_max(average, delta_global, num_iters, rank)

        end if

        call array_copy(old(1:Mp, 1:Np), new(1:Mp, 1:Np))

        num_iters = num_iters + 1

    end do

    if (rank .eq. 0) then

        write(*,*) "Done after", num_iters, "! Copying back..."

    end if

    ! Step 4: copy old array back to buf

    call array_copy(buf(1:Mp, 1:Np), old(1:Mp, 1:Np))

    ! Now gather from buf back to buf

    call gather_data(buf, Mp * Np, masterbuf, Mp * Np, comm)

    ! Write buf to image

    if (rank .eq. 0) then

        call pgmwrite('write_192x128.pgm', masterbuf)

    end if

    deallocate(new)
    deallocate(edge)
    deallocate(old)
    deallocate(buf)

    deallocate(nbrs)

    call MPI_FINALIZE(ierr)

contains 

subroutine get_average(array, M, N, comm, average)

    double precision, dimension(:,:), intent(in)    :: array 

    integer,                          intent(in)    :: comm
    integer,                          intent(in)    :: M
    integer,                          intent(in)    :: N

    double precision,                 intent(out)   :: average

    integer                                         :: ierr

    double precision                                :: local_sum
    
    local_sum = sum(array)

    call MPI_ALLREDUCE(local_sum, average, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)

    average = local_sum / ( M * N )

end subroutine get_average


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


subroutine get_global_delta(delta, comm, global_delta)

    
    double precision,                 intent(in)    :: delta

    integer,                          intent(in)    :: comm 

    double precision,                 intent(out)   :: global_delta

    integer                                         :: ierr

    call MPI_ALLREDUCE(delta, global_delta, 1, MPI_DOUBLE_PRECISION, &
                       MPI_MAX, comm, ierr)

end subroutine get_global_delta


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


subroutine get_array_element_1d(new, old, edge, i, j)
    
    double precision, dimension(:,:), intent(inout) :: new

    double precision, dimension(:,:), intent(in)    :: old
    double precision, dimension(:,:), intent(in)    :: edge

    integer,                          intent(in)    :: i
    integer,                          intent(in)    :: j

    double precision                                :: frac = .25d0
    
    double precision                                :: temp_sum
    
    temp_sum = old(i-1, j) + old(i+1, j) + old(i, j+1) +  old(i, j-1)

    new(i, j) = frac * ( temp_sum - edge(i, j) )

end subroutine get_array_element_1d


subroutine array_copy(new, old)

    double precision, dimension(:,:), intent(inout) :: new
    double precision, dimension(:,:), intent(in)    :: old

    integer                                         :: dims_new(2)
    integer                                         :: dims_old(2)

    integer                                         :: i, j

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


subroutine initialise_standard_topology_1d(pool_size, old_comm, new_comm, new_rank, nbrs)
    
    use neighbour_indexes                   ! Provides standardised values of left, right

    integer, intent(in)                     :: pool_size
    integer, intent(in)                     :: old_comm

    integer, intent(out)                    :: new_comm
    integer, intent(out)                    :: new_rank
    integer, intent(inout), dimension(:)    :: nbrs

    integer, parameter                      :: num_dims = 1
    integer                                 :: dims(1)
    integer                                 :: ierr
    integer                                 :: x_dir 
    integer                                 :: displacement

    logical                                 :: periodic(1)
    logical                                 :: reorder(1)

    reorder = .false.

    ! Set dims to zero and periodic to false

    do i = 1, num_dims

        periodic(i) = .false.
        dims(i)     = 0

    end do

    x_dir = 0
    displacement = 1
    
    ! Now create the topology

    call MPI_DIMS_CREATE(pool_size, num_dims, dims, ierr)
    call MPI_CART_CREATE(comm, num_dims, dims, periodic, reorder, cart_comm, ierr)
    call MPI_COMM_RANK(cart_comm, rank, ierr)
    call MPI_CART_SHIFT(cart_comm, x_dir, displacement, nbrs(left), nbrs(right), ierr)

end subroutine initialise_standard_topology_1d


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


subroutine mpi_send_data_1d(to_send, size_to_send, to_recv, size_to_recv, comm)

    double precision, dimension(:,:), intent(in)    :: to_send
    
    integer,                          intent(in)    :: size_to_send
    
    double precision, dimension(:,:), intent(inout) :: to_recv

    integer,                          intent(in)    :: size_to_recv
    integer,                          intent(in)    :: comm

    call MPI_SCATTER(to_send, size_to_send, MPI_DOUBLE_PRECISION, to_recv, &
                     size_to_recv, MPI_DOUBLE_PRECISION, 0, comm, ierr)

end subroutine mpi_send_data_1d


subroutine print_onrank0(msg, rank)

    character(*), intent(in)    :: msg
    integer,      intent(in)    :: rank

    if (rank .eq. 0) then

        print *, msg

    end if

end subroutine print_onrank0


subroutine send_halos_1d(old, Np, M, nbrs, cart_comm)

    use neighbour_indexes

    double precision, dimension(0:,0:), intent(inout)   :: old

    integer,                          intent(in)        :: Np
    integer,                          intent(in)        :: M
    integer,          dimension(:),   intent(in)        :: nbrs
    integer,                          intent(in)        :: cart_comm

    integer, dimension(mpi_status_size)                 :: recv_status      ! Receive status, non-blocking send

    integer                                             :: ierr

    ! Send/recv right

    call MPI_SENDRECV(old(1, Np), M, MPI_DOUBLE, nbrs(right), 0, &
    old(1,0), M, MPI_DOUBLE, nbrs(left), 0, cart_comm, &
    recv_status, ierr)

    ! Send/recv left

    call MPI_SENDRECV(old(1, 1), M, MPI_DOUBLE, nbrs(left), 0, &
    old(1, Np+1), M, MPI_DOUBLE, nbrs(right), 0, cart_comm, &
    recv_status, ierr)

end subroutine send_halos_1d



subroutine gather_data(to_send, send_size, to_recv, recv_size, comm)

    double precision, dimension(:,:), intent(in)        :: to_send

    integer,                          intent(in)        :: send_size
    
    double precision, dimension(:,:), intent(out)       :: to_recv

    integer,                          intent(in)        :: recv_size
    integer,                          intent(in)        :: comm

    integer, dimension(mpi_status_size)                 :: recv_status      ! Receive status, non-blocking send

    integer                                             :: ierr

    call MPI_GATHER(to_send, send_size, MPI_DOUBLE_PRECISION, to_recv, recv_size,&
                    MPI_DOUBLE_PRECISION, 0, comm, ierr)

end subroutine gather_data



end program

