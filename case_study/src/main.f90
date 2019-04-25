!
! Jack Tyler: Image re/de-composition code; serial version.
!
! Modified 25/4 to compute the average value of the numbers every 10% of iterations.
! Modified 25/4 to cease computation whenever the difference slips below 0.1.
!

program main

    use mpi                                                 ! Use MPI library
    use pgmio                                               ! Use .pgm reading library

    !
    ! Variable declarations for PGM file I/O
    !

    integer                             :: M                ! Number of pixels horizontally
    integer                             :: N                ! Number of pixels vertically

    integer                             :: Mp               ! Horizontal image slice
    integer                             :: Np               ! Vertical image slice

    integer                             :: i                ! Loop variables
    integer                             :: j                ! Loop variables
    integer                             :: k                ! Loop variables
    
    integer                             :: num_iters
    integer                             :: local_sum

    integer, parameter                  :: max_iters = 1000000 ! Number of iterations of jacobi
    integer, parameter                  :: P = 2            ! Number of cores

    integer, parameter                  :: check_int = 100

    !
    ! PGM data arrays
    !

    double precision, allocatable       :: old(:,:)   ! "old" array: for iteration
    double precision, allocatable       :: buf(:,:)         ! Buffer_array: initial read-in location
    double precision, allocatable       :: new_array(:,:)   ! "new" array: for iteration
    double precision, allocatable       :: edge(:,:)        ! Array containing edge values
    double precision, allocatable       :: masterbuf(:,:)   ! Data array for rank 0

    double precision, parameter         :: frac = .25d0     ! Optimisation: compute fraction

    double precision                    :: delta_global = 10! Delta flag: set unnecessarily high

    double precision                    :: delta = 100
    double precision                    :: temp

    !
    ! MPI variable initialisation
    !

    integer                             :: comm
    integer                             :: rank
    integer                             :: pool_size
    integer                             :: ierr

    !
    ! Cartesian topology
    !

    integer                             :: cart_comm        ! Cartesian topology communicator
    integer, dimension(mpi_status_size) :: recv_status      ! Receive status, non-blocking send
    integer                             :: request          ! Request wait
    integer                             :: num_dims         ! Number of dimensions
    integer, allocatable                :: dims(:)          ! Array of size num_dims: alloc later
    integer, allocatable                :: nbrs(:)        ! Array for addresses of neighbours

    integer                             :: left             ! Index of the 'left' direction
    integer                             :: right            ! Index of the 'right' direction
    
    integer, parameter                  :: x_dir = 0        ! What direction is +x?
    integer, parameter                  :: displacement = 1 ! Displacement for cart. top.

    integer                             :: average

    logical                             :: reorder          ! Are we to reorder the dims?
    logical, allocatable                :: periodic(:)      ! Logical array for periodic BCs

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
    allocate(new_array(0:Mp+1, 0:Np+1))
    allocate(edge(0:Mp+1, 0:Np+1))
    allocate(buf(Mp, Np))

    ! Initialise MPI, compute the size and the rank of the array

    call MPI_INIT(ierr)

    comm = MPI_COMM_WORLD

    call MPI_COMM_RANK(comm, rank, ierr)
    call MPI_COMM_SIZE(comm, pool_size, ierr)

    ! If size =/ P, then exit

    if (pool_size .ne. P) then

        call MPI_FINALIZE(ierr)
        error stop "Number of processors does not match the expected number of processes."

    end if

    ! Step 0: Initialise the cartesian topology

    ! Initialise dims and periodic

    allocate(periodic(num_dims))
    allocate(dims(num_dims))
    allocate(nbrs(num_dims*2))

    left = 1
    right = 2
    reorder = .false.

    ! Set dims to zero and periodic to false

    do i = 1, num_dims

        periodic(i) = .false.
        dims(i)     = 0

    end do

    ! Now create the topology

    call MPI_DIMS_CREATE(pool_size, num_dims, dims, ierr)
    call MPI_CART_CREATE(comm, num_dims, dims, periodic, reorder, cart_comm, ierr)
    call MPI_COMM_RANK(cart_comm, rank, ierr)
    call MPI_CART_SHIFT(cart_comm, x_dir, displacement, nbrs(left), nbrs(right), ierr)

    ! Step 1: read the edges data file into the buffer

    if (rank .eq. 0) then

        call pgmread('edge192x128.pgm', masterbuf)

    end if

    ! Now we've read it on the master thread, time to scatter to all procs

    call MPI_SCATTER(masterbuf, Mp * Np, MPI_DOUBLE_PRECISION, buf, Mp * Np, MPI_DOUBLE_PRECISION, 0, comm, ierr)

    ! Step 2: loop over M, N

    if (rank .eq. 0) then

        write(*,*) "Setting up arrays..."

    end if

    do j = 1, Np

        do i = 1, Mp

            edge(i, j) = buf(i, j)                          ! Read into buf
            old(i, j)  = 255                          ! Set old array to white (255)

        end do

    end do

    if (rank .eq. 0) then

        write(*,*) "Iterating..."

    end if

    ! Step 3: Iterate through our edges

    num_iters = 0

    do while (delta_global > 0.1 .and. num_iters < max_iters)

        ! We first need to send the halos from the left and right processes

        ! Send right

        call MPI_ISSEND(old(1, Np), M, MPI_DOUBLE, nbrs(right), 0, cart_comm, request, ierr)
        call MPI_RECV(old(1,0), M, MPI_DOUBLE, nbrs(left), 0, cart_comm, recv_status, ierr)

        call MPI_WAIT(request, recv_status, ierr)

        ! Send left

        call MPI_ISSEND(old(1, 1), M, MPI_DOUBLE, nbrs(left), 0, cart_comm, request, ierr)
        call MPI_RECV(old(1, Np+1), M, MPI_DOUBLE, nbrs(right), 0, cart_comm, recv_status, ierr)

        call MPI_WAIT(request, recv_status, ierr)

        if (mod(num_iters, check_int) .eq. 0) then

             delta = 0

             do j = 1, Np

                 do i = 1, Mp ! Column major

                     new_array(i, j) = frac * ( old(i-1, j) + &
                         old(i+1, j) + old(i, j+1) + &
                         old(i, j-1) - edge(i, j) )

                     ! Compute the local delta: get the difference between the old and new array
                     temp = abs(new_array(i, j) - old(i, j))

                     if (temp > delta) then

                         delta = temp

                     end if

                 end do

             end do

             ! Bias delta againt size of array

             call MPI_ALLREDUCE(delta, delta_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
            
            if (rank .eq. 0) then

                print *, "Total delta:", delta_global, "number of iterations", num_iters
                print *, "Average is:", average

            end if


        else
         
            do j = 1, Np

                 do i = 1, Mp ! Column major

                 new_array(i, j) = frac * ( old(i-1, j) + &
                     old(i+1, j) + old(i, j+1) + &
                     old(i, j-1) - edge(i, j) )


                 end do

            end do

        end if

        old(1:Mp, 1:Np) = new_array(1:Mp, 1:Np)      ! Set old = new, excluding halos

        num_iters = num_iters + 1

    end do

    if (rank .eq. 0) then

        write(*,*) "Done! After", num_iters, "Copying back..."

    end if

    ! Step 4: copy old array back to buf

    buf(1:Mp, 1:Np) = old(1:Mp, 1:Np)

    ! Now gather from buf back to buf

    call MPI_GATHER(buf, Mp * Np, MPI_DOUBLE_PRECISION, masterbuf, Mp * Np, MPI_DOUBLE_PRECISION, 0, comm, ierr)

    ! Write buf to image

    if (rank .eq. 0) then

        call pgmwrite('write_192x128.pgm', masterbuf)

        write(*,*) "Done"

    end if

    deallocate(new_array)
    deallocate(edge)
    deallocate(old)
    deallocate(buf)

    deallocate(periodic)
    deallocate(dims)
    deallocate(nbrs)

    call MPI_FINALIZE(ierr)

contains 


subroutine get_average(array, M, N, comm, average)

    double precision, dimension(:,:), intent(in)    :: array 

    integer,                          intent(in)    :: comm
    integer,                          intent(in)    :: M
    integer,                          intent(in)    :: N

    double precision,                 intent(out)   :: average

    integer                                         :: i, j
    integer                                         :: ierr

    double precision                                :: arr_sum
    
    local_sum = sum(array)

    call MPI_ALLREDUCE(local_sum, average, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)

    average = average / ( M * N )

end subroutine get_average


subroutine get_max(array, comm, max)

    double precision, dimension(:,:), intent(in)    :: array

    integer,                          intent(in)    :: comm 
    
    double precision,                 intent(out)   :: max

    double precision                                :: local_max

    integer                                         :: ierr

    local_max = max(array)

    call MPI_ALLREDUCE(local_max, max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm, ierr)

end subroutine get_max

end program main

