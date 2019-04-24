!
! Jack Tyler: Image re/de-composition code; serial version.
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

    integer, parameter                  :: num_iters = 10000 ! Number of iterations of jacobi
    integer, parameter                  :: P = 2            ! Number of cores

    !
    ! PGM data arrays
    !

    double precision, allocatable       :: old(:,:)   ! "old" array: for iteration
    double precision, allocatable       :: buf(:,:)         ! Buffer_array: initial read-in location
    double precision, allocatable       :: new_array(:,:)   ! "new" array: for iteration
    double precision, allocatable       :: edge(:,:)        ! Array containing edge values
    double precision, allocatable       :: masterbuf(:,:)   ! Data array for rank 0

    double precision, parameter         :: frac = .25d0     ! Optimisation: compute fraction

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

    write(*,*) "Setting up arrays..."

    do j = 1, Np

        do i = 1, Mp

            edge(i, j) = buf(i, j)                          ! Read into buf
            old(i, j)  = 255                          ! Set old array to white (255)

        end do

    end do

    write(*,*) "Iterating..."

    ! Step 3: Iterate through our edges

    do k = 1, num_iters

        ! We first need to send the halos from the left and right processes

        ! Send right

        call MPI_ISSEND(old(1, Np), M, MPI_DOUBLE, nbrs(right), 0, cart_comm, request, ierr)
        call MPI_RECV(old(1,0), M, MPI_DOUBLE, nbrs(left), 0, cart_comm, recv_status, ierr)

        call MPI_WAIT(request, recv_status, ierr)

        ! Send left

        call MPI_ISSEND(old(1, 1), M, MPI_DOUBLE, nbrs(left), 0, cart_comm, request, ierr)
        call MPI_RECV(old(1, Np+1), M, MPI_DOUBLE, nbrs(right), 0, cart_comm, recv_status, ierr)

        call MPI_WAIT(request, recv_status, ierr)

        do j = 1, Np

            do i = 1, Mp ! Column major

                new_array(i, j) = frac * ( old(i-1, j) + &
                    old(i+1, j) + old(i, j+1) + &
                    old(i, j-1) - edge(i, j) )

            end do

        end do

        old(1:Mp, 1:Np) = new_array(1:Mp, 1:Np)      ! Set old = new, excluding halos

    end do

    write(*,*) "Done! Copying back..."

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

end program main

