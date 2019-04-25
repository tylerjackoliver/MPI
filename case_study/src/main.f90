!
! Jack Tyler: Image re/de-composition code; serial version.
!
! Modified 25/4 to compute the average value of the numbers every 10% of iterations.
! Modified 25/4 to cease computation whenever the difference slips below 0.1.
!

program main

    use neighbour_indexes
    use utility_functions                                   ! Contains non-mpi subroutiens used below
    use mpi                                                 ! Use MPI library
    use mpi_wrappers                                        ! Contains MPI wrappers overloaded to number of dimensions
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
    double precision, allocatable       :: new      (:,:)   ! "new" array: for iteration
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

    call mpi_initialise_standard_topology(num_dims, pool_size, comm, cart_comm, rank, nbrs)

    ! Step 1: read the edges data file into the buffer

    if (rank .eq. 0) then

        call pgmread('edge192x128.pgm', masterbuf)

    end if

    ! Now we've read it on the master thread, time to scatter to all procs

    call mpi_send_data(num_dims, masterbuf, Mp*Np, buf, Mp * Np, cart_comm)

    ! Step 2: Copy arrays

    call print_onrank0("Initialising arrays...", rank)
    call array_copy(edge(1:Mp, 1:Np), buf(1:Mp, 1:Np))
    call array_init(old, 255.d0)

    ! Step 3: Iterate through our edges

    call print_onrank0("Iterating...", rank)

    num_iters = 0

    do while (delta_global > 0.1 .and. num_iters < max_iters)

        ! We first need to send the halos to/from the left and right processes

        call mpi_send_halos(num_dims, old, Np, M, nbrs, cart_comm)

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

end program main

