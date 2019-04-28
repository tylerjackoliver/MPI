!
! Jack Tyler: Image re/de-composition code; serial version.
!
! Modified 25/4 to compute the average value of the numbers every 10% of iterations.
! Modified 25/4 to cease computation whenever the difference slips below 0.1.
!

program main

#ifdef PARALLEL1d

    use wrappers_1d

#elif PARALLEL2d
    
    !     use wrappers_2d
        use mpi 
        use neighbour_indexes
        use problem_constants
        use pgmio
        use utility_functions

#elif SERIAL

    use serial

#endif

    implicit none

    !
    ! Variable declarations for PGM file I/O
    !

    integer                             :: i                ! Loop variables
    integer                             :: j                ! Loop variables
    
    integer                             :: num_iters

    !
    ! PGM data arrays
    !

    double precision, allocatable       :: old(:,:)         ! "old" array: for iteration
    double precision, allocatable       :: buf(:,:)         ! Buffer_array: initial read-in location
    double precision, allocatable       :: new(:,:)         ! "new" array: for iteration
    double precision, allocatable       :: edge(:,:)        ! Array containing edge values
    double precision, allocatable       :: masterbuf(:,:)   ! Data array for rank 0

    double precision, parameter         :: frac = .25d0     ! Optimisation: compute fraction

    double precision                    :: delta_global = 10! Delta flag: set unnecessarily high

    double precision                    :: local_delta = 100

    !
    ! MPI variable initialisation
    !

    integer                             :: rank
    integer                             :: pool_size
    integer                             :: ierr
    integer                             :: cart_comm

    !
    ! Cartesian topology
    !

    integer, allocatable                :: dims(:)          ! Array of size num_dims: alloc later
    integer, allocatable                :: nbrs(:)        ! Array for addresses of neighbours

    double precision                    :: average

    !
    ! Two-dimensional variables
    !

    integer, allocatable                :: counts(:)
    integer, allocatable                :: displacements(:)

    integer                             :: master_type
    integer                             :: subarray_type
    integer                             :: v_halo_type
    integer                             :: h_halo_type

    integer                             :: Mp, Np

    !
    ! Execute program
    !

    ! Define x and y dimensions

    ! Now allocate data arrays

    allocate(nbrs(num_dims*2))
    allocate(dims(num_dims))

    ! Initialise MPI, compute the size and the rank of the array

    write(*,*) "Here, before mpi_initialise"

    call mpi_initialise(num_dims, pool_size, rank, nbrs, dims, cart_comm)

    allocate(masterbuf(M, N))
    
    allocate(old(0:Mp+1, 0:Np+1))
    allocate(new(0:Mp+1, 0:Np+1))
    allocate(edge(0:Mp+1, 0:Np+1))
    allocate(buf(Mp, Np))

    ! Step 1: read the edges data file into the buffer

    if (rank .eq. 0) then

        call pgmread('edge192x128.pgm', masterbuf)

    end if

    ! Now we've read it on the master thread, time to scatter to all procs
    !
    ! The master buffer at least changes, which is good

    call mpi_send_data(num_dims, masterbuf, Mp*Np, buf, Mp * Np, cart_comm)

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

        call mpi_send_halos(num_dims, old, Np, M, dims, nbrs, cart_comm)

        do j = 1, Np

            do i = 1, Mp ! Column major

                new(i, j) = frac * ( old(i-1, j) + &
                    old(i+1, j) + old(i, j+1) + &
                    old(i, j-1) - edge(i, j) )

            end do

        end do

        if ((mod(num_iters, check_int) .eq. 0)) then

            call mpi_get_average(new, M, N, cart_comm, average)

            call util_get_local_delta(new, old, local_delta)
            call mpi_get_global_delta(num_dims, local_delta, cart_comm, delta_global)

            call util_print_average_max(average, delta_global, num_iters, rank)

        end if

        call util_array_copy(old(1:Mp, 1:Np), new(1:Mp, 1:Np))

        num_iters = num_iters + 1

    end do

    if (rank .eq. 0) then

        write(*,*) "Done! After", num_iters, "Copying back..."

    end if

    ! Step 4: copy old array back to buf

    call util_array_copy(buf(1:Mp, 1:Np), old(1:Mp, 1:Np))

    ! Now gather from buf back to masterbuf

    call mpi_gather_data(num_dims, buf, Mp * Np, masterbuf, Mp * Np, cart_comm)

    ! Write buf to image

    if (rank .eq. 0) then

        call pgmwrite('write_192x128.pgm', masterbuf)

        write(*,*) "Done"

    end if

    deallocate(new)
    deallocate(edge)
    deallocate(old)
    deallocate(buf)

    deallocate(dims)
    deallocate(nbrs)

    call mpi_finalize(ierr)

    contains

    subroutine mpi_define_vectors(num_dims, dims, pool_size, Mp, Np)

        integer,                    intent(in)  :: num_dims
        integer,                    intent(in)  :: pool_size
        integer,                    intent(in)  :: Mp
        integer,                    intent(in)  :: Np

        integer, dimension(:),      intent(in)  :: dims

        integer, allocatable                    :: sizes(:)
        integer, allocatable                    :: subsizes(:)
        integer, allocatable                    :: starts(:)

        integer                                 :: i
        integer                                 :: ierr
        integer                                 :: temp_type

        integer(kind=mpi_address_kind)          :: start
        integer(kind=mpi_address_kind)          :: integer_extent
        integer(kind=mpi_address_kind)          :: lb           ! Lower bound
        integer(kind=mpi_address_kind)          :: double_extent
        integer(kind=mpi_address_kind)          :: tot_extent

        ! Allocate arrays based on num_dims

        allocate(sizes(num_dims))
        allocate(subsizes(num_dims))
        allocate(starts(num_dims))

        ! Create the block type, which will allocate the subarrays for each processor
        ! from the master array (i.e. an Mp x Np subarray of an (Mp+2) x (Np+2) array)

        ! Original size of masterbuf

        sizes(1) = Mp + 2
        sizes(2) = Np + 2

        ! New size of the subarray in masterbuf

        subsizes(1) = Mp
        subsizes(2) = Np

        ! Starting location of the subarray
        
        starts(1) = 1
        starts(2) = 1

        call MPI_TYPE_CREATE_SUBARRAY(num_dims, sizes, subsizes, starts, &
                                      MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, subarray_type, ierr)
    
        ! Create the master block, whih provides an array to convert directly from the master buffer
        ! to the working thread

        ! Original size of the masterbuf

        sizes(1) = Mp * dims(1)
        sizes(2) = Np * dims(2)

        ! New size of the subarray in the masterbuf
        
        subsizes(1) = Mp 
        subsizes(2) = Np

        ! Starting location of the subarray

        starts(1) = 0
        starts(2) = 0

        call MPI_TYPE_CREATE_SUBARRAY(num_dims, sizes, subsizes, starts, &
                                      MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, temp_type, ierr)

        ! Vertical halo vectors
        ! Np block, one element per block, Mp+2 elements between blocks

        call MPI_TYPE_VECTOR(Np, 1 , Mp+2, MPI_DOUBLE_PRECISION, v_halo_type, ierr)

        ! Horizontal halo vectors
        ! 1 block, Mp elements in each block, Mp elements between each block

        call MPI_TYPE_VECTOR(1, Mp, Mp, MPI_DOUBLE_PRECISION, h_halo_type, ierr)

        ! So that we can use SCATTERV in the two-dimensional topology, we need to
        ! resize the original master block to fit properly into each quadrant of the
        ! domain. We can do this using MPI_CREATE_RESIZED after getting the length
        ! from MPI_TYPE_EXTENT
        !

        call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION, lb, double_extent, ierr)
        
        start = 0
        tot_extent = Mp * double_extent
        
        call MPI_TYPE_CREATE_RESIZED(temp_type, start, tot_extent, &
                                     master_type, ierr)

        ! Commit these new datatypes so that we can use them

        call MPI_TYPE_COMMIT(master_type, ierr)
        call MPI_TYPE_COMMIT(subarray_type, ierr)
        call MPI_TYPE_COMMIT(v_halo_type, ierr)
        call MPI_TYPE_COMMIT(h_halo_type, ierr)

        deallocate(sizes)
        deallocate(subsizes)
        deallocate(starts)

    end subroutine mpi_define_vectors

    subroutine compute_counts_displacements(num_dims, dims, pool_size, Np, counts, displacements)

        integer,               intent(in)       :: num_dims

        integer, dimension(:), intent(in)       :: dims
        
        integer,               intent(in)       :: pool_size
        integer,               intent(in)       :: Np

        integer, dimension(:), intent(inout)    :: counts
        integer, dimension(:), intent(inout)    :: displacements

        integer                                 :: baseline
        integer                                 :: i

        baseline = 1

        do i = 1, pool_size

                counts(i) = 1
                displacements(i) = (baseline-1) + mod(i-1, dims(1))

                if (mod(i, dims(1)) .eq. 0) then
                
                    baseline = baseline + Np * dims(1)

                end if

        end do

    end subroutine compute_counts_displacements


    subroutine mpi_initialise(num_dims, pool_size, rank, nbrs, dims, cart_comm)

        integer, intent(in)                     :: num_dims

        integer, intent(out)                    :: pool_size
        integer, intent(out)                    :: rank

        integer, dimension(:), intent(inout)    :: nbrs
        integer, dimension(:), intent(inout)    :: dims

        integer,               intent(out)      :: cart_comm

        integer                                 :: i
        integer                                 :: ierr
        integer                                 :: comm

        ! integer, allocatable                    :: counts(:)

        call MPI_INIT(ierr)

        comm = MPI_COMM_WORLD

        call MPI_COMM_SIZE(comm, pool_size, ierr)

        call mpi_initialise_standard_topology(num_dims, dims, cart_comm, nbrs, rank)

        allocate(counts(pool_size))
        allocate(displacements(pool_size))

        call compute_counts_displacements(num_dims, dims, pool_size, Np, counts, displacements)

        call mpi_define_vectors(num_dims, dims, pool_size, Mp, Np)

    end subroutine mpi_initialise


    subroutine mpi_initialise_standard_topology(num_dims, dims, cart_comm, nbrs, rank)

        ! Creates a 2D cartesian communicator with the appropriate dimensions,
        ! compute the neighbours for each process and computes the MP and NP
        ! At the end calls the routine to create the derived datatypes
        
        integer,               intent(in)       :: num_dims
        integer, dimension(:), intent(inout)    :: dims
        integer,               intent(out)      :: cart_comm
        integer, dimension(:), intent(inout)    :: nbrs

        integer,               intent(out)      :: rank

        logical                                 :: periodic(2)
        logical                                 :: reorder

        integer                                 :: ierr
        integer                                 :: x_dir
        integer                                 :: y_dir
        integer                                 :: displacement
        integer                                 :: i

        integer                                 :: pool_size
        
        call MPI_COMM_SIZE(MPI_COMM_WORLD, pool_size, ierr)

        do i = 1, pool_size

            dims(i) = 0
            periodic(i) = .false.

        end do

        reorder      = .false.
        x_dir        = 0
        y_dir        = 1
        displacement = 1

        ! Create cartesian topology 2D

        call MPI_DIMS_CREATE(pool_size, num_dims, dims, ierr)

        call MPI_CART_CREATE(MPI_COMM_WORLD, num_dims, dims, periodic, reorder, cart_comm, ierr)
        
        call MPI_COMM_RANK(cart_comm, rank, ierr)
        
        ! Get neighbours

        call MPI_CART_SHIFT(cart_comm, x_dir, displacement, nbrs(left), nbrs(right), ierr)
        call MPI_CART_SHIFT(cart_comm, y_dir, displacement, nbrs(down), nbrs(up), ierr)
    
        ! Compute the new array dimensions

        if ( (mod(M, dims(1)) .ne. 0) .or. (mod(N, dims(2)) .ne. 0) ) then

            call MPI_FINALIZE(ierr)
            error stop "Error: M or N is not divisible by the number of dimensions."

        end if

        Mp = M / dims(1)
        Np = N / dims(2)
    
    end subroutine mpi_initialise_standard_topology


    subroutine mpi_send_data(num_dims, to_send, size_to_send, to_recv, size_to_recv, comm)
                
        integer,                          intent(in)    :: num_dims 

        double precision, dimension(:,:), intent(in)    :: to_send
        
        integer,                          intent(in)    :: size_to_send
        
        double precision, dimension(:,:), intent(inout) :: to_recv

        integer,                          intent(in)    :: size_to_recv
        integer,                          intent(in)    :: comm

        integer                                         :: ierr

        call MPI_Scatterv(to_send, counts, displacements, master_type, &
                to_recv, size_to_send, MPI_DOUBLE_PRECISION, 0, cart_comm, ierr)

    end subroutine mpi_send_data


    subroutine mpi_send_halos(num_dims, old, Np, M, dims, nbrs, cart_comm)

        integer,                            intent(in) :: num_dims
        
        double precision, dimension(0:,0:), intent(inout)   :: old

        integer,                          intent(in)        :: Np
        integer,                          intent(in)        :: M
        integer,          dimension(:),   intent(in)        :: nbrs
        integer,                          intent(in)        :: cart_comm

        integer,          dimension(:),   intent(in)        :: dims

        integer, dimension(mpi_status_size)                 :: recv_status      ! Receive status, non-blocking send

        integer                                             :: ierr

        ! Send halos up, receive from down

        call MPI_SENDRECV(old(1, Np), 1, h_halo_type, nbrs(up), 0, &
        old(1,0), 1, h_halo_type, nbrs(down), 0, cart_comm, &
        recv_status, ierr)

        ! Send halos down, receive from up

        call MPI_SENDRECV(old(1, 1), 1, h_halo_type, nbrs(down), 0, &
        old(1, Np+1), 1, h_halo_type, nbrs(up), 0, cart_comm, &
        recv_status, ierr)

        ! Send halos right, receive from left

        call MPI_SENDRECV(old(Mp, 1), 1, v_halo_type, nbrs(right), 0, &
        old(0, 1), 1, v_halo_type, nbrs(left), 0, cart_comm, recv_status, ierr)

        ! Send halos left, receive from right

        call MPI_SENDRECV(old(1,1), 1, v_halo_type, nbrs(left), 0, &
        old(Mp+1, 1), 1, v_halo_type, nbrs(right), 0, cart_comm, recv_status, ierr)


    end subroutine mpi_send_halos


    subroutine mpi_get_average(array, M, N, comm, average)

        double precision, dimension(:,:), intent(in)    :: array 

        integer,                          intent(in)    :: comm
        integer,                          intent(in)    :: M
        integer,                          intent(in)    :: N

        double precision,                 intent(out)   :: average

        integer                                         :: ierr

        double precision                                :: local_sum
        
        local_sum = sum(array)

        call MPI_ALLREDUCE(local_sum, average, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)

        average = average / ( M * N )

    end subroutine mpi_get_average


    subroutine mpi_get_global_delta(num_dims, delta, comm, global_delta)

        integer,                          intent(in)    :: num_dims

        double precision,                 intent(in)    :: delta

        integer,                          intent(in)    :: comm 

        double precision,                 intent(out)   :: global_delta

        integer                                         :: ierr

        if (num_dims .eq. 0) then

            global_delta = delta

        else

            call MPI_ALLREDUCE(delta, global_delta, 1, MPI_DOUBLE_PRECISION, &
                            MPI_MAX, comm, ierr)

        end if

    end subroutine mpi_get_global_delta

    subroutine mpi_gather_data(num_dims, to_send, send_size, to_recv, recv_size, comm)

        integer,                          intent(in)        :: num_dims

        double precision, dimension(:,:), intent(in)        :: to_send

        integer,                          intent(in)        :: send_size
        
        double precision, dimension(:,:), intent(out)       :: to_recv

        integer,                          intent(in)        :: recv_size
        integer,                          intent(in)        :: comm

        integer, dimension(mpi_status_size)                 :: recv_status      ! Receive status, non-blocking send

        integer                                             :: ierr

        call MPI_GATHERV(to_send, 1, subarray_type, to_recv, counts, displacements, master_type, &
                            0, comm, ierr)

    end subroutine mpi_gather_data


end program main

