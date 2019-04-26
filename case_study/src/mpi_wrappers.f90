module mpi_wrappers

    use mpi 
    use problem_constants
    use neighbour_indexes

    implicit none

    ! MPI implementation constants

    integer             :: master_type
    integer             :: subarray_type
    integer             :: v_halo_type
    integer             :: h_halo_type
    integer             :: Mp               ! Horizontal image slice
    integer             :: Np               ! Vertical image slice
    
    integer             :: comm
    integer             :: cart_comm

    contains

    subroutine mpi_initialise_standard_topology(num_dims, pool_size, old_comm, new_comm, new_rank, dims, nbrs)
        
        use neighbour_indexes                   ! Provides standardised values of left, right

        integer, intent(in)                     :: num_dims
        integer, intent(in)                     :: pool_size
        integer, intent(in)                     :: old_comm

        integer, intent(out)                    :: new_comm
        integer, intent(out)                    :: new_rank
        integer, intent(out), allocatable       :: dims(:)

        integer, intent(inout), dimension(:)    :: nbrs

        integer                                 :: ierr
        integer                                 :: x_dir 
        integer                                 :: displacement
        integer                                 :: i

        logical                                 :: periodic(1)
        logical                                 :: reorder(1)

        allocate(dims(num_dims))
        reorder = .false.

        if (num_dims .eq. 0) then   ! running serially

            new_rank = 0            ! Fake that we're on rank 0
            new_comm = -1           ! Makes sure we get an error if we try to MPI
            nbrs     = 0            ! Makes sure we get an error if we try to MPI

        elseif (num_dims .eq. 1) then

            call initialise_standard_topology_1d(pool_size, old_comm, new_comm, new_rank, nbrs)
        
        elseif (num_dims .eq. 2) then

            call initialise_standard_topology_2d(pool_size, old_comm, new_comm, new_rank, dims, nbrs)

        end if

    end subroutine mpi_initialise_standard_topology


    subroutine initialise_standard_topology_1d(pool_size, old_comm, new_comm, new_rank, nbrs)
        
        use neighbour_indexes                   ! Provides standardised values of left, right

        integer, intent(in)                     :: pool_size
        integer, intent(in)                     :: old_comm

        integer, intent(out)                    :: new_comm
        integer, intent(out)                    :: new_rank
        integer, intent(inout), dimension(:)    :: nbrs

        integer, parameter                      :: num_dims = 1
        integer, intent(in), dimension(:)       :: dims
        integer                                 :: ierr
        integer                                 :: x_dir 
        integer                                 :: displacement
        integer                                 :: i

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
        call MPI_CART_CREATE(old_comm, num_dims, dims, periodic, reorder, new_comm, ierr)
        call MPI_COMM_RANK(new_comm, new_rank, ierr)
        call MPI_CART_SHIFT(new_comm, x_dir, displacement, nbrs(left), nbrs(right), ierr)

    end subroutine initialise_standard_topology_1d


    subroutine initialise_standard_topology_2d(pool_size, old_comm, new_comm, new_rank, dims, nbrs)
        
        use neighbour_indexes                   ! Provides standardised values of left, right

        integer, intent(in)                     :: pool_size
        integer, intent(in)                     :: old_comm

        integer, intent(out)                    :: new_comm
        integer, intent(out)                    :: new_rank
        integer, intent(out)                    :: dims(2)

        integer, intent(inout), dimension(:)    :: nbrs

        integer, parameter                      :: num_dims = 2
        integer                                 :: ierr
        integer                                 :: x_dir 
        integer                                 :: y_dir
        integer                                 :: displacement_x
        integer                                 :: displacement_y
        integer                                 :: i

        logical                                 :: periodic(2)
        logical                                 :: reorder

        ! Initialise dimensions, periodic and reorder

        do i = 1, num_dims

            dims(i) = 0
            periodic(i) = .false.

        end do

        reorder = .false.
    
        x_dir = 0
        y_dir = 1

        displacement_x = 1
        displacement_y = 1

        ! For 2-d topology, nbrs is a 1x4 array, so need to call
        ! the cartesian topology twice

        call MPI_DIMS_CREATE(pool_size, num_dims, dims, ierr)
        call MPI_CART_CREATE(old_comm, num_dims, dims, periodic, reorder, new_comm, ierr)
        call MPI_COMM_RANK(new_comm, new_rank, ierr)
        call MPI_CART_SHIFT(new_comm, x_dir, displacement_x, nbrs(left), nbrs(right), ierr)
        call MPI_CART_SHIFT(new_comm, y_dir, displacement_y, nbrs(down), nbrs(up), ierr)

    end subroutine initialise_standard_topology_2d


    subroutine mpi_send_data(num_dims, to_send, size_to_send, to_recv, size_to_recv, comm)
        
        integer,                          intent(in)    :: num_dims 

        double precision, dimension(:,:), intent(in)    :: to_send
        
        integer,                          intent(in)    :: size_to_send
        
        double precision, dimension(:,:), intent(inout) :: to_recv

        integer,                          intent(in)    :: size_to_recv
        integer,                          intent(in)    :: comm

        integer                                         :: ierr

        if (num_dims .eq. 0) then

            to_recv = to_send

        else if (num_dims .eq. 1) then

            call mpi_send_data_1d(to_send, size_to_send, to_recv, size_to_recv, comm)

        elseif (num_dims .eq. 2) then

            call mpi_send_data_2d(to_send, size_to_send, to_recv, size_to_recv, comm)

        end if

    end subroutine


    subroutine mpi_send_data_1d(to_send, size_to_send, to_recv, size_to_recv, comm)

        double precision, dimension(:,:), intent(in)    :: to_send
        
        integer,                          intent(in)    :: size_to_send
        
        double precision, dimension(:,:), intent(inout) :: to_recv

        integer,                          intent(in)    :: size_to_recv
        integer,                          intent(in)    :: comm

        integer                                         :: ierr

        call MPI_SCATTER(to_send, size_to_send, MPI_DOUBLE_PRECISION, to_recv, &
                        size_to_recv, MPI_DOUBLE_PRECISION, 0, comm, ierr)

    end subroutine mpi_send_data_1d


    subroutine mpi_send_data_2d(to_send, size_to_send, to_recv, size_to_recv, comm)

        double precision, dimension(:,:), intent(in)    :: to_send
        
        integer,                          intent(in)    :: size_to_send
        
        double precision, dimension(:,:), intent(inout) :: to_recv

        integer,                          intent(in)    :: size_to_recv
        integer,                          intent(in)    :: comm

        integer                                         :: ierr

        call MPI_Scatterv(source, counts, displs, MASTER_BLOCK_T, &
                          dest, Mp*Np, MPI_REALNUMBER, 0, cartcomm,ierr)

    end subroutine mpi_send_data_2d

    
    subroutine mpi_send_halos(num_dim, old, Np, M, dims, nbrs, cart_comm)

        integer,                            intent(in)      :: num_dim

        double precision, dimension(0:,0:), intent(inout)   :: old

        integer,                            intent(in)      :: Np
        integer,                            intent(in)      :: M
        integer,          dimension(:),     intent(in)      :: dims
        integer,          dimension(:),     intent(in)      :: nbrs
        integer,                            intent(in)      :: cart_comm

        if (num_dim .eq. 0) then

            ! Do nothing - old is already what it should be

        elseif (num_dim .eq. 1) then

            call send_halos_1d(old, Np, M, nbrs, cart_comm)

        elseif (num_dim .eq. 2) then

            call send_halos_2d(old, Np, M, dims, nbrs, cart_comm)

        end if

    end subroutine mpi_send_halos


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


    subroutine send_halos_2d(old, Np, M, dims, nbrs, cart_comm, buf)

        double precision, dimension(0:,0:), intent(inout)   :: old

        integer,                          intent(in)        :: Np
        integer,                          intent(in)        :: M
        integer,          dimension(:),   intent(in)        :: nbrs
        integer,                          intent(in)        :: cart_comm

        integer, dimension(mpi_status_size)                 :: recv_status      ! Receive status, non-blocking send

        integer                                             :: ierr
        integer                                             :: pool_size
        integer                                             :: num_dims

        ! For 2d, num_dims = 2

        num_dims = 2

        call get_comm_size(cart_comm, pool_size)

        call MPI_Scatterv(old, counts, displacements, master_type, &
                          buf, Mp*Np, MPI_DOUBLE_PRECISION, 0, cart_comm, ierr)


    end subroutine send_halos_2d


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


    ! subroutine get_array_element_1d(new, old, edge, i, j)
        
    !     double precision, dimension(:,:), intent(inout) :: new

    !     double precision, dimension(:,:), intent(in)    :: old
    !     double precision, dimension(:,:), intent(in)    :: edge

    !     integer,                          intent(in)    :: i
    !     integer,                          intent(in)    :: j

    !     double precision                                :: frac = .25d0
        
    !     double precision                                :: temp_sum
        
    !     temp_sum = old(i-1, j) + old(i+1, j) + old(i, j+1) +  old(i, j-1)

    !     new(i, j) = frac * ( temp_sum - edge(i, j) )

    ! end subroutine get_array_element_1d


    subroutine mpi_get_global_delta(num_dims, delta, comm, global_delta)

        integer,                          intent(in)    :: num_dims
    
        double precision,                 intent(in)    :: delta
    
        integer,                          intent(in)    :: comm 
    
        double precision,                 intent(out)   :: global_delta
    
        integer                                         :: ierr
    
        if (num_dims .eq. 0) then

            global_delta = delta

        elseif (num_dims .eq. 1) then

            call MPI_ALLREDUCE(delta, global_delta, 1, MPI_DOUBLE_PRECISION, &
                            MPI_MAX, comm, ierr)

        elseif (num_dims .eq. 2) then

            print *, "TBD."

        end if
    
    end subroutine mpi_get_global_delta


    subroutine util_get_average(array, M, N, comm, average)

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

    end subroutine util_get_average
        

    subroutine mpi_gather_data(num_dims, to_send, send_size, to_recv, recv_size, comm)

        integer,                          intent(in)        :: num_dims

        double precision, dimension(:,:), intent(in)        :: to_send

        integer,                          intent(in)        :: send_size
        
        double precision, dimension(:,:), intent(out)       :: to_recv

        integer,                          intent(in)        :: recv_size
        integer,                          intent(in)        :: comm

        integer, dimension(mpi_status_size)                 :: recv_status      ! Receive status, non-blocking send

        integer                                             :: ierr

        if (num_dims .eq. 0) then ! Serial case

            to_recv = to_send

        end if

        if (num_dims .eq. 1) then

            call MPI_GATHER(to_send, send_size, MPI_DOUBLE_PRECISION, to_recv, recv_size,&
                        MPI_DOUBLE_PRECISION, 0, comm, ierr)

        elseif (num_dims .eq. 2) then

            print *, "Pass"

        end if

    end subroutine mpi_gather_data


    subroutine mpi_define_vectors(num_dims, dims, pool_size, Mp, Np)

        integer,                    intent(in)  :: num_dims
        integer,                    intent(in)  :: pool_size
        integer,                    intent(in)  :: Mp
        integer,                    intent(in)  :: Np

        integer, dimension(:),      intent(in)  :: dims

        integer,                    intent(out) :: subarray_type
        integer,                    intent(out) :: master_type
        integer,                    intent(out) :: v_halo_type
        integer,                    intent(out) :: h_halo_type

        integer, allocatable,       intent(out) :: cnts
        integer, allocatable,       intent(out) :: dsplcmnts

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
        extent = Mp * double_extent
        
        call MPI_TYPE_CREATE_RESIZED(long_type, start, extent, &
                                     master_type,ierr)

        ! Commit these new datatypes so that we can use them

        call MPI_TYPE_COMMIT(master_type, ierr)
        call MPI_TYPE_COMMIT(subarray_type, ierr)
        call MPI_TYPE_COMMIT(v_halo_type, ierr)
        call MPI_TYPE_COMMIT(h_halo_type, ierr)

        deallocate(sizes(num_dims))
        deallocate(subsizes(num_dims))
        deallocate(starts(num_dims))

    end subroutine mpi_define_vectors


    subroutine compute_counts_displacements(num_dims, dims, pool_size, Np, counts, displacements)

        integer,               intent(in)       :: num_dims

        integer, dimension(:), intent(in)       :: dims
        
        integer,               intent(in)       :: pool_size
        integer,               intent(in)       :: Np

        integer, allocatable                    :: counts(:)
        integer, allocatable                    :: displacements(:)

        integer                                 :: baseline

        allocate(cnts(num_dims))
        allocate(dsplcmnts(num_dims))

        baseline = 1

        do i = 1, pool_size

                counts(i) = 1
                displs(i) = (baseline-1) + mod(i-1, dims(1))

                if (mod(i, dims(1)) .eq. 0) then
                
                    baseline = baseline + Np * dims(1)

                end if

        end do

    end subroutine compute_counts_displacements


    subroutine get_comm_size(comm, size)

        integer, intent(in) :: comm
        integer, intent(out):: size 

        integer             :: ierr

        call MPI_COMM_SIZE(comm, size, ierr)

    end subroutine get_comm_size


    subroutine mpi_initialise(num_dims, size, rank)

        integer, intent(in)     :: num_dims

        integer, intent(out)    :: size
        integer, intent(out)    :: rank

        integer                 :: ierr

        call MPI_INIT(ierr)

        comm = MPI_COMM_WORLD

        call get_comm_size(comm, size)

        ! If size =/ P, then exit
    
        if (pool_size .ne. P) then
    
            call MPI_FINALIZE(ierr)
            error stop "Number of processors does not match the expected number of processes."
    
        end if
        
        ! Define local problem array sizes

        if (rank .eq. 0) then

            Mp = M 
            Np = N

        elseif (rank .eq. 1) then

            Mp = M
            Np = ceiling(dble(N/P))                                ! Assumes N/P is perfect

        elseif (rank .eq. 2) then

            Mp = M/dims(1)
            Np = N/dims(2)

        end if


        call util_print_welcomemessage(rank)
        call mpi_initialise_standard_topology(num_dims, size, comm, cart_comm, rank, dims, nbrs)

        call compute_counts_displacements(num_dims, dims, size, Np, counts, displacements)
        call mpi_define_vectors(num_dims, dims, pool_size, Mp, Np)

end module mpi_wrappers

