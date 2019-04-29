module wrappers_2d

    use problem_constants
    use mpi
    use utility_functions
    use neighbour_indexes
    use pgmio     

    implicit none

    PUBLIC

    integer                             :: Mp
    integer                             :: Np

    !
    ! MPI variable initialisation
    !

    integer                             :: master_type
    integer                             :: subarray_type
    integer                             :: v_halo_type
    integer                             :: h_halo_type
    integer                             :: cart_comm

    integer                             :: rank
    integer                             :: pool_size
    integer                             :: ierr

    integer, dimension(8)               :: request

    !
    ! Cartesian topology
    !

    integer, allocatable                :: dims(:)          ! Array of size num_dims: alloc later
    integer, allocatable                :: nbrs(:)        ! Array for addresses of neighbours

    double precision                    :: average

    integer, allocatable                :: counts(:)
    integer, allocatable                :: displacements(:)

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


    subroutine initialise(pool_size, rank, nbrs, dims, cart_comm)

        integer, intent(out)                    :: pool_size
        integer, intent(out)                    :: rank

        integer, dimension(:), intent(inout)    :: nbrs
        integer, dimension(:), intent(inout)    :: dims

        integer,               intent(out)      :: cart_comm

        integer                                 :: i
        integer                                 :: ierr
        integer                                 :: comm

        ! Allocate data arrays

        call MPI_INIT(ierr)

        comm = MPI_COMM_WORLD

        call MPI_COMM_SIZE(comm, pool_size, ierr)

        call mpi_initialise_standard_topology(num_dims, dims, cart_comm, nbrs, rank)

        call util_print_welcomemessage(rank)

        allocate(counts(pool_size))
        allocate(displacements(pool_size))

        call compute_counts_displacements(num_dims, dims, pool_size, Np, counts, displacements)

        call mpi_define_vectors(num_dims, dims, pool_size, Mp, Np)

        call util_print_onrank0("Computing on a two-dimensional grid.", rank)

    end subroutine initialise


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

        ! 28/04 JT: Dimensions shifted to be consistent with Fortran order

        call MPI_CART_CREATE(MPI_COMM_WORLD, num_dims, (/dims(2), dims(1)/), periodic, reorder, cart_comm, ierr)
        
        call MPI_COMM_RANK(cart_comm, rank, ierr)
        
        ! Get neighbours

        call MPI_CART_SHIFT(cart_comm, y_dir, displacement, nbrs(left), nbrs(right), ierr)
        call MPI_CART_SHIFT(cart_comm, x_dir, displacement, nbrs(down), nbrs(up), ierr)
    
        ! Compute the new array dimensions

        if ( (mod(M, dims(1)) .ne. 0) .or. (mod(N, dims(2)) .ne. 0) ) then

            call MPI_FINALIZE(ierr)
            error stop "Error: M or N is not divisible by the number of dimensions."

        end if

        Mp = M / dims(1)
        Np = N / dims(2)
    
    end subroutine mpi_initialise_standard_topology


    subroutine send_data(to_send, size_to_send, to_recv, size_to_recv, comm)
                
        double precision, dimension(:,:), intent(in)    :: to_send
        
        integer,                          intent(in)    :: size_to_send
        
        double precision, dimension(:,:), intent(inout) :: to_recv

        integer,                          intent(in)    :: size_to_recv
        integer,                          intent(in)    :: comm

        integer                                         :: ierr


        call MPI_Scatterv(to_send, counts, displacements, master_type, &
                to_recv, size_to_send, MPI_DOUBLE_PRECISION, 0, cart_comm, ierr)

    end subroutine send_data


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
        integer, dimension(MPI_STATUS_SIZE,8)               :: stats
        ! ! Send halos up, receive from down

        ! call MPI_SENDRECV(old(1, Np), 1, h_halo_type, nbrs(up), 0, &
        ! old(1,0), 1, h_halo_type, nbrs(down), 0, cart_comm, &
        ! recv_status, ierr)

        ! ! Send halos down, receive from up

        ! call MPI_SENDRECV(old(1, 1), 1, h_halo_type, nbrs(down), 0, &
        ! old(1, Np+1), 1, h_halo_type, nbrs(up), 0, cart_comm, &
        ! recv_status, ierr)

        ! ! Send halos right, receive from left

        ! call MPI_SENDRECV(old(Mp, 1), 1, v_halo_type, nbrs(right), 0, &
        ! old(0, 1), 1, v_halo_type, nbrs(left), 0, cart_comm, recv_status, ierr)

        ! ! Send halos left, receive from right

        ! call MPI_SENDRECV(old(1,1), 1, v_halo_type, nbrs(left), 0, &
        ! old(Mp+1, 1), 1, v_halo_type, nbrs(right), 0, cart_comm, recv_status, ierr)

        call MPI_Issend(old(Mp,1),1, v_halo_type, nbrs(right) ,0,cart_comm,request(1),ierr)
        call MPI_Issend(old(1,1),1, v_halo_type, nbrs(left)   ,0,cart_comm,request(3),ierr)
        call MPI_Issend(old(1,1), 1, h_halo_type, nbrs(down),0,cart_comm,request(7),ierr)
        call MPI_Issend(old(1,Np), 1, h_halo_type, nbrs(up) ,0,cart_comm,request(5),ierr)
        
        call MPI_Irecv(old(Mp+1,1),1, v_halo_type, nbrs(right) ,0,cart_comm,request(4),ierr)
        call MPI_Irecv(old(0,1)   ,1, v_halo_type, nbrs(left) ,0,cart_comm,request(2),ierr)
        call MPI_Irecv(old(1,0),1, h_halo_type, nbrs(down), 0,cart_comm,request(6),ierr)
        call MPI_Irecv(old(1,Np+1)   ,1, h_halo_type, nbrs(up) ,0,cart_comm,request(8),ierr)
        
        call MPI_Waitall(8,request,stats,ierr)

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

    subroutine gather_data(num_dims, to_send, send_size, to_recv, recv_size, comm)

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

    end subroutine gather_data


    subroutine check_criterion(new, old, cart_comm, num_dims, num_iters, average, delta_global)

        double precision, dimension(:,:), intent(in)    :: new
        double precision, dimension(:,:), intent(in)    :: old

        integer,                          intent(in)    :: cart_comm
        integer,                          intent(in)    :: num_dims
        integer,                          intent(in)    :: num_iters

        double precision,                 intent(out)   :: average
        double precision,                 intent(out)   :: delta_global

        double precision                                :: local_delta

        call mpi_get_average(new, M, N, cart_comm, average)
        call util_get_local_delta(new, old, local_delta)
        call mpi_get_global_delta(num_dims, local_delta, cart_comm, delta_global)
        call util_print_average_max(average, delta_global, num_iters, rank)

    end subroutine check_criterion


    subroutine util_printfinish(num_iters, time_start, time_finish, rank)

        integer,            intent(in)  :: num_iters
        
        double precision,   intent(in)  :: time_start
        double precision,   intent(in)  :: time_finish

        integer,            intent(in)  :: rank

        integer                         :: ierr

        if (num_iters .eq. max_iters) then

            call MPI_FINALIZE(ierr)
            error stop "Maximum number of iterations reached. Stopping..."
        
        end if

        if (rank .eq. 0) then

            print '(A)',          "  ****************************************************************"
            print '(A, F8.4, A)', "             Finished computing in ", time_finish - time_start, " seconds."
            print '(A)',          "  ****************************************************************"

        end if

    end subroutine

    subroutine finalise()

        integer :: ierr

        ! Free the types we created

        call MPI_TYPE_FREE(master_type, ierr)
        call MPI_TYPE_FREE(subarray_type, ierr)
        call MPI_TYPE_FREE(h_halo_type, ierr)
        call MPI_TYPE_FREE(v_halo_type, ierr)

        call MPI_FINALIZE(ierr)

    end subroutine finalise


end module wrappers_2d