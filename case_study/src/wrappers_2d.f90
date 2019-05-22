module wrappers_2d

    ! Modules

    use problem_constants                                           ! Defines fixed problem variables
    use mpi                                                         ! Provides distributed parallelism
    use utility_functions                                           ! Provides common non-MPI subroutines
    use neighbour_indexes                                           ! Provides standard nbrs array indices
    use pgmio                                                       ! Provides IO routines for PGM files

    implicit none

    PUBLIC

    integer                             :: Mp                       ! Horizontal chunk per processor
    integer                             :: Np                       ! Vertical chunk per processor

    !
    ! MPI variable initialisation
    !

    integer                             :: master_type              ! New MPI type to read from M x N data array
    integer                             :: subarray_type            ! MPI type to read to Mp x Np data array (resized)
    integer                             :: v_halo_type              ! MPI vector to send vertical halos
    integer                             :: h_halo_type              ! MPI vector to send horizontal halos
    integer                             :: cart_comm                ! New commmunicator for cartesian topology

    integer                             :: rank                     ! Processor rank
    integer                             :: pool_size                ! Size of the worker pool
    integer                             :: ierr                     ! Integer error code

    integer, dimension(8)               :: request                  ! Send/receive request status

    !
    ! Cartesian topology
    !

    integer, allocatable                :: dims(:)                  ! Array holding topology dimensions: 2 * num_dims
    integer, allocatable                :: nbrs(:)                  ! Array for addresses of neighbours: 2 * dims
    integer, allocatable                :: counts(:)                ! Array to hold the number of elements for the resized subarray
    integer, allocatable                :: displacements(:)         ! Array to hold the displacements for the resized subarray

    double precision                    :: average                  ! Average value of the image data array

    contains


    subroutine mpi_define_vectors(num_dims, dims, Mp, Np)

        !
        ! MPI_DEFINE_VECTORS
        ! ~~~~~~~~~~~~~~~~~~
        !
        ! Creates and commits the supplementary MPI types used in the program.
        !
        ! Inputs
        ! ~~~~~~
        ! num_dims: Number of dimensions in the problem. Type: integer.
        ! dims: Assigned dimensions.                     Type: 1-d integer array.
        ! pool_size: Number of workers in the pool.      Type: integer.
        ! Mp: Horizontal chunk size.                     Type: integer.
        ! Np: Vertical chunk size.                       Type: integer.
        !

        integer,                    intent(in)  :: num_dims
        integer,                    intent(in)  :: Mp
        integer,                    intent(in)  :: Np

        integer, dimension(:),      intent(in)  :: dims

        integer, allocatable                    :: sizes(:)
        integer, allocatable                    :: subsizes(:)
        integer, allocatable                    :: starts(:)

        integer                                 :: ierr
        integer                                 :: temp_type

        integer(kind=mpi_address_kind)          :: start
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
    
        ! Create the temporary master block, whih provides an array to convert directly from the master buffer
        ! to the working thread. This will soon be resized into the true master block.

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

        !
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

        ! Deallocate arrays we're no longer usingcounts

        deallocate(sizes)
        deallocate(subsizes)
        deallocate(starts)

    end subroutine mpi_define_vectors

    subroutine compute_counts_displacements(dims, pool_size, Np, counts, displacements)
        
        !
        ! COMPUTE_COUNTS_DISPLACEMENTS
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        !
        ! Initialises the counts and displacements arrays for use in 
        ! SCATTERV and GATHERV subroutines.
        !
        ! Inputs
        ! ~~~~~~
        ! num_dims: Number of dimensions in the problem. Type: integer.
        ! dims: Assigned dimensions.                     Type: 1-d integer array.
        ! pool_size: Number of workers in the pool.      Type: integer.
        ! Np: Vertical chunk size.                       Type: integer.
        !
        ! Outputs
        ! ~~~~~~
        ! counts: Specifies the number of elements to
        !         send to each processor.                Type: 1-d integer array.
        ! displacements: displacement from which to take
        !                the outgoing data to process.   Type: 1-d integer array.
        !

        integer, dimension(:), intent(in)       :: dims
        
        integer,               intent(in)       :: pool_size
        integer,               intent(in)       :: Np

        integer, dimension(:), intent(inout)    :: counts
        integer, dimension(:), intent(inout)    :: displacements

        integer                                 :: baseline
        integer                                 :: i

        baseline = 1

        do i = 1, pool_size

                ! Counts is always one

                counts(i) = 1

                displacements(i) = (baseline-1) + mod(i-1, dims(1))

                if (mod(i, dims(1)) .eq. 0) then
                
                    baseline = baseline + Np * dims(1)

                end if

        end do

    end subroutine compute_counts_displacements


    subroutine initialise(pool_size, rank, nbrs, dims, cart_comm)

        !
        ! INITIALISE prepares the environment for the MPI runtime:
        ! it checks the size of the pool is as expected, initialises
        ! a cartesian topolog, computes the counts and displacements and
        ! initialises the vectors and custom types to be used.
        !
        ! Inputs & Outputs
        ! ~~~~~~~~~~~~~~~~
        ! nbrs: Array containing the locations of 
        !       each processors neigbours.          Type: 1-d integer array.
        ! dims: Array containing the dimensions of
        !       the new cartesian topology.         Type: 1-d integer array.     
        ! 
        ! Outputs
        ! ~~~~~~~
        ! pool_size: Number of workers in the pool. Type: integer.
        ! rank: identifier of the current
        !       processor.                          Type: integer.
        ! cart_comm: New communicator for the
        !            cartesian topology.            Type: integer
        !

        integer, intent(out)                    :: pool_size
        integer, intent(out)                    :: rank

        integer, dimension(:), intent(inout)    :: nbrs
        integer, dimension(:), intent(inout)    :: dims

        integer,               intent(out)      :: cart_comm

        integer                                 :: comm                             ! Old communicator

        ! Initialise MPI

        call MPI_INIT(ierr)

        comm = MPI_COMM_WORLD

        ! Get the size of the pool

        call MPI_COMM_SIZE(comm, pool_size, ierr)

        ! Check that we have the correct number of threads

        if (pool_size .ne. P) then

            call MPI_FINALIZE(ierr)
            error stop "Incorrect number of processors specified."

        end if


        ! Initialise a 2-d cartesian topology

        call mpi_initialise_standard_topology(num_dims, dims, cart_comm, nbrs, rank)

        ! Print a welcome message

        call util_print_welcomemessage(rank)

        ! Allocate arrays used in counts and displacements.

        allocate(counts(pool_size))
        allocate(displacements(pool_size))

        call compute_counts_displacements(dims, pool_size, Np, counts, &
                                         displacements)

        ! Define our custom mpi datatypes

        call mpi_define_vectors(num_dims, dims, Mp, Np)

        ! Print a wonderfully informative message.

        call util_print_onrank0("Computing on a two-dimensional grid.", rank)

    end subroutine initialise


    subroutine mpi_initialise_standard_topology(num_dims, dims, cart_comm, nbrs, rank)

        ! MPI_INITIALISE_STANDARD_TOPOLOGY
        !
        ! Creates a 2D cartesian communicator with the appropriate dimensions,
        ! computes the neighbours for each process and computes MP and NP for
        ! this problem.
        !
        ! Inputs
        ! ~~~~~~
        ! num_dims: Number of dimensions in the problem. Type: integer
        !
        ! Inputs & Outputs
        ! ~~~~~~~~~~~~~~~~
        ! dims: array containing the length in each 
        !       dimension.                               Type: 1-d integer array
        ! nbrs: Array containing the locations of 
        !       each processors neigbours.               Type: 1-d integer array.
        !
        ! Outputs
        ! ~~~~~~~
        ! cart_comm: New communicator for the cartesian
        !            topology.                           Type: integer
        ! rank: Identifier for this processor.           Type: integer
        
        integer,               intent(in)       :: num_dims

        integer, dimension(:), intent(inout)    :: dims

        integer,               intent(out)      :: cart_comm

        integer, dimension(:), intent(inout)    :: nbrs

        integer,               intent(out)      :: rank

        logical                                 :: periodic(2)                          ! Tracks whether we want periodic boundaries              
        logical                                 :: reorder                              ! Tracks whether we want to re-order ranks for speed

        integer                                 :: x_dir                                ! x-direction in our cartesian shift
        integer                                 :: y_dir                                ! y-direction in our cartesian shift
        integer                                 :: displacement                         ! displacement: step to take in cartesian shift
        integer                                 :: i                                    ! Loop counter
        integer                                 :: pool_size                            ! Number of workers in the pool
        
        ! Get the size of the pool

        call MPI_COMM_SIZE(MPI_COMM_WORLD, pool_size, ierr)

        ! Initialise dims to zero so that it may be populated, and periodic to false
        ! so we have non-periodic boundary conditions on the nodes

        do i = 1, P

            dims(i) = 0
            periodic(i) = .false.

        end do

        ! Initialise the remainder of our variables

        reorder      = .false.
        x_dir        = 0
        y_dir        = 1
        displacement = 1

        ! Create 2D cartesian topology

        call MPI_DIMS_CREATE(P, num_dims, dims, ierr)

        ! Create the cartesian topology: dimensions are shifted to be consistent
        ! with Fortran array ordering.

        call MPI_CART_CREATE(MPI_COMM_WORLD, num_dims, (/dims(2), dims(1)/), &
                             periodic, reorder, cart_comm, ierr)
        
        ! Get this processor's rank

        call MPI_COMM_RANK(cart_comm, rank, ierr)
        
        ! Get neighbours - up and down

        call MPI_CART_SHIFT(cart_comm, y_dir, displacement, nbrs(left), nbrs(right), ierr)
        call MPI_CART_SHIFT(cart_comm, x_dir, displacement, nbrs(down), nbrs(up), ierr)
    
        ! Compute the new array dimensions

        print *, "the dimensions are:", dims(1), dims(2)


        if ( (mod(M, dims(1)) .ne. 0) .or. (mod(N, dims(2)) .ne. 0) ) then

            call MPI_FINALIZE(ierr)
            print *, "Error: M or N is not divisible by the number of dimensions."
            error stop "Error: M or N is not divisible by the number of dimensions."

        end if

        Mp = M / dims(1)
        Np = N / dims(2)
    
    end subroutine mpi_initialise_standard_topology


    subroutine send_data(to_send, size_to_send, to_recv, comm)

        !
        ! send_data controls the initial transfer of data from the master buffer
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

        integer                                         :: ierr

        ! Call SCATTERV: assumes counts and displacements have been previously defined
        ! and populated.    

        call MPI_Scatterv(to_send, counts, displacements, master_type, &
                to_recv, size_to_send, MPI_DOUBLE_PRECISION, 0, comm, ierr)

    end subroutine send_data


    subroutine mpi_send_halos(old, nbrs, cart_comm)

        !
        ! mpi_send_halos transfers the halo regions of the image between each
        ! processor
        !
        ! Inputs:
        ! ~~~~~~~
        ! Np: Height of image chunk. Type: integer
        ! M: Width of original image. Type: integer
        ! nbrs: Array of cartesian neighbours. Type: 1-d integer array
        ! cart_comm: Cartesian communicator to use for the transfers. Type: integer
        ! dims: Array of dimensions of the cartesian topology. Type: 1-d integer array.
        !
        ! Inputs & Outputs
        ! ~~~~~~~~~~~~~~~~
        ! old: Array to send data to and from. Type: 2-d double precision, indexed from 0
        !

        double precision, dimension(0:,0:), intent(inout)   :: old

        integer,          dimension(:),   intent(in)        :: nbrs
        integer,                          intent(in)        :: cart_comm

        integer                                             :: ierr
        integer, dimension(MPI_STATUS_SIZE,8)               :: stats

        integer, dimension(MPI_STATUS_SIZE)                 :: stat

        ! Send right

        call MPI_Issend(old(Mp, 1), 1, v_halo_type, nbrs(right),  0, cart_comm, request(1), ierr)
        
        ! Send left

        call MPI_Issend(old(1, 1),  1, v_halo_type, nbrs(left),   0, cart_comm, request(3), ierr)
        
        ! Send down

        call MPI_Issend(old(1, 1),  1, h_halo_type, nbrs(down),   0, cart_comm, request(7), ierr)
        
        ! Send up

        call MPI_Issend(old(1, Np), 1, h_halo_type, nbrs(up),     0, cart_comm, request(5), ierr)

        ! Receive from right

        call MPI_Irecv(old(Mp+1, 1), 1, v_halo_type, nbrs(right), 0, cart_comm, request(4), ierr)
        
        ! Receive from left
        
        call MPI_Irecv(old(0, 1),    1, v_halo_type, nbrs(left),  0, cart_comm, request(2), ierr)
       
        ! Receive from below

        call MPI_Irecv(old(1, 0),    1, h_halo_type, nbrs(down),  0, cart_comm, request(6), ierr)
        
        ! Receive from above
        
        call MPI_Irecv(old(1, Np+1), 1, h_halo_type, nbrs(up),    0, cart_comm, request(8), ierr)
        
        ! Wait for all the transfers to complete

        call MPI_Waitall(8, request, stats, ierr)

    end subroutine mpi_send_halos


    subroutine mpi_get_average(array, M, N, comm, average)

        ! mpi_get_average computes the average value of the elements in array
        !
        ! Inputs
        ! ~~~~~~
        ! array: Array to compute the average of elements in. Type: 2-d array of double precision
        ! comm: Cartesian communicator to reduce the data from
        ! M: Width of original image. Type: integer
        ! N: Height of original image. Type: integer
        ! 
        ! Outputs
        ! ~~~~~~~
        ! Average: average value of elements in array. Type: double precision.
        !

        double precision, dimension(:,:), intent(in)    :: array 

        integer,                          intent(in)    :: comm
        integer,                          intent(in)    :: M
        integer,                          intent(in)    :: N

        double precision,                 intent(out)   :: average

        integer                                         :: ierr

        double precision                                :: local_sum
        

        ! Get the local sum of all the values in the local versions of arrays

        local_sum = sum(array)

        ! Reduce all the local sums onto every processor, taking the sum of every value

        call MPI_ALLREDUCE(local_sum, average, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)

        ! Normalise to the size of the matrix: integer division should be promoted to double precision by all
        ! compilers

        average = average / ( M * N )

    end subroutine mpi_get_average


    subroutine mpi_get_global_delta(delta, comm, global_delta)

        !
        ! mpi_get_global_delta computes the maximum change between two iterations
        ! of the algorithm across all processors.
        !
        ! Inputs
        ! ~~~~~~
        !
        ! delta: Local value of the change between two iterations. Type: double precision
        ! comm: Communicator to reduce on. Type: integer
        !
        ! Outputs
        ! ~~~~~~~
        ! global_delta: Maximum delta across all processors. Type: double precision
        !

        double precision,                 intent(in)    :: delta

        integer,                          intent(in)    :: comm 

        double precision,                 intent(out)   :: global_delta

        integer                                         :: ierr

        call MPI_ALLREDUCE(delta, global_delta, 1, MPI_DOUBLE_PRECISION, &
                            MPI_MAX, comm, ierr)

    end subroutine mpi_get_global_delta


    subroutine gather_data(to_send, to_recv, comm)

        !
        ! gather_data retrieves the data from all other processors at the end of
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

        integer                                             :: ierr

        call MPI_GATHERV(to_send, 1, subarray_type, to_recv, counts, displacements, master_type, &
                            0, comm, ierr)

    end subroutine gather_data


    subroutine check_criterion(new, old, cart_comm, num_iters, average, delta_global)

        !
        ! check_criterion is a wrapping function that allows the determination of the global
        ! delta, local delta, average value of an array, and also prints a status message.
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

        double precision                                :: local_delta

        call mpi_get_average(new, M, N, cart_comm, average)
        call util_get_local_delta(new, old, local_delta)
        call mpi_get_global_delta(local_delta, cart_comm, delta_global)
        call util_print_average_max(average, delta_global, num_iters, rank)

    end subroutine check_criterion


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

        !
        ! finalise completes the MPI runtime and frees the created types.
        !

        integer :: ierr

        ! Free the types we created

        call MPI_TYPE_FREE(master_type, ierr)
        call MPI_TYPE_FREE(subarray_type, ierr)
        call MPI_TYPE_FREE(h_halo_type, ierr)
        call MPI_TYPE_FREE(v_halo_type, ierr)

        call MPI_FINALIZE(ierr)

    end subroutine finalise


end module wrappers_2d
