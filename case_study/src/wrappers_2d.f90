module wrappers_2d

    use problem_constants
    use mpi
    use utility_functions
    use neighbour_indexes
    use pgmio     

    implicit none

    integer :: Mp
    integer :: Np

    integer :: master_type
    integer :: subarray_type
    integer :: v_halo_type
    integer :: h_halo_type
    
    PUBLIC :: Mp, Np

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

        integer, intent(in)     :: num_dims

        integer, intent(out)    :: pool_size
        integer, intent(out)    :: rank

        integer, dimension(:), intent(inout)  :: nbrs
        integer, dimension(:), intent(inout)  :: dims

        integer,               intent(out)    :: cart_comm

        integer                 :: i
        integer                 :: ierr
        integer                 :: comm

        call MPI_INIT(ierr)

        comm = MPI_COMM_WORLD

        call get_comm_size(comm, pool_size)

        allocate(counts(pool_size))

        if ( (mod(M, dims(1)) .ne 0) .or. (mod(N, dims(2)) .ne. 0) ) then

            call MPI_FINALIZE(ierr)
            error stop "Error: M or N is not divisible by the number of dimensions."

        end if

        Mp = M/dims(1)
        Np = N/dims(2)

        call compute_counts_displacements(num_dims, dims, pool_size, Np, counts, displacements)
        call mpi_define_vectors(num_dims, dims, pool_size, Mp, Np)

    end subroutine mpi_initialise

end module wrappers_2d