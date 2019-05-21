program main

    use mpi

    implicit none

    integer                             :: rank, pool_size, comm, request
    integer                             :: i, j

    integer                             :: left, right
    integer                             :: to_send, local_count, recvd_data
    integer                             :: ierr

    double precision                    :: time1, time2, tot_time

    integer, dimension(MPI_STATUS_SIZE) :: stat, send_stat

    !
    ! Data variables for cartesian topology stuff
    !


    ! Array of dimension lengths (only one here)
    integer                             :: dims(1)
    ! Number of dimensions
    integer                             :: ndims
    ! Logicals: periodic BCs? Reorder the grid?
    logical                             :: periods(1), reorder
    ! New pool for topology
    integer                             :: new_comm
    ! Number of steps l/r
    integer                             :: displacement
    call MPI_INIT(ierr)

    ! Define the communicator

    comm = MPI_COMM_WORLD

    ! Get the size and rank of the pool

    call MPI_COMM_RANK(comm, rank, ierr)
    call MPI_COMM_SIZE(comm, pool_size, ierr)

    ! Initialise our dimensions array

    dims(1) = 0 ! 0 => mpi should fill in the value properly

    ! Initialise other cartesian problem variables

    ndims = 1
    periods(1) = .true. ! Want periodic BCs
    reorder = .false. ! Don't change things to make them more efficient
    displacement = 1

    ! Create our dimensions

    call MPI_DIMS_CREATE(pool_size, ndims, dims, ierr)

    ! Create a new pool for our cartesian topology: attaches topology information to the grid
    ! Our old communicator is just MPI_COMM_WORLD

    call MPI_CART_CREATE(comm, ndims, dims, periods, reorder, new_comm, ierr)

    ! Get the rank in the new communicator

    call MPI_COMM_RANK(new_comm, rank, ierr)

    call MPI_CART_SHIFT(new_comm, 0, displacement, left, right, ierr)

    ! Now go through and do the whole message passing thing

    to_send = rank

    call MPI_BARRIER(new_comm, ierr)

    time1 = mpi_wtime()

    do i = 1, pool_size

        call MPI_ISsend(to_send, 1, MPI_INTEGER, right, 0, new_comm, &
            request, ierr)

        call MPI_Recv(recvd_data, 1, MPI_INTEGER, left, 0, new_comm, &
            stat, ierr)

        call MPI_WAIT(request, send_stat, ierr)

        local_count = local_count + recvd_data

        to_send = recvd_data

    end do

    call MPI_BARRIER(new_comm, ierr)

    time2 = mpi_wtime()

    print *, "I am rank", rank, "and the sum is", local_count

    if (rank .eq. 0) then

        print *, "Execution time:", time2-time1

    end if

    call MPI_FINALIZE(ierr)

end program main
