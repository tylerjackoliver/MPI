program main

    use mpi

    implicit none

    integer                         :: comm, rank, pool_size
    integer                         :: ierr
    integer                         :: left, right
    integer, &
        dimension(MPI_STATUS_SIZE)  :: stat, send_stat

    integer                         :: i,j

    integer                         :: local_count
    integer                         :: recvd_data

    integer                         :: request
    integer                         :: to_send

    double precision                :: time1, time2, tot_time

    ! Initialise MPI

    call MPI_INIT(ierr)

    comm = MPI_COMM_WORLD

    ! Get the size of the pool and rank of the proc

    call MPI_COMM_SIZE(MPI_COMM_WORLD, pool_size, ierr)

    ! Get the rank for each processor

    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

    ! Initialise ring counter

    local_count = 0
    recvd_data = 0

    to_send = (rank + 1) ** 2

    ! Allocate array that contains the address of the destination

    left = rank - 1

    if (rank .eq. 0) then

        left = pool_size - 1

    end if

    right = rank + 1

    if ((rank + 1) .eq. pool_size) then

        right = 0

    end if

    ! Go round the ring P-1 times

    call MPI_BARRIER(comm, ierr)

    time1 = mpi_wtime()

    do i = 1,pool_size

            call MPI_Sendrecv(to_send, 1, MPI_INTEGER, &
                right, 0, recvd_data, 1, MPI_INTEGER, &
                left, 0, comm, stat, ierr)

            local_count = local_count + recvd_data

            to_send = recvd_data

    end do

    call MPI_BARRIER(comm, ierr)

    time2 = mpi_wtime()

    if (rank .eq. 0) then

        print *, "Time required: ", time2-time1

    end if

    ! Print all the messages

    ! Finalise MPI

    call MPI_FINALIZE(ierr)

end program main
