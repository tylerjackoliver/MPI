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

    integer                         :: sizeof_int

    double precision                :: time1, time2, tot_time

    !
    ! Data buffers
    !

    integer, parameter              :: bufsize = 1000

    character                       :: buffer(bufsize)

    ! Initialise MPI

    call MPI_INIT(ierr)

    comm = MPI_COMM_WORLD

    ! Get the size of the pool and rank of the proc

    call MPI_COMM_SIZE(MPI_COMM_WORLD, pool_size, ierr)

    ! Get the rank for each processor

    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

    ! Create a buffer space: first get size of integer

    call MPI_TYPE_EXTENT(MPI_INTEGER, sizeof_int, ierr)

    ! Now create buffer space

    call MPI_BUFFER_ATTACH(buffer, bufsize, ierr)

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

            call MPI_Bsend(to_send, 1, MPI_INTEGER, &
                right, 0, comm, ierr)

            call MPI_Recv(recvd_data, 1, MPI_INTEGER, &
                left, 0, comm, stat, ierr)

            ! Since we've used non-blocking receive, we need to implement a wait

!            call MPI_WAIT(request, send_stat, ierr)

            local_count = local_count + recvd_data

            to_send = recvd_data

    end do

    call MPI_BARRIER(comm, ierr)

    time2 = mpi_wtime()

    if (rank .eq. 0) then

        print *, "Time required: ", time2-time1

    end if

    ! Detach the buffer

!    call MPI_BUFFER_DETACH(1000, ierr)

    ! Finalise MPI

    call MPI_FINALIZE(ierr)

end program main
