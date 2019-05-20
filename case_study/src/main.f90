!
! Jack Tyler: Image re/de-composition code; serial version.
!
! Modified 25/4 to compute the average value of the numbers every 10% of iterations.
! Modified 25/4 to cease computation whenever the difference slips below 0.1.
!
! Modified 20/5 to modularise source code.
!

program main

#ifdef PARALLEL2d

    use wrappers_2d

#elif SERIAL 

    use mpi,    only :: mpi_wtime
    use serial

#endif

    implicit none

    !
    ! PGM data arrays
    !

    double precision, allocatable       :: old(:,:)         ! "old" array: for iteration
    double precision, allocatable       :: buf(:,:)         ! Buffer_array: initial read-in location
    double precision, allocatable       :: new(:,:)         ! "new" array: for iteration
    double precision, allocatable       :: edge(:,:)        ! Array containing edge values
    double precision, allocatable       :: masterbuf(:,:)   ! Data array for rank 0

    !
    ! Algorithm variables
    !

    double precision, parameter         :: frac = .25d0     ! Optimisation: compute fraction

    double precision                    :: delta_global = 10! Delta flag: set unnecessarily high
    double precision                    :: local_delta = 100

    integer                             :: i                ! Loop variables
    integer                             :: j                ! Loop variables
    integer                             :: num_iters

    !
    ! Timing variables
    !

    double precision                    :: time_start
    double precision                    :: time_finish

    !
    ! Execute program
    !

    !
    ! In order to initialise our MPI runtime, topology and custom datatypes, we first need
    ! to allocate nbrs and dims to be assigned in the function.
    !
    ! If running serially, these arrays are not needed, hence the #ifdef.
    !

#ifdef PARALLEL2d

    allocate(nbrs(num_dims*2))
    allocate(dims(num_dims))

#endif

    !
    ! In parallel runtimes: initialise MPI, cartesian topology (2D), get the rank and size of the new pool.
    !                       Get arrangement of processors in the new topology(dims) and their communicator (cart_comm)
    !                       Print welcome message.
    !
    ! In serial runtimes:   Print welcome message. 
    !                       Set rank = 0.
    !                       Set all other variables = -1.
    !

    call initialise(pool_size, rank, nbrs, dims, cart_comm)

    !
    ! Initialise data arrays for the algorithm, based on the size of the problem.
    ! Serially, Mp = M, Np = N.
    !
    ! Mp, Np, M, N provided by the problem_constants module.
    !

    allocate(masterbuf(M, N))
    allocate(old(0:Mp+1, 0:Np+1))
    allocate(new(0:Mp+1, 0:Np+1))
    allocate(edge(Mp, Np))
    allocate(buf(Mp, Np)) ! Potentially not needed

    !
    ! Read the input file given by fname into masterbuf: in parallel runs, done only on rank 0.
    !
    ! fname provided by the problem_constants module.
    !

    if (rank .eq. 0) then

        call pgmread(fname, masterbuf)

    end if

    !
    ! In parallel runtimes, we now need to send the data to all the other processors:
    ! go from masterbuf into edge.
    !
    ! We do this using SCATTERV for the parallel case using custom vector and subarray types.
    ! For serial, we just set edge = masterbuf.
    !

    call send_data(masterbuf, Mp*Np, edge, Mp * Np, cart_comm)

    ! Set all entries in the old array to 255, as part of the algorithm.

    old(:,:) = 255

    call util_print_onrank0("Iterating...", rank)

    ! Initialise the iteration counter

    num_iters = 0

    ! We use a stopping criterion for the algorithm below:
    ! we want to proceed while the criterion is less than our threshold, and
    ! we haven't hit the maximum allowed number of iterations.
    !
    ! In order to properly time the execution, we'll also use an MPI_BARRIER
    ! and mpi_wtime() here. mpi_wtime() is used even if we're running
    ! serially.

#ifdef PARALLEL2d

    call MPI_BARRIER(cart_comm, ierr)

#endif

    time_start = mpi_wtime()

    do while (delta_global > stopping_criterion .and. num_iters < max_iters)

#ifdef PARALLEL2d        

        !
        ! We must first transfer our 'halo' regions - our overlaps - when running
        ! in parallel.
        !

        call mpi_send_halos(old, Np, M, dims, nbrs, cart_comm)

#endif

        !
        ! Perform the main computation loop: update new based on old and edge.
        !

        do j = 1, Np

            do i = 1, Mp ! Column major

                new(i, j) = frac * ( old(i-1, j) + &
                    old(i+1, j) + old(i, j+1) + &
                    old(i, j-1) - edge(i, j) )

            end do

        end do

        !
        ! If we're at an iteration specified to be one of our checking iterations
        ! (we compute the stopping criterion here only, to save computational time),
        ! then compute the criterion, and print a progress message.
        !

        if ((mod(num_iters, check_int) .eq. 0)) then

            call check_criterion(new, old, cart_comm, num_dims, num_iters, average, delta_global)

        end if

        !
        ! Now copy new into old for the next loop.
        ! util_array_copy performs additional error checking.
        !

        call util_array_copy(old(1:Mp, 1:Np), new(1:Mp, 1:Np))

        ! Update the iteration counter.

        num_iters = num_iters + 1

    end do

    !
    ! Now that this loop has completed, we want to re-time.
    !

#ifdef PARALLEL2d

    call MPI_BARRIER(cart_comm, ierr)

#endif

    time_finish = mpi_wtime()

    ! Call a subroutine that prints the appropriate status message for the completion of the loop.
    ! If num_iters = max_iters, an error is thrown. Else, a print of the time is given.

    call util_printfinish(num_iters, time_start, time_finish, rank)

    ! Now we need to gather the data from all processors if running in parallel, or
    ! copy our data from old to masterbuf if running serially.

    call gather_data(num_dims, old, Mp * Np, masterbuf, Mp * Np, cart_comm)

    ! Lastly, write masterbuf back to the image file. We only want rank 0 to do this,
    ! so enclose in an if.

    if (rank .eq. 0) then

        call pgmwrite('write_192x128.pgm', masterbuf)

    end if

    ! Deallocate variables used in all cases of the program.

    deallocate(new)
    deallocate(edge)
    deallocate(old)
    deallocate(buf)

    ! Deallocate variables used only in the parallel cases.

#ifdef PARALLEL2d

    deallocate(dims)
    deallocate(nbrs)

#endif

    call finalise()

end program main

