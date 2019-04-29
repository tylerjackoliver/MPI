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
    
    use wrappers_2d

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
    ! Execute program
    !

    ! Define x and y dimensions

    ! Now allocate data arrays

    allocate(nbrs(num_dims*2))
    allocate(dims(num_dims))

    ! Initialise MPI, compute the size and the rank of the array

    call initialise(pool_size, rank, nbrs, dims, cart_comm)

    allocate(masterbuf(M, N))
    allocate(old(0:Mp+1, 0:Np+1))
    allocate(new(0:Mp+1, 0:Np+1))
    allocate(edge(Mp, Np))
    allocate(buf(Mp, Np))

    ! Step 1: read the edges data file into the buffer

    if (rank .eq. 0) then

        call pgmread(fname, masterbuf)

    end if

    ! Now we've read it on the master thread, time to scatter to all procs
    !

    call send_data(masterbuf, Mp*Np, edge, Mp * Np, cart_comm)

    ! Step 2: loop over M, N

    if (rank .eq. 0) then

        write(*,*) "Setting up arrays..."

    end if

    old(:,:) = 255

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

    call mpi_gather_data(num_dims, old, Mp * Np, masterbuf, Mp * Np, cart_comm)

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


end program main

