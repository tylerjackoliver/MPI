!
! Jack Tyler: Image re/de-composition code; serial version.
!


program main

    use mpi                                                 ! Use MPI library
    use pgmio                                               ! Use .pgm reading library

    !
    ! Variable declarations for PGM file I/O
    !

    integer                             :: M                ! Number of pixels horizontally
    integer                             :: N                ! Number of pixels vertically

    integer                             :: i                ! Loop variables
    integer                             :: j                ! Loop variables
    integer                             :: k                ! Loop variables

    integer, parameter                  :: num_iters = 1000 ! Number of iterations of jacobi

    double precision, allocatable       :: old_array(:,:)   ! "old" array: for iteration
    double precision, allocatable       :: buf(:,:)         ! Buffer_array: initial read-in location
    double precision, allocatable       :: new_array(:,:)   ! "new" array: for iteration
    double precision, allocatable       :: edge(:,:)        ! Array containing edge values

    double precision, parameter         :: frac = .25d0     ! Optimisation: compute fraction

    !
    ! Execute program
    !

    ! Define x and y dimensions

    M = 768
    N = 768

    ! Now allocate data arrays

    allocate(old_array(0:M+1, 0:N+1))
    allocate(new_array(0:M+1, 0:N+1))
    allocate(edge(0:M+1, 0:N+1))
    allocate(buf(M, N))

    ! Step 1: read the edges data file into the buffer

    call pgmread('edge768x768.pgm', buf)

    ! Step 2: loop over M, N

    write(*,*) "Setting up arrays..."

    do j = 1, N

        do i = 1, M

            edge(i, j) = buf(i, j)                          ! Read into buf
            old_array(i, j)  = 255                          ! Set old array to white (255)

        end do

    end do

    write(*,*) "Iterating..."

    ! Step 3: Iterate through our edges

    do k = 1, num_iters

        do j = 1, N

            do i = 1, M ! Column major

                new_array(i, j) = frac * ( old_array(i-1, j) + &
                    old_array(i+1, j) + old_array(i, j+1) + &
                    old_array(i, j-1) - edge(i, j) )

            end do

        end do

        old_array(1:M, 1:N) = new_array(1:M, 1:N)          ! Set old = new, excluding halos

    end do

    write(*,*) "Done! Copying back..."

    ! Step 4: copy old array back to buf

    buf(1:M, 1:N) = old_array(1:M, 1:N)

    ! Write buf to image

    call pgmwrite('write_192x128.pgm', buf)

    write(*,*) "Done"

    deallocate(new_array)
    deallocate(edge)
    deallocate(old_array)
    deallocate(buf)

end program main

