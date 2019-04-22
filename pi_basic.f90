!
!
! Fortran program for the parallel computation of pi.
!
!


program hello

  use mpi

  implicit none

  integer            :: ierr                        ! MPI error integer
  integer, parameter :: N = 840                     ! Order of approximation

  integer            :: rank, comm_size             ! Each proc's rank and the total size
  integer            :: local_set                   ! Number of iterations to perform locally
  integer            :: i                           ! Loop variable

  double precision   :: local_value                 ! Local value of pi
  double precision   :: pi_approx                   ! Approximation to Pi: master thread only

  ! Initialise MPI

  call MPI_Init(ierr)

  do i = 1, N

    pi_approx = pi_approx + 1.d0/(1.d0 + ( (i-.5d0)/N) ** 2)

  end do

  pi_approx = 4.d0 * pi_approx / N

  ! Get each PROC's rank

  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

  ! Print status message

  write(*,*) "I am process identifier", rank, "and I have pi to be", pi_approx

  ! Finalise MPI

  call MPI_Finalize(ierr)

end program hello
