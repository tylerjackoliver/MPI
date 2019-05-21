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
  integer            :: local_set(2)                ! Number of iterations to perform locally
  integer            :: i                           ! Loop variable
  integer, &
      dimension(MPI_STATUS_SIZE) :: stat

  double precision   :: local_value                 ! Local value of pi
  double precision   :: pi_approx                   ! Approximation to Pi: master thread only
  double precision   :: comm_size_quotient          ! Slight optimisation
  double precision   :: received_pi

  ! Initialise MPI

  call MPI_Init(ierr)

  local_value = 0

  ! Get each PROC's rank and the total size of the pool

  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, comm_size, ierr)

  ! Based on its rank, divide up the work

  comm_size_quotient = ceiling(dble(N / comm_size))

  local_set(1) = (rank * comm_size_quotient) + 1
  local_set(2) = min(dble((rank+1) * comm_size_quotient), dble(N))


  do i = local_set(1), local_set(2)

    local_value = local_value + 1.d0/(1.d0 + ( (i-.5d0)/N) ** 2)

  end do

  ! Now receive everything: must quicker to do MPI_REDUCE but oh well

  if (rank .eq. 0) then

      pi_approx = local_value

      do i = 1, comm_size - 1

        call MPI_RECV(received_pi, 1, MPI_DOUBLE_PRECISION, i, 0, MPI_COMM_WORLD, stat, ierr)
        pi_approx = pi_approx + received_pi

      end do

      pi_approx = 4.d0 * pi_approx / N

      ! Print status message
      
      write(*,*) "I am the master of all: I have pi to be", pi_approx

  else

      call MPI_Ssend(local_value, 1, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, ierr)

  end if

  ! Finalise MPI

  call MPI_Finalize(ierr)

end program hello
