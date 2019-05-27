!
! A simple solution to the Case Study exercise from the EPCC MPI course.
!
! Communications is done using the sendrecv routine; a proper solution
! would use non-blocking communications (ie some combination of issend/recv
! and ssend/irecv).
!
! Note that the special rank of MPI_PROC_NULL is a "black hole" for
! communications (similar to /dev/null in Unix). The use of this value
! for processes off the edges of the image means we do not need any
! additional logic to ensure that processes at the edges do not attempt
! to send to or receive from invalid ranks (ie rank = -1 and rank = P).
!

program casestudy

  use mpi
  use pgmio
 
  implicit none

  integer, parameter :: M=768, N=768

  integer, parameter :: P = 2

  integer, parameter :: MAXITER   = 1500
  integer, parameter :: PRINTFREQ =  200

  integer, parameter :: MP=M, NP=N/P

  real, dimension(0:MP+1, 0:NP+1) :: new, old, edge

  real, dimension(MP,NP) :: buf
  real, dimension(M, N ) :: masterbuf

  integer, parameter :: maxlen = 32

  character*(maxlen) :: filename

  integer :: i, j, iter

  integer :: comm, rank, size, ierr, uprank, dnrank

  integer, dimension(MPI_STATUS_SIZE) :: status

  call mpi_init(ierr)

  comm = MPI_COMM_WORLD

  call mpi_comm_size(comm, size, ierr)
  call mpi_comm_rank(comm, rank, ierr)

  if (size .ne. P) then
    if (rank .eq. 0) write(*,*) 'ERROR: size = ', size, ', P = ', P
    call mpi_finalize(ierr)
    stop
  end if

  if (rank .eq. 0) then

    write(*,*) 'Processing ', M, ' x ' , N, ' image on ', P, ' processes'
    write(*,*) 'Number of iterations = ', MAXITER

    filename = 'edge768x768.pgm'

    write(*,*)
    write(*,*) 'Reading ', filename
    write(*,*)

    call pgmread(filename, masterbuf)

  end if

  call mpi_scatter(masterbuf, MP*NP, MPI_REAL, &
                   buf,       MP*NP, MPI_REAL, &
                   0, comm, ierr)

  do j = 1, NP
    do i = 1, MP

      edge(i,j) = buf(i,j)

    end do
  end do


  do j = 0, NP+1
    do i = 0, MP+1

      old(i,j) = 255.0

    end do
  end do

  uprank = rank + 1
  if (uprank .ge. P) uprank = MPI_PROC_NULL

  dnrank = rank - 1
  if (dnrank .lt. 0) dnrank = MPI_PROC_NULL

  do iter = 1, MAXITER
    
    if (mod(iter,PRINTFREQ) .eq. 0)  then

      if (rank .eq. 0) write(*,*) 'Iteration ', iter

    end if

    call mpi_sendrecv(old(1,NP), MP, MPI_REAL, uprank, 1, &
                      old(1,0),  MP, MPI_REAL, dnrank, 1, &
                      comm, status, ierr)
    call mpi_sendrecv(old(1,1),    MP, MPI_REAL, dnrank, 1, &
                      old(1,NP+1), MP, MPI_REAL, uprank, 1, &
                      comm, status, ierr)

    do j = 1, NP
      do i = 1, MP

        new(i,j) = 0.25*(old(i+1,j)+old(i-1,j)+old(i,j+1)+old(i,j-1) &
                         - edge(i,j))

      end do
    end do

    do j = 1, NP
      do i = 1, MP

        old(i,j) = new(i,j)

      end do
    end do

  end do

  if (rank .eq. 0) then
    write(*,*)
    write(*,*) 'Finished ', iter-1, ' iterations'
  end if

!  Gather data

  do j = 1, NP
    do i = 1, MP

      buf(i,j) = old(i,j)

    end do
  end do

  call mpi_gather(buf,       MP*NP, MPI_REAL, &
                  masterbuf, MP*NP, MPI_REAL, &
                  0, comm, ierr)

  if (rank .eq. 0) then

    filename='image768x768.pgm'
    write(*,*)
    write(*,*) 'Writing ', filename
    call pgmwrite(filename, masterbuf)

  end if

  call mpi_finalize(ierr)

end program casestudy
