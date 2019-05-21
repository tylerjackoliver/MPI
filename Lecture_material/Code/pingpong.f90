    
      program main

          use mpi

          implicit none

          integer                               :: ierr
          integer                               :: comm, rank, comm_size
          integer                               :: i
          integer, dimension(MPI_STATUS_SIZE)   :: stat
          integer                               :: extent
          integer                               :: length

          double precision, allocatable         :: data_arr(:)
          integer                               :: numiter

          double precision                      :: tstart, time, tstop, totmess

          character(len=10)                     :: temp_char10


          call mpi_init(ierr)

          ! Portability: specify the communicator now

          comm = MPI_COMM_WORLD

          ! Get rank

          call MPI_COMM_RANK(comm, rank, ierr)

          ! Get command line argument for length of array

          if (rank .eq. 0) then

              call getarg(1, temp_char10)
              read(temp_char10, *) length

              call getarg(2, temp_char10)
              read(temp_char10, *) numiter

          end if

          ! Broadcast everything to the other procs

          call MPI_BCAST(length, 1, MPI_INTEGER, 0, comm, ierr)
          call MPI_BCAST(numiter, 1, MPI_INTEGER, 0, comm, ierr)

          ! Now check the other procs out of the game

          if (rank .gt. 1) then

              write(*, '(A, I2, A)'), "I am rank", rank, &
                                      "and I'm checking out, yo."

          end if

          ! Allocate the data arrays

          allocate(data_arr(length))

          do i = 1,length

            data_arr(i) = rank + 10.d0

          end do

          ! Barrier and start timing

          call MPI_BARRIER(comm, ierr)

          tstart = MPI_WTIME()

          do i = 1, numiter

            if (rank .eq. 0) then

                call MPI_SSEND(data_arr(1), length, &
                    MPI_DOUBLE_PRECISION, 1, 0, comm, ierr)
                call MPI_RECV(data_arr(1), length, &
                    MPI_DOUBLE_PRECISION, 1, 0, comm, stat, ierr)

            elseif (rank .eq. 1) then

                call MPI_RECV(data_arr(1), length, &
                    MPI_DOUBLE_PRECISION, 0, 0, comm, stat, ierr)
                call MPI_SSEND(data_arr(1), length, &
                    MPI_DOUBLE_PRECISION, 0, 0, comm, ierr)

            end if

          end do

          tstop = MPI_WTIME()

          time = tstop - tstart

          call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, extent, ierr)

          if (rank .eq. 0) then

              totmess = 2.d0 * extent * length / (1024.d0) * numiter &
                  /1024.d0 ! Convert bytes -> MB

              write (*,*) "Ping-pong of twice", extent*length, &
                  "bytes, for", numiter, "times."
              write(*,*) "Total compute time is", time, "s."
              write(*,*) "Latency (time per msg) is", &
                  time/numiter * .5d0, "s."
              write(*,*) "Bandwidth (msg. per time)", totmess/time, &
                  "MB/s."

              if (time .lt. 1.d0) then

                  write(*,*) "Watch it, the timing is < 1s!"

              end if

          end if

          deallocate(data_arr)

          call mpi_finalize(ierr)

      end program main
