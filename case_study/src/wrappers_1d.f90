module wrappers_1d

    use problem_constants
    use mpi
    use utility_functions
    use neighbour_indexes
    use pgmio                                               ! Use .pgm reading library

    implicit none

    integer :: Mp = M
    integer :: Np = ceiling(dble(N/P))                                ! Assumes N/P is perfect 

    PUBLIC :: Mp, Np

    contains


        subroutine mpi_send_data(num_dims, to_send, size_to_send, to_recv, size_to_recv, comm)
                
            integer,                          intent(in)    :: num_dims 

            double precision, dimension(:,:), intent(in)    :: to_send
            
            integer,                          intent(in)    :: size_to_send
            
            double precision, dimension(:,:), intent(inout) :: to_recv

            integer,                          intent(in)    :: size_to_recv
            integer,                          intent(in)    :: comm

            integer                                         :: ierr

            if (num_dims .eq. 0) then

                to_recv = to_send

            else if (num_dims .eq. 1) then

                call mpi_send_data_1d(to_send, size_to_send, to_recv, size_to_recv, comm)

            elseif (num_dims .eq. 2) then

                call mpi_send_data_2d(to_send, size_to_send, to_recv, size_to_recv, comm)

            end if

        end subroutine


        subroutine mpi_send_data_1d(to_send, size_to_send, to_recv, size_to_recv, comm)

            double precision, dimension(:,:), intent(in)    :: to_send
            
            integer,                          intent(in)    :: size_to_send
            
            double precision, dimension(:,:), intent(inout) :: to_recv

            integer,                          intent(in)    :: size_to_recv
            integer,                          intent(in)    :: comm

            integer                                         :: ierr

            call MPI_SCATTER(to_send, size_to_send, MPI_DOUBLE_PRECISION, to_recv, &
                            size_to_recv, MPI_DOUBLE_PRECISION, 0, comm, ierr)

            ! call mpi_send_data(num_dims, masterbuf, Mp*Np, buf, Mp * Np, cart_comm)

        end subroutine mpi_send_data_1d


        subroutine mpi_send_data_2d(to_send, size_to_send, to_recv, size_to_recv, comm)

            double precision, dimension(:,:), intent(in)    :: to_send
            
            integer,                          intent(in)    :: size_to_send
            
            double precision, dimension(:,:), intent(inout) :: to_recv

            integer,                          intent(in)    :: size_to_recv
            integer,                          intent(in)    :: comm

            integer                                         :: ierr

            ! call MPI_Scatterv(to_send, counts, displacements, master_type, &
            !                   to_recv, Mp*Np, MPI_DOUBLE_PRECISION, 0, comm, ierr)

        end subroutine mpi_send_data_2d


        subroutine get_comm_size(comm, size)

            integer, intent(in) :: comm
            integer, intent(out):: size 

            integer             :: ierr

            call MPI_COMM_SIZE(comm, size, ierr)

        end subroutine get_comm_size


        subroutine mpi_initialise(num_dims, pool_size, rank, nbrs, dims, cart_comm)

            integer, intent(in)     :: num_dims

            integer, intent(out)    :: pool_size
            integer, intent(out)    :: rank

            integer, dimension(:), intent(inout)  :: nbrs
            integer, dimension(:), intent(inout)  :: dims

            integer,               intent(out)    :: cart_comm

            integer                 :: i
            integer                 :: ierr
            integer                 :: comm

            call MPI_INIT(ierr)

            comm = MPI_COMM_WORLD

            call get_comm_size(comm, pool_size)

            ! If size =/ P, then exit

            if (pool_size .ne. P) then

                call MPI_FINALIZE(ierr)
                error stop "Number of processors does not match the expected number of processes."

            end if
            
            ! Define local problem array size

            call util_print_welcomemessage(rank)

            Mp = M
            Np = ceiling(dble(N/P))                                ! Assumes N/P is perfect
            
            call initialise_standard_topology_1d(num_dims, dims, cart_comm, nbrs, rank)

        end subroutine

        subroutine initialise_standard_topology_1d(num_dims, dims, cart_comm, nbrs, rank)
            ! initialise_new_standard_topology_1d(num_dims, dims, cart_comm, nbrs, rank)

            use neighbour_indexes

            integer,               intent(in)       :: num_dims
            integer, dimension(:), intent(inout)    :: dims
            integer,               intent(out)      :: cart_comm
            integer, dimension(:), intent(inout)    :: nbrs

            integer,               intent(out)      :: rank

            logical                                 :: periodic(1)
            logical                                 :: reorder

            integer                                 :: ierr
            integer                                 :: x_dir
            integer                                 :: displacement
            integer                                 :: i

            integer                                 :: pool_size

            reorder      = .false.
            x_dir        = 0
            displacement = 1

            ! Set dims to zero and periodic to false

            do i = 1, num_dims

                periodic(i) = .false.
                dims(i)     = 0

            end do

            call get_comm_size(MPI_COMM_WORLD, pool_size)

            ! Now create the topology

            call MPI_DIMS_CREATE(pool_size, num_dims, dims, ierr)
            call MPI_CART_CREATE(MPI_COMM_WORLD, num_dims, dims, periodic, reorder, cart_comm, ierr)
            call MPI_COMM_RANK(cart_comm, rank, ierr)
            call MPI_CART_SHIFT(cart_comm, x_dir, displacement, nbrs(left), nbrs(right), ierr)

        end subroutine initialise_standard_topology_1d


        ! subroutine mpi_initialise_standard_topology(num_dims, pool_size, old_comm, new_comm, new_rank, dims, nbrs)
                
        !     integer, intent(in)                     :: num_dims
        !     integer, intent(in)                     :: pool_size
        !     integer, intent(in)                     :: old_comm

        !     integer, intent(out)                    :: new_comm
        !     integer, intent(out)                    :: new_rank
        !     integer, intent(inout), dimension(:)    :: dims

        !     integer, intent(inout), dimension(:)    :: nbrs

        !     integer                                 :: ierr

        !     if (num_dims .eq. 0) then   ! running serially

        !         new_rank = 0            ! Fake that we're on rank 0
        !         new_comm = -1           ! Makes sure we get an error if we try to MPI
        !         nbrs     = 0            ! Makes sure we get an error if we try to MPI

        !     elseif (num_dims .eq. 1) then

        !         call initialise_standard_topology_1d(pool_size, old_comm, new_comm, new_rank, dims, nbrs)
            
        !     elseif (num_dims .eq. 2) then

        !         ! call initialise_standard_topology_2d(pool_size, old_comm, new_comm, new_rank, dims, nbrs)

        !     end if

        ! end subroutine mpi_initialise_standard_topology


        subroutine mpi_send_halos(num_dim, old, Np, M, dims, nbrs, cart_comm)

            integer,                            intent(in)      :: num_dim

            double precision, dimension(0:,0:), intent(inout)   :: old

            integer,                            intent(in)      :: Np
            integer,                            intent(in)      :: M
            integer,          dimension(:),     intent(in)      :: dims
            integer,          dimension(:),     intent(in)      :: nbrs
            integer,                            intent(in)      :: cart_comm

            if (num_dim .eq. 0) then

                ! Do nothing - old is already what it should be

            elseif (num_dim .eq. 1) then

                call send_halos_1d(old, Np, M, nbrs, cart_comm)

            elseif (num_dim .eq. 2) then

                ! call send_halos_2d(old, Np, M, dims, nbrs, cart_comm)

            end if

        end subroutine mpi_send_halos


        subroutine send_halos_1d(old, Np, M, nbrs, cart_comm)

            use neighbour_indexes

            double precision, dimension(0:,0:), intent(inout)   :: old

            integer,                          intent(in)        :: Np
            integer,                          intent(in)        :: M
            integer,          dimension(:),   intent(in)        :: nbrs
            integer,                          intent(in)        :: cart_comm

            integer, dimension(mpi_status_size)                 :: recv_status      ! Receive status, non-blocking send

            integer                                             :: ierr

            ! Send/recv right


            call MPI_SENDRECV(old(1, Np), M, MPI_DOUBLE, nbrs(up), 0, &
            old(1,0), M, MPI_DOUBLE, nbrs(down), 0, cart_comm, &
            recv_status, ierr)

            ! Send/recv left

            call MPI_SENDRECV(old(1, 1), M, MPI_DOUBLE, nbrs(down), 0, &
            old(1, Np+1), M, MPI_DOUBLE, nbrs(up), 0, cart_comm, &
            recv_status, ierr)

        end subroutine send_halos_1d

        subroutine mpi_get_global_delta(num_dims, delta, comm, global_delta)

            integer,                          intent(in)    :: num_dims

            double precision,                 intent(in)    :: delta

            integer,                          intent(in)    :: comm 

            double precision,                 intent(out)   :: global_delta

            integer                                         :: ierr

            if (num_dims .eq. 0) then

                global_delta = delta

            else

                call MPI_ALLREDUCE(delta, global_delta, 1, MPI_DOUBLE_PRECISION, &
                                MPI_MAX, comm, ierr)

            end if

        end subroutine mpi_get_global_delta


        subroutine mpi_get_average(array, M, N, comm, average)

            double precision, dimension(:,:), intent(in)    :: array 

            integer,                          intent(in)    :: comm
            integer,                          intent(in)    :: M
            integer,                          intent(in)    :: N

            double precision,                 intent(out)   :: average

            integer                                         :: ierr

            double precision                                :: local_sum
            
            local_sum = sum(array)

            call MPI_ALLREDUCE(local_sum, average, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)

            average = average / ( M * N )

        end subroutine mpi_get_average


        subroutine mpi_gather_data(num_dims, to_send, send_size, to_recv, recv_size, comm)

            integer,                          intent(in)        :: num_dims

            double precision, dimension(:,:), intent(in)        :: to_send

            integer,                          intent(in)        :: send_size
            
            double precision, dimension(:,:), intent(out)       :: to_recv

            integer,                          intent(in)        :: recv_size
            integer,                          intent(in)        :: comm

            integer, dimension(mpi_status_size)                 :: recv_status      ! Receive status, non-blocking send

            integer                                             :: ierr

            if (num_dims .eq. 0) then ! Serial case

                to_recv = to_send

            end if

            if (num_dims .eq. 1) then

                call MPI_GATHER(to_send, send_size, MPI_DOUBLE_PRECISION, to_recv, recv_size,&
                            MPI_DOUBLE_PRECISION, 0, comm, ierr)

            end if

        end subroutine mpi_gather_data


        subroutine finalize_runtime()

            integer                 :: ierr 

            call MPI_FINALIZE(ierr)

        end subroutine finalize_runtime

end module wrappers_1d


