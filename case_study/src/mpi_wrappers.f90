module mpi_wrappers

    use mpi 
    use neighbour_indexes

    implicit none

    contains

    subroutine mpi_initialise_standard_topology(num_dims, pool_size, old_comm, new_comm, new_rank, nbrs)
        
        use neighbour_indexes                   ! Provides standardised values of left, right

        integer, intent(in)                     :: num_dims
        integer, intent(in)                     :: pool_size
        integer, intent(in)                     :: old_comm

        integer, intent(out)                    :: new_comm
        integer, intent(out)                    :: new_rank
        integer, intent(inout), dimension(:)    :: nbrs

        integer                                 :: dims(1)
        integer                                 :: ierr
        integer                                 :: x_dir 
        integer                                 :: displacement
        integer                                 :: i

        logical                                 :: periodic(1)
        logical                                 :: reorder(1)

        reorder = .false.

        if (num_dims .eq. 1) then

            call initialise_standard_topology_1d(pool_size, old_comm, new_comm, new_rank, nbrs)
        
        elseif (num_dims .eq. 2) then

            call initialise_standard_topology_2d(pool_size, old_comm, new_comm, new_rank, nbrs)

        end if

    end subroutine mpi_initialise_standard_topology


    subroutine initialise_standard_topology_1d(pool_size, old_comm, new_comm, new_rank, nbrs)
        
        use neighbour_indexes                   ! Provides standardised values of left, right

        integer, intent(in)                     :: pool_size
        integer, intent(in)                     :: old_comm

        integer, intent(out)                    :: new_comm
        integer, intent(out)                    :: new_rank
        integer, intent(inout), dimension(:)    :: nbrs

        integer, parameter                      :: num_dims = 1
        integer                                 :: dims(1)
        integer                                 :: ierr
        integer                                 :: x_dir 
        integer                                 :: displacement
        integer                                 :: i

        logical                                 :: periodic(1)
        logical                                 :: reorder(1)

        reorder = .false.

        ! Set dims to zero and periodic to false

        do i = 1, num_dims

            periodic(i) = .false.
            dims(i)     = 0

        end do

        x_dir = 0
        displacement = 1
        
        ! Now create the topology

        call MPI_DIMS_CREATE(pool_size, num_dims, dims, ierr)
        call MPI_CART_CREATE(old_comm, num_dims, dims, periodic, reorder, new_comm, ierr)
        call MPI_COMM_RANK(new_comm, new_rank, ierr)
        call MPI_CART_SHIFT(new_comm, x_dir, displacement, nbrs(left), nbrs(right), ierr)

    end subroutine initialise_standard_topology_1d


    subroutine initialise_standard_topology_2d(pool_size, old_comm, new_comm, new_rank, nbrs)
        
        use neighbour_indexes                   ! Provides standardised values of left, right

        integer, intent(in)                     :: pool_size
        integer, intent(in)                     :: old_comm

        integer, intent(out)                    :: new_comm
        integer, intent(out)                    :: new_rank
        integer, intent(inout), dimension(:)    :: nbrs

        integer, parameter                      :: num_dims = 1
        integer                                 :: dims(1)
        integer                                 :: ierr
        integer                                 :: x_dir 
        integer                                 :: displacement
        integer                                 :: i

        logical                                 :: periodic(1)
        logical                                 :: reorder(1)

        reorder = .false.

        print *, "Pass."

        new_comm = 1
        new_rank = 1

    end subroutine initialise_standard_topology_2d


    subroutine mpi_send_data(num_dims, to_send, size_to_send, to_recv, size_to_recv, comm)
        
        integer,                          intent(in)    :: num_dims 

        double precision, dimension(:,:), intent(in)    :: to_send
        
        integer,                          intent(in)    :: size_to_send
        
        double precision, dimension(:,:), intent(inout) :: to_recv

        integer,                          intent(in)    :: size_to_recv
        integer,                          intent(in)    :: comm

        integer                                         :: ierr

        if (num_dims .eq. 1) then

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

    end subroutine mpi_send_data_1d


    subroutine mpi_send_data_2d(to_send, size_to_send, to_recv, size_to_recv, comm)

        double precision, dimension(:,:), intent(in)    :: to_send
        
        integer,                          intent(in)    :: size_to_send
        
        double precision, dimension(:,:), intent(inout) :: to_recv

        integer,                          intent(in)    :: size_to_recv
        integer,                          intent(in)    :: comm

        integer                                         :: ierr

        print *, "Pass"

    end subroutine mpi_send_data_2d

    
    subroutine mpi_send_halos(num_dim, old, Np, M, nbrs, cart_comm)

        integer,                            intent(in)      :: num_dim

        double precision, dimension(0:,0:), intent(inout)   :: old

        integer,                            intent(in)      :: Np
        integer,                            intent(in)      :: M
        integer,          dimension(:),     intent(in)      :: nbrs
        integer,                            intent(in)      :: cart_comm

        if (num_dim .eq. 1) then

            call send_halos_1d(old, Np, M, nbrs, cart_comm)

        elseif (num_dim .eq. 2) then

            call send_halos_2d(old, Np, M, nbrs, cart_comm)

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

        call MPI_SENDRECV(old(1, Np), M, MPI_DOUBLE, nbrs(right), 0, &
        old(1,0), M, MPI_DOUBLE, nbrs(left), 0, cart_comm, &
        recv_status, ierr)

        ! Send/recv left

        call MPI_SENDRECV(old(1, 1), M, MPI_DOUBLE, nbrs(left), 0, &
        old(1, Np+1), M, MPI_DOUBLE, nbrs(right), 0, cart_comm, &
        recv_status, ierr)

    end subroutine send_halos_1d


    subroutine send_halos_2d(old, Np, M, nbrs, cart_comm)

        use neighbour_indexes

        double precision, dimension(0:,0:), intent(inout)   :: old

        integer,                          intent(in)        :: Np
        integer,                          intent(in)        :: M
        integer,          dimension(:),   intent(in)        :: nbrs
        integer,                          intent(in)        :: cart_comm

        integer, dimension(mpi_status_size)                 :: recv_status      ! Receive status, non-blocking send

        integer                                             :: ierr

        print *, "Pass."

    end subroutine send_halos_2d


    subroutine gather_data(to_send, send_size, to_recv, recv_size, comm)

        double precision, dimension(:,:), intent(in)        :: to_send

        integer,                          intent(in)        :: send_size
        
        double precision, dimension(:,:), intent(out)       :: to_recv

        integer,                          intent(in)        :: recv_size
        integer,                          intent(in)        :: comm

        integer, dimension(mpi_status_size)                 :: recv_status      ! Receive status, non-blocking send

        integer                                             :: ierr

        call MPI_GATHER(to_send, send_size, MPI_DOUBLE_PRECISION, to_recv, recv_size,&
                        MPI_DOUBLE_PRECISION, 0, comm, ierr)

    end subroutine gather_data


    ! subroutine get_array_element_1d(new, old, edge, i, j)
        
    !     double precision, dimension(:,:), intent(inout) :: new

    !     double precision, dimension(:,:), intent(in)    :: old
    !     double precision, dimension(:,:), intent(in)    :: edge

    !     integer,                          intent(in)    :: i
    !     integer,                          intent(in)    :: j

    !     double precision                                :: frac = .25d0
        
    !     double precision                                :: temp_sum
        
    !     temp_sum = old(i-1, j) + old(i+1, j) + old(i, j+1) +  old(i, j-1)

    !     new(i, j) = frac * ( temp_sum - edge(i, j) )

    ! end subroutine get_array_element_1d


    subroutine get_global_delta(delta, comm, global_delta)

    
        double precision,                 intent(in)    :: delta
    
        integer,                          intent(in)    :: comm 
    
        double precision,                 intent(out)   :: global_delta
    
        integer                                         :: ierr
    
        call MPI_ALLREDUCE(delta, global_delta, 1, MPI_DOUBLE_PRECISION, &
                           MPI_MAX, comm, ierr)
    
    end subroutine get_global_delta


    subroutine get_average(array, M, N, comm, average)

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

    end subroutine get_average
        

end module mpi_wrappers

    