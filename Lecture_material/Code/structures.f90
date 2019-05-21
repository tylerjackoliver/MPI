program main

    use mpi

    implicit none

    type compound

        integer             :: ival
        double precision    :: dval

    end type compound

    ! Number of variables in type
    integer, parameter      :: cnt = 2

    integer                 :: rank
    integer                 :: pool_size
    integer                 :: comm
    integer                 :: new_comm
    integer                 :: ierr

    type(compound)          :: x

    !
    ! Address stuff
    !

    ! Arrays of displacements and lengths
    ! Allocated to size cnt

    integer, allocatable    :: array_of_blocklengths(:)
    integer(kind= &
        MPI_ADDRESS_KIND),&
        allocatable         :: array_of_displacements(:)

    integer, allocatable    :: array_of_types(:)

    ! Address variable

    integer(kind= & 
        MPI_ADDRESS_KIND)   :: address
    integer(kind= &
        MPI_ADDRESS_KIND)   :: address_offset

    integer                 :: newtype ! New data type!

    ! Prepare for the declaration of a structure

    allocate(array_of_blocklengths(cnt))
    allocate(array_of_types(cnt))
    allocate(array_of_displacements(cnt)) ! Equal to 2 here, but not always!

    ! Initialise MPI

    call MPI_INIT(ierr)

    comm = MPI_COMM_WORLD

    ! Construct the blocklength arrays

    array_of_blocklengths(0) = 1
    array_of_blocklengths(1) = 1

    ! Construct the array of types

    array_of_types(1) = MPI_INTEGER
    array_of_types(2) = MPI_DOUBLE

    !
    ! Now we need to construct the displacements array
    ! This is the offset between each variable in the structure
    ! We have cnt elements in the structure, so we bias
    ! to element 1.
    !
    ! Because MPI is a library and not a compiler,
    ! we need to tell MPI what the compiler has done, essentially
    !

    call MPI_GET_ADDRESS(x%ival, address_offset, ierr)
    call MPI_GET_ADDRESS(x%dval, address, ierr)

    array_of_displacements(1) = 0
    array_of_displacements(2) = address - address_offset

    ! Create the structure now

    call MPI_TYPE_CREATE_STRUCT(cnt, array_of_blocklengths, &
        array_of_displacements, array_of_types, newtype, ierr)
    
    ! Commit the datatype so we can use it!

    call MPI_TYPE_COMMIT(newtype, ierr)

    call MPI_FINALIZE(ierr)

end program main
