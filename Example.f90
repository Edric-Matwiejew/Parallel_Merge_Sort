!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!132

module Rand

    use, intrinsic :: iso_fortran_env, sp => real32, dp => real64

    implicit none

    contains

    subroutine  Seed_Random_Number(seed)

        ! Seeds the intrinsic fotran function 'random_number'. To to ensure that each local array segment for this example contains
        ! a unique set of values.

        integer, intent(in) :: seed

        integer :: seed_length
        integer, dimension(:), allocatable :: seed_array

        call random_seed(size=seed_length)
        allocate(seed_array(seed_length))
        seed_array = seed
        call random_seed(put=seed_array)
        deallocate(seed_array)

    end subroutine Seed_Random_Number



end module Rand

program Example

    use :: Rand
    use :: MPI
    use :: Parallel_Merge, only: Parallel_Merge_Sort
    use :: Parallel_Merge_Indexed, only: Parallel_Merge_Sort_Indexed

    implicit none

    real(dp), dimension(:), allocatable :: local_array, collated_array
    integer, dimension(:), allocatable :: local_indices, collated_indices
    integer, dimension(:), allocatable :: block_lens, disps
    integer :: local_size = 10 ! each node will have an array segment of size 10
    integer :: total_size 

    integer :: i

    ! MPI Environment
    integer :: flock ! the number of MPI nodes
    integer :: rank
    integer :: root = 0
    integer :: ierr

    ! Initialized MPI environment

    call MPI_init(ierr)
   
    ! Determine total flock size and  node rank
    call MPI_comm_size(MPI_comm_world, flock, ierr)
    call MPI_comm_rank(MPI_comm_world, rank, ierr)

    total_size = local_size*flock 

    allocate(local_array(local_size), local_indices(local_size))
    allocate(block_lens(flock), disps(flock))

    if (rank == root) then
        ! These arrays  will contain the final sorted array and array indicies. 
        allocate(collated_indices(total_size), collated_array(total_size))
    endif    

    ! In this example each node has the same block length, this is not required. block_len may vary per node.
    block_lens = 10

    do i = 1, flock
        disps(i) = (i - 1)*local_size
    enddo

    call Seed_Random_Number(rank*local_size + i)

    ! Populating local array segements
    do i = 1, local_size
        call random_number(local_array(i))
    enddo

    call Parallel_Merge_Sort(local_array, collated_array, root, block_lens, disps, MPI_comm_world)

    if (rank == root) then 

        write(*,*) "MPI_MERGE_SORT"
        write(*,*)

        do i = 1, size(collated_array)
            write(*,*) collated_array(i)
        enddo

        write(*,*)

    endif
    
    call mpi_barrier(MPI_comm_world, ierr)

    ! Populate_local array segment and local_indices
    do i = 1, local_size
        call random_number(local_array(i))
        local_indices(i) = rank*local_size + i
    enddo

    call Parallel_Merge_Sort_Indexed(local_array, local_indices, collated_array, collated_indices, root, block_lens, disps, & 
            & MPI_comm_world)

    if (rank == root) then 

        write(*,*) "MPI_MERGE_SORT_INDEXED"
        write(*,*)

        do i = 1, size(collated_array)
            write(*,*) collated_indices(i), collated_array(i)
        enddo

    endif

    call MPI_finalize(ierr)     

end program Example
