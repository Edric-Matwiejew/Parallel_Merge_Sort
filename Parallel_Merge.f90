!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!132

module Parallel_Merge

    use, intrinsic :: iso_fortran_env, sp => real32, dp => real64
    use :: MPI

    implicit none

    private

    public :: Parallel_Merge_Sort

    contains

    subroutine Parallel_Merge_Sort(local_array, collated_array, root, block_lens, disps, MPI_communicator)

        ! MPI_Merge_Sort_Indexed takes a real valued array distributed across an mpi communicator. The array segements are sorted 
        ! locally using a merge sort and are gathered by the root process. The root process merges the sorted segements into the 
        ! final sorted array.

        real(dp), intent(inout), dimension(:) :: local_array ! array to be sorted
        real(dp), intent(inout), dimension(:) :: collated_array ! total sorted array to be formed at the root node 
        integer, intent(in) :: root ! node at which to form the final sorted array
        integer, intent(in), dimension(:) :: block_lens ! the size of the array at each MPI node, ordered by accending rank
        integer, intent(in), dimension(:) :: disps !displacement of each array section from the global array start (0 - indexed) 
        integer, intent(in) :: MPI_communicator

        ! MPI environment
        integer :: flock ! the number of MPI nodes
        integer :: rank
        integer :: ierr

        call MPI_comm_size(MPI_comm_world, flock, ierr)
        call MPI_comm_rank(MPI_comm_world, rank, ierr)

        call Merge_Sort(local_array, 1, size(local_array))

        call MPI_gatherv(local_array, size(local_array), MPI_double, collated_array, block_lens, disps, MPI_double, root, &
                & MPI_communicator, ierr)

        if (rank == root) then
            call Merge_Ordered_Lists(collated_array, flock)
        endif

    end subroutine Parallel_Merge_Sort


    subroutine Insertion_Sort(array)

        real(dp), intent(inout), dimension(:) :: array

        real(dp) :: temp
        integer :: i, j

        do i = 2, size(array)

            temp = array(i)
            j = i - 1

            do while (j >= 1)

                if (array(j) <= temp) exit

                    array(j + 1) = array(j)
                    j = j - 1

            enddo

            array(j + 1) = temp

        enddo

    end subroutine

    subroutine Merge(array, start, mid, finish)

        real(dp), intent(inout), dimension(:) :: array
        integer, intent(in) :: start, mid, finish

        real(dp), dimension(:), allocatable :: temp
        integer :: i, j, k

        allocate(temp(finish - start + 1))

        i = start
        j = mid + 1
        k = 1

        do while (i <= mid .and. j <= finish)

            if (array(i) <= array(j)) then

                temp(k) = array(i)

                k = k + 1
                i = i + 1

            else

                temp(k) = array(j)

                k = k + 1
                j = j+ 1

            endif

        enddo

        do while (i <= mid)

            temp(k) = array(i)

            k = k + 1
            i = i + 1

        enddo

        do while (j <= finish)

            temp(k) = array(j)

            k = k + 1
            j = j + 1

        enddo

        do i = start, finish
            array(i) = temp(i - start + 1)
        enddo

    end subroutine Merge

    recursive subroutine Merge_Sort(array, start, finish)

    ! Merge_Sort_Indexed performs a merge sort of the local input 'array'. Array segments of sizes less than 128 are sorted by an 
    ! insertion sort so as to prevent over-partitioning and limit openMP task creation.


        real(dp), intent(inout), dimension(:) :: array
        integer, intent(in) :: start, finish ! specifies the reigon of array to be sorted.

        integer :: mid

        if (start < finish) then
            if (finish - start >= 512) then

                mid = (start + finish) / 2

                !$omp taskgroup

                !$omp task shared(array) untied
                call Merge_Sort(array, start, mid)
                !$omp end task

                !$omp task shared(array) untied
                call Merge_Sort(array, mid + 1, finish)
                !$omp end task

                !$omp end taskgroup

                call Merge(array, start, mid, finish)

            else
                call Insertion_Sort(array(start:finish))
            endif
        endif

    end subroutine Merge_Sort

    subroutine Merge_ordered_Lists(array, flock)

    ! Merge_Ordered_Lists constructs a table 'Bounds' which tracks the lower and upper bounds of each array fragment. This is used
    ! guide merge sorts on pairs of the collated sorted array partitions in roughly accending size.

        real(dp), intent(inout), dimension(:) :: array
        integer, intent(in) :: flock

        integer :: partition, lower, middle, upper
        integer, dimension(:,:), allocatable :: Bounds

        integer :: i, indx_1, indx_2

        if (flock == 1) return

        partition = size(array)/flock

        allocate(Bounds(flock, 2))

        do i = 1, flock
            Bounds(i, : ) = [(i -1)*partition + 1, i*partition]
        enddo

        Bounds(flock, 2) = size(array)

        i = 1

        do
            lower = 0
            middle = 0
            upper = 0

            do while ((lower == 0) .and. (i < flock))

                if (Bounds(i, 1) /= 0) then

                    lower = Bounds(i, 1)
                    middle = Bounds(i, 2)

                    indx_1 = i
                    i = i + 1
                    exit

                else

                    i = i + 1

                endif

            enddo

            if ((lower == lbound(array, 1)) .and. (middle == ubound(array, 1))) exit

            do while ((upper == 0) .and. (i <= flock))

                if (Bounds(i, 2) /= 0) then

                    upper = Bounds(i, 2)
                    
                    indx_2 = i
                    i = i + 1
                    exit

                else

                    i = i + 1

                endif

            enddo

            if ((lower /= 0) .and. (upper /= 0)) then

                !$omp task shared(array, Bounds) firstprivate(upper, middle, lower) untied

                    call Merge(array, lower, middle, upper)

                    Bounds(indx_1, :) = [lower, upper]
                    Bounds(indx_2, :) = 0

                !$omp end task

            else

                !$omp taskwait
                i = 1

            endif

        enddo

    end subroutine Merge_ordered_Lists


end module Parallel_Merge
