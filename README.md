# Parallel_Merge_Sort #
OpenMP + MPI implementations of the classic merge sort algorithm, written in Fortran.

 ## Overview ##
An array distributed across an MPI communicator is first merged sorted locally. Those segments are then gathered and merged at the root note. MPI handles the partitioning of each array while openMP is used to parallelize the sort at the individual node level using the OpenMP 'task' construct. This approach is most ideal when developing for an OpenMP + MPI hybrid environment whereby each MPI node constitutes a separate CPU socket. </P>

Parallel_Merge_Indexed.mod sorts a secondary integer array the primary array re-ordering. This is intended for users who wish to keep track of the initial array indices. 

## Build Instructions ##
Parallel_Merge_Sort was developed using Intel compilers. Though, to the best of my knowledge, no Intel specific functionality has been used. To compile the module files and Example.f90 simply type 'mpiifort -fopenmp Parallel_Merge.f90 Parallel_Merge_Sort.f90 Example.f90 -o Example' into your terminal of choice.

## Useage ##
Subroutines are called using an interface which is based on the MPI standard. In particular the input and output structure should be familiar to those who have made use of 'MPI_GatherV'. These subroutines must called from within an initialized MPI environment. 

### Parallel_Merge.mod ###

Parallel_Merge_Sort(local_array, collated_array, root, block_lens, disps, MPI_communicator)

#### Inputs ####
* local_array (double precision real array)
  * section of the MPI-partitioned array
* root (integer)
  * MPI node at which to collate the completed sort
* block_lens (integer array of MPI group size)
  * The number of local_array elements at each node, arranged rank
* disps (integer array of MPI group size)
  * The displacement of each local_array segment from the start of collated_array, arranged by rank
* MPI_communicator
  * the MPI communicator group, typically 'MPI_COMM_WORLD' unless otherwise defined by the user
  
  #### Outputs ####
* collated_array (double precision real array)
  * where the total sorted array will be stored (relevant to the root process only)
 
### Parallel_Merge_Indexed.mod ###

Parallel_Merge_Sort(local_array, local_indices, collated_array, collated_indices root, block_lens, disps, MPI_communicator)

#### Inputs ####
* local_array (double precision real array)
  * section of the MPI-partitioned array
* local_indices (integer array of size local_array)
  * indices or other integers corresponding to the values stored in local_array
* root (integer)
  * MPI node at which to collate the completed sort
* block_lens (integer array of MPI group size)
  * The number of local_array elements at each node, arranged rank
* disps (integer array of MPI group size)
  * The displacement of each local_array segment from the start of collated_array, arranged by rank
* MPI_communicator
  * the MPI communicator group, typically 'MPI_COMM_WORLD' unless otherwise defined by the user
  
  #### Outputs ####
* collated_array (double precision real array)
  * where the total sorted array will be stored (relevant to the root process only)
* collated_indices (integer array of size collated_array)
  * where the total sorted indices will be stored (relevant to the root process only)
