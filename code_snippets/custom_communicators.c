
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char** argv) {

    unsigned world_rank = 0, local_rank = 0;
    unsigned world_count = 0, local_count = 0;
    unsigned num_of_groups = 8, group = 0;
    
    //New communicator
    MPI_Comm COMM_LOCAL;
    
    //Basic stuff
    MPI_Init(&argc, &argv);    
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_count);
    
    //Figure out which node go where and their local rank.
    group = world_rank % num_of_groups;
    local_rank = world_rank / num_of_groups;
    
    //Add processors from COMM_WORLD to their respective COMM_LOCAL
    MPI_Comm_split(MPI_COMM_WORLD, group, local_rank, &COMM_LOCAL);
    
    //Get number of processors in each communicator
    MPI_Comm_size(COMM_LOCAL, &local_count);
    
    //use if statement to control which processors print output
    //if(group == 2), if(local_rank == 3), etc
    if(1)
        printf("Group %d of %d, Local rank: %d of %d, Global rank: %d of %d\n", group, num_of_groups, local_rank, local_count, world_rank, world_count);
    
    MPI_Finalize();
    return (EXIT_SUCCESS);
}

