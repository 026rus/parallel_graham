
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char** argv) {

    unsigned rank = 0, count = 0, i = 0;
    
    //Basic stuff
    MPI_Init(&argc, &argv);    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &count);
    
    //Define point structure
    typedef struct {
        float x, y;
    } point;
    
    //Define point structure for use with MPI communicators
    MPI_Datatype point_type; //Name of derived typed
    MPI_Type_contiguous(2, MPI_FLOAT, &point_type); //Define derived type
    MPI_Type_commit(&point_type); //Derived type must be committed to be used
    
    point the_point;
    
    //Seed node 0, all other nodes are set to origin
    if(rank == 0) {
        the_point.x = 2.0341;
        the_point.y = 3.3041;
    }
    else {
        the_point.x = 0.0;
        the_point.y = 0.0;    
    }
    
    //Print point values before communication
    if(rank == 0)
        printf("Values before broadcast\n");
    MPI_Barrier(MPI_COMM_WORLD);
    printf("Node: %d of %d, point x: %f y: %f\n", rank, count, the_point.x, the_point.y);
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
    
    //Broadcast the point
    MPI_Bcast(&the_point, 1, point_type, 0, MPI_COMM_WORLD);
    
    //Print point values after communication
    if(rank == 0)
        printf("Values after broadcast\n");
    MPI_Barrier(MPI_COMM_WORLD);
    printf("Node: %d of %d, point x: %f y: %f\n", rank, count, the_point.x, the_point.y);
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
    
    
    MPI_Finalize();
    return (EXIT_SUCCESS);
}