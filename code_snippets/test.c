#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <mpi/mpi.h>

/******************************************************************************/
typedef struct
{
	float x;
	float y;
} point;
/******************************************************************************/
bool isInside(point a, point b, point c, point d);
void collectPoints(point * &, unsigned rank)
/******************************************************************************/
int main(int argc, char** argv)
{ 
    unsigned rank = 0, count = 0, i = 0;
    
    //Basic stuff
    MPI_Init(&argc, &argv);    
	MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &count);

	/*  ---------------------------  */
    //Define point structure for use with MPI communicators
    MPI_Datatype point_type; //Name of derived typed
    MPI_Type_contiguous(2, MPI_FLOAT, &point_type); //Define derived type
    MPI_Type_commit(&point_type); //Derived type must be committed to be used
	/*  ---------------------------  */
	
	point p[5];
	/* --------------  */
	if(rank == 0)
	{
		p[0].x = 1.0;
		p[0].y = 7.0;
		p[1].x = 9.0;
		p[1].y = 3.0;
		p[2].x = 7.0;
		p[2].y = 6.0;
		p[3].x = 3.0;
		p[3].y = 4.0;
		p[4].x = 0.0;
		p[4].y = 0.0;
	}
	else
	{
		p[0].x = 0.0;
		p[0].y = 0.0;
		p[1].x = 0.0;
		p[1].y = 0.0;
		p[2].x = 0.0;
		p[2].y = 0.0;
		p[3].x = 0.0;
		p[3].y = 0.0;
		p[4].x = 0.0;
		p[4].y = 0.0;

	}
	/* --------------  */
	double start = MPI_Wtime();

    MPI_Bcast(&p, 5, point_type, 0, comm);

 	double finish = MPI_Wtime();
 	if (!rank)
 		printf ("Run time: %g" , finish-start); 


	/////////////////////////////////	
	// MPI_Type_free(&mpi_point);
	MPI_Finalize();
	return 0;
}
/******************************************************************************/
void collectPoints(point * &, unsigned rank)
{
	
}
/******************************************************************************/
bool isInside(point a, point b, point c, point origen)
{
	bool acd = ((origen.y - a.y)*(c.x - a.x)) > ((c.y - a.y)*(origen.x - a.x));
	bool bcd = ((origen.y - b.y)*(c.x - b.x)) > ((c.y - b.y)*(origen.x - b.x));
	bool abd = ((origen.y - a.y)*(b.x - a.x)) > ((b.y - a.y)*(origen.x - a.x));
	bool abc = ((c.y - a.y)*(b.x - a.x)) > ((b.y - a.y)*(c.x - a.x));

	return ( acd != bcd )&&( abc != abd);
}
