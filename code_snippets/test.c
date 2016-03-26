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
void collectPoints(point *a, int *, unsigned rank);
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
	
	point *p;
	int numpoints;
	/* --------------  */
	if(rank == 0)
	{
		collectPoints(p, &numpoints,  rank);
	}


	double start = MPI_Wtime();
    MPI_Bcast(&p, 5, point_type, 0, comm);



 	double finish = MPI_Wtime();
	/* --------------  */

 	if (!rank)
 		printf ("Run time: %g" , finish-start); 

	/////////////////////////////////	
	MPI_Type_free(&point_type);
	MPI_Finalize();
	return 0;
}
/******************************************************************************/
void collectPoints(point *a, int *size,  unsigned rank)
{
	FILE *infile;

	infile = fopen("test_cases/64_int_redius_100.txt", "r");
	if(infile == NULL)
		printf ( "Cannot read file !!!");
	int numpoints;
	float x,
		  y;
	if (fscanf(infile, "%d", &numpoints) != 1)
		printf ( "Cannot read file !!!");

	a = malloc(numpoints * sizeof(point) );
	
	int i =0;
	for (i = 0; i < numpoints; ++i)
	{
		if (fscanf (infile, "%f %f", &x, &y) == 2)
		{
			a[i].x = x;
			a[i].y = y;
		}
	}
	fclose(infile);
	*size = numpoints;
	
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
