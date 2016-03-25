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
	
	point a;
	point b;
	point c;
	point d;
	point f;
	/* --------------  */
	if(rank == 0)
	{
		a.x = 1.0;
		a.y = 7.0;

		b.x = 9.0;
		b.y = 3.0;

		c.x = 7.0;
		c.y = 6.0;

		d.x = 3.0;
		d.y = 4.0;

		f.x = 0.0;
		f.y = 0.0;
	}
	else
	{
		a.x = 0.0;
		a.y = 0.0;

		b.x = 0.0;
		b.y = 0.0;

		c.x = 0.0;
		c.y = 0.0;

		d.x = 0.0;
		d.y = 0.0;

		f.x = 0.0;
		f.y = 0.0;

	}
	/* --------------  */
	double start = MPI_Wtime();

    MPI_Bcast(&a, 1, point_type, 0, comm);
    MPI_Bcast(&b, 1, point_type, 0, comm);
    MPI_Bcast(&c, 1, point_type, 0, comm);
    MPI_Bcast(&d, 1, point_type, 0, comm);
    MPI_Bcast(&f, 1, point_type, 0, comm);

	if (rank == 1)
	{
		if (isInside(a, b, c, f))
			printf("YES\n");
		else
			printf("NO!!!\n");
    	fflush(stdout);
		if (isInside(a, b, f, c))
			printf("YES\n");
		else
			printf("NO!!!\n");
    	fflush(stdout);
		if (isInside(b, a, c, f))
			printf("YES\n");
		else
			printf("NO!!!\n");
    	fflush(stdout);
		if (isInside(b, a, f, c))
			printf("YES\n");
		else
			printf("NO!!!\n");
    	fflush(stdout);
		if (isInside(a, b, d, f))
			printf("YES\n");
		else
			printf("NO!!!\n");
    	fflush(stdout);
		if (isInside(a, b, f, d))
			printf("YES\n");
		else
			printf("NO!!!\n");
    	fflush(stdout);
		if (isInside(b, a, d, f))
			printf("YES\n");
		else
			printf("NO!!!\n");
    	fflush(stdout);
		if (isInside(b, a, f, d))
			printf("YES\n");
		else
			printf("NO!!!\n");
    	fflush(stdout);
	
	}
 	double finish = MPI_Wtime();
 	if (!rank)
 		printf ("Run time: %g" , finish-start); 


	/////////////////////////////////	
	// MPI_Type_free(&mpi_point);
	MPI_Finalize();
	return 0;
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
