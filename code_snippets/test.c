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
void collectPoints(point **a, int *, unsigned rank);
void aki(point **a, int *size_a, point **ex, int *size_ex );
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
		collectPoints(&p, &numpoints,  rank);
		int i=0;
		for (i = 0; i < numpoints; ++i)
		{
			printf(" %f -=- %f \n", p[i].x, p[i].y);
		}
	}


	double start = MPI_Wtime();
    MPI_Bcast(&numpoints, 1, MPI_INT, 0, comm);
	
	if (rank != 0)
	{
		p = malloc(numpoints * sizeof(point) );
		printf("rank %d numpoint = %d\n", rank, numpoints);
	}

    MPI_Bcast(p, numpoints, point_type, 0, comm);

	if (rank == 2 )
	{
		printf("rank %d numpoint = %d\n", rank, numpoints);
		printf(" %f -=- %f \n", p[0].x, p[0].y);
		point *mx;
		int size_mx = 0;

		aki(&p, &numpoints, &mx, &size_mx);

	}


 	double finish = MPI_Wtime();
	/* --------------  */

 	if (!rank)
 		printf ("Run time: %g\n" , finish-start); 

	/////////////////////////////////	
	MPI_Type_free(&point_type);
	MPI_Finalize();
	return 0;
}
/******************************************************************************/
void aki(point **a, int *size_a, point **ex, int *size_ex )
{
	int i=0;
	point *temp= *a;
	point xmax = temp[0],
		  ymax = temp[0],
		  xmin = temp[0],
		  ymin = temp[0];

	point	origen;
	origen.x =0;
	origen.y =0;

	for (i=0; i<*size_a; i++)
	{
		if (temp[i].x > xmax.x ) xmax = temp[i];
		if (temp[i].x < xmin.x ) xmin = temp[i];

		if (temp[i].y > ymax.y ) ymax = temp[i];
		if (temp[i].y < ymin.y ) ymin = temp[i];
	}
	int new_size_a	=0;	
	for (i=0; i<*size_a; i++)
	{
		if ( 	!( (temp[i].x == xmax.x) && (temp[i].y == xmax.y) )
			&& 	!( (temp[i].x == ymax.x) && (temp[i].y == ymax.y) )
			&& 	!( (temp[i].x == xmin.x) && (temp[i].y == xmin.y) )
			&& 	!( (temp[i].x == ymin.x) && (temp[i].y == ymin.y) ) )
		{
			if ( (isInside(xmax, ymax, temp[i], origen))
				&&	(isInside(xmin, ymax, temp[i], origen))
				&&	(isInside(xmin, ymin, temp[i], origen))
				&&	(isInside(xmax, ymin, temp[i], origen)) )
			{
				printf ("Inside %d  = %f ; %f \n", i+1 , (*a)[i].x, (*a)[i].y);
			}
			else
				new_size_a++;
		}
			
	}

// 	printf ("X max  = %f ; %f \n", xmax.x, xmax.y);
// 	printf ("X min  = %f ; %f \n", xmin.x, xmin.y);
// 
// 	printf ("Y max  = %f ; %f \n", ymax.x, ymax.y);
// 	printf ("Y min  = %f ; %f \n", ymin.x, ymin.y);

}
/******************************************************************************/
void collectPoints(point **a, int *size,  unsigned rank)
{
	FILE *infile;

	// infile = fopen("../test_cases/64_int_radius_100.txt", "r");
	infile = fopen("../test_cases/10_int_radius_10.txt", "r");
	if(infile == NULL)
		printf ( "Cannot read file !!!");
	int numpoints;
	float x,
		  y;
	if (fscanf(infile, "%d", &numpoints) != 1)
		printf ( "Cannot read file !!!");

	point *arr = malloc(numpoints * sizeof(point) );
	
	int i =0;
	for (i = 0; i < numpoints; ++i)
	{
		if (fscanf (infile, "%f %f", &x, &y) == 2)
		{
			arr[i].x = x;
			arr[i].y = y;
		}
	}
	fclose(infile);
	*size = numpoints;
	*a = arr;
}
/******************************************************************************/
bool isInside(point a, point b, point c, point origen)
{
	bool acd = ((origen.y - a.y)*(c.x - a.x)) > ((c.y - a.y)*(origen.x - a.x));
	bool bcd = ((origen.y - b.y)*(c.x - b.x)) > ((c.y - b.y)*(origen.x - b.x));
	bool abd = ((origen.y - a.y)*(b.x - a.x)) > ((b.y - a.y)*(origen.x - a.x));
	bool abc = ((c.y - a.y)*(b.x - a.x)) > ((b.y - a.y)*(c.x - a.x));

	return !( ( acd != bcd )&&( abc != abd) );
}
