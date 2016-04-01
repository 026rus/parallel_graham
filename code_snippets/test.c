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
point getOrigen(point x1, point y1, point x2, point y2);
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
point getOrigen(point x1, point y1, point x2, point y2)
{
	point a;
	a.x = 0;
	a.y = 0; 
	return a;
}
void aki(point **a, int *size_a, point **ex, int *size_ex )
{
	int i=0;
	int new_size_a	=0;	
	point *temp;
	point xmax = (*a)[0],
		  ymax = (*a)[0],
		  xmin = (*a)[0],
		  ymin = (*a)[0];

	for (i=0; i<*size_a; i++)
	{
		if ((*a)[i].x > xmax.x ) xmax = (*a)[i];
		if ((*a)[i].x < xmin.x ) xmin = (*a)[i];

		if ((*a)[i].y > ymax.y ) ymax = (*a)[i];
		if ((*a)[i].y < ymin.y ) ymin = (*a)[i];
	}

	point origen =  getOrigen(xmax, ymax, xmin, ymin);

	for (i=0; i<*size_a; i++)
	{
		if ( 	!( ((*a)[i].x == xmax.x) && ((*a)[i].y == xmax.y) )
			&& 	!( ((*a)[i].x == ymax.x) && ((*a)[i].y == ymax.y) )
			&& 	!( ((*a)[i].x == xmin.x) && ((*a)[i].y == xmin.y) )
			&& 	!( ((*a)[i].x == ymin.x) && ((*a)[i].y == ymin.y) ) )
		{
			if ( (isInside(xmax, ymax, (*a)[i], origen))
				&&	(isInside(xmin, ymax, (*a)[i], origen))
				&&	(isInside(xmin, ymin, (*a)[i], origen))
				&&	(isInside(xmax, ymin, (*a)[i], origen)) )
			{
				printf ("Inside %d  = %f ; %f \n", i+1 , (*a)[i].x, (*a)[i].y);
				(*a)[i] = origen;
			}
			else
				new_size_a++;
		}
	}

	 temp = malloc(new_size_a * sizeof(point) );
	 int c=0;
	 for (i=0; i<*size_a; i++)
	 {
	     if ( ((*a)[i].x != origen.x) && ((*a)[i].y != origen.y) )
	     {
	    	temp[c] = (*a)[i];
	    	c++;
	    	if (c == new_size_a) break;
	     }
	 }
	 free (*a);
	 *a= temp;
	 *size_a=new_size_a;
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
