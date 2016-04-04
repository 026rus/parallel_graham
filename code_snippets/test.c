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
void aki(point **a, int *size_a, point **ex, int *size_ex,unsigned rank, unsigned count, MPI_Comm comm);
point getOrigen(point x1, point y1, point x2, point y2);
bool isNotMPoint(point a, point **ex, int *size_ex);
void aki_outside(point **a, int *size_a, point **ex, int *size_ex);
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
		collectPoints(&p, &numpoints,  rank);

	double start = MPI_Wtime();
    MPI_Bcast(&numpoints, 1, MPI_INT, 0, comm);
	
	if (rank != 0)
	{
		p = malloc(numpoints * sizeof(point) );
		// printf("rank %d numpoint = %d\n", rank, numpoints);
	}

    MPI_Bcast(p, numpoints, point_type, 0, comm);

	point *mx;
	int size_mx = 0;
	aki(&p, &numpoints, &mx, &size_mx, rank, count, comm);

	if (rank == 2 )
	{

		printf ("Inside %d num max min : %d \n",  numpoints, size_mx);
		for (i=0; i< numpoints; i++)
		{
			printf ("%f \t<=>\t %f\n",  p[i].x, p[i].y);
		}
		for (i=0; i< size_mx; i++)
		{
			printf ("\t\t\t Here is max points %d: %f \t<=>\t %f\n",i,  mx[i].x, mx[i].y);
		}

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
void aki(point **a, int *size_a, point **ex, int *size_ex,unsigned rank, unsigned count, MPI_Comm comm)
{
	int i=0;
	int new_size_a	=0;	
	*size_ex = 8;
	point *temp;
	/*  ---------------------------  */
    //Define point structure for use with MPI communicators
    MPI_Datatype point_type; //Name of derived typed
    MPI_Type_contiguous(2, MPI_FLOAT, &point_type); //Define derived type
    MPI_Type_commit(&point_type); //Derived type must be committed to be used
	/*  ---------------------------  */
	/* I think it will be beter to make array here, but for now it eill work */

	/* make sure there is nothing in the array */
	if (*ex == NULL)
		free (*ex);

	/* Allocating space for array */
	*ex = malloc(*size_ex * sizeof(point));

	for (i=0; i<*size_ex; i++)
		(*ex)[i] = (*a)[0];

	/* finding 8 point for aki  */
	for (i=0; i<*size_a; i++)
	{
		if ( (*a)[i].x < (*ex)[0].x ) 								(*ex)[0] = (*a)[i];
		if (((*a)[i].x - (*a)[i].y) < ((*ex)[1].x - (*ex)[1].y) )   (*ex)[1] = (*a)[i];
		if ( (*a)[i].y > (*ex)[2].y ) 								(*ex)[2] = (*a)[i];
		if (((*a)[i].x + (*a)[i].y) > ((*ex)[3].x + (*ex)[3].y) )   (*ex)[3] = (*a)[i];
		if ( (*a)[i].x > (*ex)[4].x ) 								(*ex)[4] = (*a)[i];
		if (((*a)[i].x - (*a)[i].y) > ((*ex)[5].x - (*ex)[5].y) )   (*ex)[5] = (*a)[i];
		if ( (*a)[i].y < (*ex)[6].y ) 								(*ex)[6] = (*a)[i];
		if (((*a)[i].x + (*a)[i].y) < ((*ex)[7].x + (*ex)[7].y) )   (*ex)[7] = (*a)[i];
	}

	/* Brotcast the aki points to all cpus */
    // MPI_Bcast((*ex), *size_ex, point_type, 0, comm);

	/* Here need to Scater points and run rest of the code on separet cpus */
	aki_outside(a, size_a, ex, size_ex);
	
}
/******************************************************************************/
void aki_outside(point **a, int *size_a, point **ex, int *size_ex)
{
	int i=0, 
		new_size_a = 0;
	/* geting center coordinates for the poygon */
	point origen =  getOrigen((*ex)[4], (*ex)[2], (*ex)[0], (*ex)[6]);
	/* Here need to Scater points and run rest of the code on separet cpus */
	// printf ("Origen x = %f, \t Y= %f\n", origen.x, origen.y);

	
	// aki_outside(point **a, int *size_a, point **ex, int *size_ex);
	/* Here need to Scater points and run rest of the code on separet cpus */
	for (i=0; i<*size_a; i++)
	{
		/* Very messy way to make sure I do not checking points to itself */
		if ( isNotMPoint( (*a)[i], ex, size_ex) )
		{
			int c=0;
			bool check = true;
			for (c=0; c<*size_ex; c++)
			{
				if 	(isInside((*ex)[c], (*ex)[c+1], (*a)[i], origen))
					check = true;
				else 
				{
					check = false;
					break;
				}
			}
			if (check) (*a)[i] = origen;
			else new_size_a++;
		}
	}
	/* making room for my points  */
	point *temp = malloc(new_size_a * sizeof(point) );
	int c=0;
	for (i=0; i<*size_a; i++)
	{
	    if (!( ((*a)[i].x == origen.x) && ((*a)[i].y == origen.y) ))
	    {
			if ( isNotMPoint( (*a)[i], ex, size_ex) )
	   		{
	    	   	temp[c] = (*a)[i];
	   			c++;
	   			if (c == new_size_a) break;

	    	}
		}
	}
	/* Here need to gether all points back and point to right pointer */
	/* free spase used by olde points  */
	free (*a);
	*a= temp;
	*size_a=new_size_a;
}
/******************************************************************************/
bool isNotMPoint(point a, point **ex, int *size_ex)
{
		return 	!( (a.x == (*ex)[0].x) && (a.y == (*ex)[0].y) )
			&& 	!( (a.x == (*ex)[1].x) && (a.y == (*ex)[1].y) )
			&& 	!( (a.x == (*ex)[2].x) && (a.y == (*ex)[2].y) )
			&& 	!( (a.x == (*ex)[3].x) && (a.y == (*ex)[3].y) ) 
			&& 	!( (a.x == (*ex)[4].x) && (a.y == (*ex)[4].y) ) 
			&& 	!( (a.x == (*ex)[5].x) && (a.y == (*ex)[5].y) ) 
			&& 	!( (a.x == (*ex)[6].x) && (a.y == (*ex)[6].y) ) 
			&& 	!( (a.x == (*ex)[7].x) && (a.y == (*ex)[7].y) ) ;
}
/******************************************************************************/
point getOrigen(point p1, point p2, point p3, point p4)
{
	point a;
	float A = (p1.x*p2.y-p2.x*p1.y)+
		(p2.x*p3.y-p3.x*p2.y)+
		(p3.x*p4.y-p4.x*p3.y)+
		(p4.x*p1.y-p2.x*p4.y);

	float s1 =  ( (p1.x+p2.x) * (p1.x*p2.y - p2.x*p1.y) ) + 
				( (p2.x+p3.x) * (p2.x*p3.y - p2.x*p3.y) ) + 
				( (p3.x+p4.x) * (p3.x*p4.y - p3.x*p4.y) );

	float s2 =  ( (p1.y+p2.y) * (p1.y*p2.x - p2.y*p1.x) ) + 
				( (p2.y+p3.y) * (p2.y*p3.x - p2.y*p3.x) ) + 
				( (p3.y+p4.y) * (p3.y*p4.x - p3.y*p4.x) );

	a.x = (s1)/ (6*A);
	a.y = (s2)/ (6*A);
	return a;
}
/******************************************************************************/
void collectPoints(point **a, int *size,  unsigned rank)
{
	FILE *infile;

//	infile = fopen("../test_cases/64_int_radius_100.txt", "r");
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
