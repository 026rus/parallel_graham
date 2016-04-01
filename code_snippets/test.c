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
// 		int i=0;
//		for (i = 0; i < numpoints; ++i)
//		{
//			printf(" %f -=- %f \n", p[i].x, p[i].y);
//		}
	}


	double start = MPI_Wtime();
    MPI_Bcast(&numpoints, 1, MPI_INT, 0, comm);
	
	if (rank != 0)
	{
		p = malloc(numpoints * sizeof(point) );
		// printf("rank %d numpoint = %d\n", rank, numpoints);
	}

    MPI_Bcast(p, numpoints, point_type, 0, comm);

	if (rank == 2 )
	{
		point *mx;
		int size_mx = 0;

		aki(&p, &numpoints, &mx, &size_mx);

		printf ("Inside %d num max min : %d \n",  numpoints, size_mx);
		for (i=0; i< numpoints; i++)
		{
			printf ("%f \t<=>\t %f\n",  p[i].x, p[i].y);
		}
		// for (i=0; i< size_mx; i++)
		// {
		// 	printf ("\t\t\t Here is max points %d: %f \t<=>\t %f\n",i,  mx[i].x, mx[i].y);
		// }

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
	int new_size_a	=0;	
	*size_ex = 8;
	point *temp;
	point xmax 		= (*a)[0],
		  ymax 		= (*a)[0],
		  xmin 		= (*a)[0],
		  ymin 		= (*a)[0],
		  summax 	= (*a)[0],
		  summin 	= (*a)[0],
		  difmax 	= (*a)[0],
		  difmin 	= (*a)[0];

	for (i=0; i<*size_a; i++)
	{
		if ((*a)[i].x > xmax.x ) xmax = (*a)[i];
		if ((*a)[i].x < xmin.x ) xmin = (*a)[i];

		if ((*a)[i].y > ymax.y ) ymax = (*a)[i];
		if ((*a)[i].y < ymin.y ) ymin = (*a)[i];

		if ( ((*a)[i].x + (*a)[i].y) > (summax.x + summax.y) )  summax = (*a)[i];
		if ( ((*a)[i].x + (*a)[i].y) < (summin.x + summin.y) )  summin = (*a)[i];

		if ( ((*a)[i].x - (*a)[i].y) > (difmax.x - difmax.y) )  difmax = (*a)[i];
		if ( ((*a)[i].x - (*a)[i].y) < (difmin.x - difmin.y) )  difmin = (*a)[i];
	}

	if (*ex == NULL)
		free (*ex);
	*ex = malloc(*size_ex * sizeof(point));
	(*ex)[0] = xmax;
	(*ex)[1] = ymax;
	(*ex)[2] = xmin;
	(*ex)[3] = ymin;
	(*ex)[4] = summax;
	(*ex)[5] = summin;
	(*ex)[6] = difmax;
	(*ex)[7] = difmin;
/* ********************************************  */
	printf ("xmax : %f ; %f\n", xmax.x, xmax.y);
	printf ("ymax : %f ; %f\n", ymax.x, ymax.y);
	printf ("xmin : %f ; %f\n", xmin.x, xmin.y);
	printf ("ymin : %f ; %f\n", ymin.x, ymin.y);
	printf ("summax : %f ; %f\n", summax.x, summax.y);
	printf ("summin : %f ; %f\n", summin.x, summin.y);
	printf ("difmax : %f ; %f\n", difmax.x, difmax.y);
	printf ("difmin : %f ; %f\n", difmin.x, difmin.y);
/* ********************************************  */
	point origen =  getOrigen(xmax, ymax, xmin, ymin);
	printf ("Origen x = %f, \t Y= %f\n", origen.x, origen.y);

	for (i=0; i<*size_a; i++)
	{
		if ( 	!( ((*a)[i].x == xmax.x) && ((*a)[i].y == xmax.y) )
			&& 	!( ((*a)[i].x == ymax.x) && ((*a)[i].y == ymax.y) )
			&& 	!( ((*a)[i].x == xmin.x) && ((*a)[i].y == xmin.y) )
			&& 	!( ((*a)[i].x == ymin.x) && ((*a)[i].y == ymin.y) ) 
			&& 	!( ((*a)[i].x == summax.x) && ((*a)[i].y == summax.y) ) 
			&& 	!( ((*a)[i].x == summin.x) && ((*a)[i].y == summin.y) ) 
			&& 	!( ((*a)[i].x == difmax.x) && ((*a)[i].y == difmax.y) ) 
			&& 	!( ((*a)[i].x == difmin.x) && ((*a)[i].y == difmin.y) ) 
			)
		{
			if ( 	(isInside(xmin, 	difmin, (*a)[i], origen))
				&&	(isInside(difmin, 	ymax, 	(*a)[i], origen))
				&&	(isInside(ymax, 	summax, (*a)[i], origen)) 
				&&	(isInside(summax, 	xmax, 	(*a)[i], origen)) 
				&&	(isInside(xmax, 	difmax, (*a)[i], origen)) 
				&&	(isInside(difmax, 	ymin, 	(*a)[i], origen)) 
				&&	(isInside(ymin, 	summin, (*a)[i], origen)) 
				&&	(isInside(summin, 	xmin, 	(*a)[i], origen)) 
				)
			{
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
	     if (!( ((*a)[i].x == origen.x) && ((*a)[i].y == origen.y) ))
	     {
			if ( 	!( ((*a)[i].x == xmax.x) && ((*a)[i].y == xmax.y) )
				&& 	!( ((*a)[i].x == ymax.x) && ((*a)[i].y == ymax.y) )
				&& 	!( ((*a)[i].x == xmin.x) && ((*a)[i].y == xmin.y) )
				&& 	!( ((*a)[i].x == ymin.x) && ((*a)[i].y == ymin.y) ) )
			{
		    	temp[c] = (*a)[i];
	    		c++;
	    		if (c == new_size_a) break;

	     	}
	 	}
	 }
	 free (*a);
	 *a= temp;
	 *size_a=new_size_a;
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

	infile = fopen("../test_cases/64_int_radius_100.txt", "r");
//	infile = fopen("../test_cases/10_int_radius_10.txt", "r");
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
