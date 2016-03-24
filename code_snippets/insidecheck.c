#include <stdio.h>
#include <stdbool.h>
#include <mpi.h>

/******************************************************************************/
struct point
{
	float x;
	float y;
};
/******************************************************************************/
bool isInside(struct point a,struct point b,struct point c,struct point d);
bool ccw(struct point a, struct point b, struct point c);
/******************************************************************************/
int main(int argc, char** argv)
{ 
	const int tag = 74775;
	int rankn;
	int np;
	int startv;
	MPI_Comm comm;
//	/*  createing a type for struct  */
//	const int tm = 2;
//	int block[2] = {1,1};
//	MPI_Datatype types[2]= {MPI_FLOAT, MPI_FLOAT};
//	MPI_Datatype mpi_point;
//	MPI_Aint  	 offset[2];
//
//	offset[0] = offsetof(point, x);
//	offset[1] = offsetof(point, y);
//
//	MPI_Type_create_struct(tm, block, offset, tyes, &mpi_point);
//	MPI_Type_commit(&mpi_point);
//	/*  ---------------------------  */
	
	MPI_Init(&argc, &argv);
	comm = MPI_COMM_WORLD;
	MPI_Comm_size(comm, &np);
	MPI_Comm_rank(comm, &rankn);
	struct point a;
	struct point b;
	struct point c;
	struct point d;
	struct point f;
	/* --------------  */
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
	/* --------------  */
	double start = MPI_Wtime();
	if (rankn == 0)
	{
		if (isInsede(a, b, c, f))
			printf("YES\n");
		else
			printf("NO!!!\n");
		if (isInsede(a, b, f, c))
			printf("YES\n");
		else
			printf("NO!!!\n");
		if (isInsede(b, a, c, f))
			printf("YES\n");
		else
			printf("NO!!!\n");
		if (isInsede(b, a, f, c))
			printf("YES\n");
		else
			printf("NO!!!\n");
	
		if (isInsede(a, b, d, f))
			printf("YES\n");
		else
			printf("NO!!!\n");
		if (isInsede(a, b, f, d))
			printf("YES\n");
		else
			printf("NO!!!\n");
		if (isInsede(b, a, d, f))
			printf("YES\n");
		else
			printf("NO!!!\n");
		if (isInsede(b, a, f, d))
			printf("YES\n");
		else
			printf("NO!!!\n");
	
	}
 	double finish = MPI_Wtime();
 	if (!rankn)
 		printf ("Run time: %g" , finish-start); 


	/////////////////////////////////	
	// MPI_Type_free(&mpi_point);
	MPI_Finalize();
	return 0;
}
/******************************************************************************/
bool isInside(struct point a, struct point b, struct point c, struct point origen)
{
	bool acd = ((origen.y - a.y)*(c.x - a.x)) > ((c.y - a.y)*(origen.x - a.x));
	bool bcd = ((origen.y - b.y)*(c.x - b.x)) > ((c.y - b.y)*(origen.x - b.x));
	bool abd = ((origen.y - a.y)*(b.x - a.x)) > ((b.y - a.y)*(origen.x - a.x));
	bool abc = ((c.y - a.y)*(b.x - a.x)) > ((b.y - a.y)*(c.x - a.x));

	return ( acd != bcd )&&( abc != abd);
}
