#include "mpi.h"
#include <stdio.h>


int** createMatrix(int size,int init)
{
  int i=0;
  int j=0;


  int **mat;
  mat = (int**)calloc(sizeof(int*),size);
  for(int i = 0; i < size; i++)
  {
      mat[i] = (int*)calloc(sizeof(int), size);
  }

  if(init == 0)
  {
    for(i=0;i<size;i++)
    {
    for(j=0;j<size;j++)
    {
          mat[i][j] = 0;

        }
      }
  }
  else
  {
    for(i=0;i<size;i++)
    {
    for(j=0;j<size;j++)
    {
          mat[i][j] = rand()%30;

        }
    }
  }
return mat;

}


void printMat(int** mat,int size)
{
  int i,j;


  for(i=0;i<size;i++)
  {
    for(j=0;j<size;j++)
      printf("%d ",mat[i][j]);
    printf("\n");
    }
}

 
int main(int argc, char *argv[])
{
	int myrank, n = 0, p = 4;
	//int p=4;
	int r = atoi("3");
	n = pow(2, r);

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	int rootp = sqrt(p * 1.0);
	MPI_Request sendreq[rootp], recvreq[rootp];

	int i = floor(myrank / rootp);
	int j = myrank % rootp;
	int nbrp = n / rootp;

    int A[n*n];
	if (myrank == 0) {
	        for (int ii=0; ii<n*n; ii++) {
	            A[ii] = ii;
	        }
	    }

	    int b[nbrp*nbrp];
	    for (int ii=0; ii<nbrp*nbrp; ii++)
	    		b[ii] = 0;

	    MPI_Datatype blocktype;
	    MPI_Datatype blocktype2;

	    MPI_Type_vector(nbrp, nbrp, n, MPI_INT, &blocktype2);
	    MPI_Type_create_resized( blocktype2, 0, sizeof(int), &blocktype);
	    MPI_Type_commit(&blocktype);

	    int disps[rootp*rootp];
	    int counts[rootp*rootp];

	    for (int ii=0; ii<rootp; ii++) {
	        for (int jj=0; jj<rootp; jj++) {
	            disps[ii*rootp+jj] = ii*n*nbrp+jj*nbrp;
	            counts [ii*rootp+jj] = 1;
	        }
	    }

	    MPI_Scatterv(A, counts, disps, blocktype, b, nbrp*nbrp, MPI_CHAR, 0, MPI_COMM_WORLD);

	    /* each proc prints it's "b" out, in order */
	    for (int proc=0; proc<p; proc++) {
	        if (proc == rank) {
	            printf("Rank = %d\n", rank);
	            if (rank == 0) {
	                printf("Global matrix: \n");
	                for (int ii=0; ii<n; ii++) {
	                    for (int jj=0; jj<n; jj++) {
	                        printf("%3d ",(int)a[ii*n+jj]);
	                    }
	                    printf("\n");
	                }
	            }
	            printf("Local Matrix:\n");
	            for (int ii=0; ii<nbrp; ii++) {
	                for (int jj=0; jj<nbrp; jj++) {
	                    printf("%3d ",(int)b[ii*nbrp+jj]);
	                }
	                printf("\n");
	            }
	            printf("\n");
	        }
	        MPI_Barrier(MPI_COMM_WORLD);
	    }

		MPI_Finalize();
		return 0;
}

