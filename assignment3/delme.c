#include "mpi.h"
#include <stdio.h>

int main(int argc, char *argv[])
{
	int myrank, n = 0, p = 4;
	//int p=4;
	int r = atoi("3");
	n = pow(2, r);
	int ii,jj;
	int proc;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	int rootp = sqrt(p * 1.0);
	MPI_Request sendreq[rootp], recvreq[rootp];

	int i ;
	int j ;
	i = floor(myrank / rootp);
	j = myrank % rootp;
	int nbrp = n / rootp;
	int** A;
	int** B;
	int** C;
	A = createContMatrix(nbrp, 0);
	B = createContMatrix(nbrp, 0);
	C = createContMatrix(nbrp, 0);


	int AA[n*n];
	int BB[n*n];
	if (myrank == 0) {
	        for(ii=0; ii<n*n; ii++) {
	            AA[ii] = ii;
	            BB[ii] = ii;
	        }
	    }

		int a[nbrp*nbrp];
		int b[nbrp*nbrp];
	    for(ii=0; ii<nbrp*nbrp; ii++)
	    		{
				a[ii] = 0;
				b[ii] = 0;
			}

	    MPI_Datatype blocktype;
	    MPI_Datatype blocktype2;

	    MPI_Type_vector(nbrp, nbrp, n, MPI_INT, &blocktype2);
	    MPI_Type_create_resized( blocktype2, 0, sizeof(int), &blocktype);
	    MPI_Type_commit(&blocktype);

	    int disps[rootp*rootp];
	    int counts[rootp*rootp];

	    for (ii=0; ii<rootp; ii++) {
	        for (jj=0; jj<rootp; jj++) {
	            disps[ii*rootp+jj] = ii*n*nbrp+jj*nbrp;
	            counts [ii*rootp+jj] = 1;
	        }
	    }

	    MPI_Scatterv(AA, counts, disps, blocktype, a, nbrp*nbrp, MPI_INT, 0, MPI_COMM_WORLD);
	    MPI_Scatterv(BB, counts, disps, blocktype, b, nbrp*nbrp, MPI_INT, 0, MPI_COMM_WORLD);

	    /* each proc prints it's "b" out, in order */
	    for (proc=0; proc<p; proc++) {
	        if (proc == myrank) {
/*	            printf("Rank = %d\n", myrank);
	            if (myrank == 0) {
	                printf("Global matrix: \n");
	                for(ii=0; ii<n; ii++) {
	                    for(jj=0; jj<n; jj++) {
	                        printf("%3d ",(int)AA[ii*n+jj]);
	                    }
	                    printf("\n");
	                }
	            }
*/
	            for(ii=0; ii<nbrp; ii++) {
	                for(jj=0; jj<nbrp; jj++) {
						A[ii][jj]=a[ii*nbrp+jj];
						B[ii][jj]=b[ii*nbrp+jj];
	                }
	            }
	        }
	        MPI_Barrier(MPI_COMM_WORLD);
	    }
        if (5 == myrank) {
			printArr(A, nbrp);
			printArr(B, nbrp);
        }
        MPI_Finalize();
		return 0;
}

