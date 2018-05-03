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

	//int i = floor(myrank / rootp);
	//int j = myrank % rootp;
	int nbrp = n / rootp;

    	int A[n*n];
	if (myrank == 0) {
	        for(ii=0; ii<n*n; ii++) {
	            A[ii] = ii;
	        }
	    }

	    int b[nbrp*nbrp];
	    for(ii=0; ii<nbrp*nbrp; ii++)
	    		b[ii] = 0;

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

	    MPI_Scatterv(A, counts, disps, blocktype, b, nbrp*nbrp, MPI_INT, 0, MPI_COMM_WORLD);

	    /* each proc prints it's "b" out, in order */
	    for (proc=0; proc<p; proc++) {
	        if (proc == myrank) {
	            printf("Rank = %d\n", myrank);
	            if (myrank == 0) {
	                printf("Global matrix: \n");
	                for(ii=0; ii<n; ii++) {
	                    for(jj=0; jj<n; jj++) {
	                        printf("%3d ",(int)A[ii*n+jj]);
	                    }
	                    printf("\n");
	                }
	            }
	            printf("Local Matrix:\n");
	            for(ii=0; ii<nbrp; ii++) {
	                for(jj=0; jj<nbrp; jj++) {
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

