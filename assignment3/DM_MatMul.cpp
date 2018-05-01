/*
 * ParMatMul.c
 *
 *  Created on: Mar 12, 2018
 *      Author: hmayani@cs.stonybrook.edu
 */

#include <mpi.h>
#include <thread>
#include <iostream>
#include <chrono>
#include<math.h>
#include <stdlib.h>
#include<stdio.h>
#include<string.h>
#include <deque>
#include <functional>
#include <queue>
#include <set>
#include <iostream>
#include <unistd.h>

using namespace std;

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


void matMulijk(int** X,int** Y,int** Z,int n)
{
      int i,j,k;

      for(i=0;i<n;i++)
      for(j=0;j<n;j++)
      for(k=0;k<n;k++)
        Z[i][j] = Z[i][j]+X[i][k]*Y[k][j];

}


void MM_rotateA_rotateB(int** A, int** B, int** C, int n, int p){
	int numBlocks = sqrt(p*1.0);
	int blockSize = n/p;
	// TODO consider non divisible results
	cout << numBlocks << "\n";
	cout << blockSize << "\n";

}

void send(int** mat, int i, int j, int rank ){
	
}

void matmul(int** Z, int** X, int** Y, int n){
	int i, j, k;

	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			for (k = 0; k < n; k++)
				Z[i][j] = Z[i][j] + X[i][k] * Y[k][j];

}
int main(int argc,char* argv[])
{

	int myrank,n=0,p=4;
	//int p=4;
	int r = atoi("3");
	n = pow(2,r);
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &p );
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
        int rootp = sqrt(p*1.0);
        MPI_Request sendreq[rootp], recvreq[rootp];



//	int coords[2]; /* 2 Dimension topology so 2 coordinates */
//	MPI_Cart_coords(cart_comm,rank,2,coords);

	int i = floor(myrank/rootp);
	int j = myrank%rootp;

	int nbrp = n/rootp;

//	cout << myrank << ": rootp " << rootp << " i: " << i << " j: " << j << " nbrt:" << nbrp << "\n";
	int** A;
	int** B;
	int** C;
	A = createMatrix(nbrp,1);
	B = createMatrix(nbrp,1);
	C = createMatrix(nbrp,0);

	int left = (rootp+j-i)%rootp;
	int up = (rootp-j+i)%rootp;
	int destA = i*rootp + left;
	int destB = up*rootp + j;

	int right = (rootp+j+i)%rootp;
        int down = (rootp+j+i)%rootp;
	int srcA = i*rootp + right;
        int srcB = down*rootp + j;

//	cout << myrank << " left: " << left << " right:" << right << ": sending to: " << destA << "  receive from:" << srcA << "\n";


	if (myrank != destA)
		MPI_Send(A, nbrp*nbrp, MPI_INT, destA, 0 ,MPI_COMM_WORLD);
	if (myrank != srcA)
		MPI_Recv(A, nbrp*nbrp, MPI_INT, srcA, 0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	
        if (myrank != destB)
                MPI_Send(B, nbrp*nbrp, MPI_INT, destB, 0 ,MPI_COMM_WORLD);

        if (myrank != srcB)
                MPI_Recv(B, nbrp*nbrp, MPI_INT, srcB, 0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);

printMat(A,nbrp);
printf("---");
/*
        int** Aik;
        int** Bkj;

        for (int l=0; l <rootp-1; l++){
                int k = (j+i+l-1)%rootp;

		left = (rootp+j-1)%rootp;
	        up = (rootp-j+i-1)%rootp;
        	destA = i*rootp + left;
        	destB = up*rootp + j;

        	right = (rootp+j+1)%rootp;
        	down = (rootp+i+1)%rootp;
        	srcA = i*rootp + right;
        	srcB = down*rootp + j;

		printMat(A,nbrp);

                matmul(C, A, B, nbrp);

}
/*                if(l < rootp){
                        send(Aik, i, left, myrank);
			if (myrank != destA)
               			MPI_Send(A, nbrp*nbrp, MPI_INT, destA, l ,MPI_COMM_WORLD);
        		if (myrank != srcA)
                		MPI_Recv(A, nbrp*nbrp, MPI_INT, srcA, l,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
                if(l < rootp){
cout << myrank << " "<< destB << " "<< srcB << "\n";
                        //send(Bkj, up, j, myrank);
			if (myrank != destB)
                		MPI_Send(B, nbrp*nbrp, MPI_INT, destB, rootp+l ,MPI_COMM_WORLD);

       			if (myrank != srcB)
               			MPI_Recv(B, nbrp*nbrp, MPI_INT, srcB, rootp+l ,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
        }
*/

/*	
	send(B,up,j, myrank);

        int** Aik;
        int** Bkj;

	for (int l=0; l <rootp; l++){
		int k = (j+i+l-1)%rootp;
		int left = (rootp+j-1)%rootp;
		int up = (rootp+i-1)%rootp;

		matmul(C, Aik, Bkj);
		if(l < rootp){
			send(Aik, i, left, myrank);
		}
		if(l < rootp){
			send(Bkj, up, j, myrank);
		}
	}
*/
	/*
 	for (int i=0; i <rootp; i++){//Par
		for (int j=0; j <rootp; j++){//Par
			int left = (rootp+j-i)%rootp;
			int up = (rootp-j+i)%rootp;

			send(A,i,left, myrank);
			send(B,up,j, myrank);
		}
	}

	for (int l=0; l <rootp; l++){
		for (int i=0; i <rootp; i++){ //Par
			for (int j=0; j <rootp; j++){//Par
				int k = (j+i+l-1)%rootp;
				int left = (rootp+j-1)%rootp;
				int up = (rootp+i-1)%rootp;

				matmul(Cij, Aik, Bkj);
				if(l < rootp){
					send(Aik, i, left, myrank);
				}
				if(l < rootp){
					send(Bkj, up, j, myrank);
				}
			}
		}
	}
*/

	if(myrank  == 0){
/*	    if((n/rootp)* (n/rootp)!=n) {
	    		printf(" of Proc must be equal to %d\nCode terminated",r);
	    		exit(0);
	    }
*/
	}


   sleep(20000);         // wait for 2 seconds before closing


	MPI_Finalize();

/*
  printf("\n-------------Serial---------------\n");
  matMulijk(X,Y,Zs,n);
  printMat(Zs,n);

  printf("\n-------------schedule Parallel---------------\n");
*/
 // MM_rotateA_rotateB(A,B, C, n, p);
  //printMat(Zp,n);

  return 0;
}

