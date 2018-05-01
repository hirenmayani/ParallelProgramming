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

void send(int** mat, int i, int j, int rank){

}

void matmul(int** C, int** A, int** B){

}
int main(int argc,char* argv[])
{

	int myrank,n=0,p=0;
	//int p=4;
	MPI_Init(&argc, &argv);
	int r = atoi("2");
	int n = pow(2,r);
	MPI_Comm_size(MPI_COMM_WORLD, &p );
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

//	int coords[2]; /* 2 Dimension topology so 2 coordinates */
//	MPI_Cart_coords(cart_comm,rank,2,coords);

	int i = floor(myrank/rootp);
	int j = myrank%rootp;

	int rootp = sqrt(p*1.0);
	int nbrp = n/rootp;

	cout << myrank << " out of "<< p << "\n";
	int** A;
	int** B;
	int** C;
	A = createMatrix(nbrp,1);
	B = createMatrix(nbrp,1);
	C = createMatrix(nbrp,0);

	int left = (rootp+j-i)%rootp;
	int up = (rootp-j+i)%rootp;

	send(A,i,left, myrank);
	send(B,up,j, myrank);

	for (int l=0; l <rootp; l++){
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




	MPI_Finalize();

/*
  printf("\n-------------Serial---------------\n");
  matMulijk(X,Y,Zs,n);
  printMat(Zs,n);

  printf("\n-------------schedule Parallel---------------\n");
*/
  MM_rotateA_rotateB(A,B, C, n, p);
  //printMat(Zp,n);

  return 0;
}

