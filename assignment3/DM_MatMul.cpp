/*
 * ParMatMul.c
 *
 *  Created on: Mar 12, 2018
 *      Author: hmayani@cs.stonybrook.edu
 */

#include < mpi.h >
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

int main(int argc,char* argv[])
{

  int r = atoi("4");
  int n = pow(2,r);
  int p = 4;

  int **A;
  int **B;
  int **C;


  A = createMatrix(n,1);
  B = createMatrix(n,1);
  C = createMatrix(n,0);

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

