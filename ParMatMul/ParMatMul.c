/*
 * ParMatMul1.c
 *
 *  Created on: Mar 12, 2018
 *      Author: hmayani@cs.stonybrook.edu
 */


#include<math.h>
#include<stdlib.h>
#include<stdio.h>
#include<string.h>
//#include <cilk/cilk.h>



void indexMM(int zi, int zj, int xi, int xj, int yi, int yj, int currSize, int** Z, int** X, int** Y)
{
    int i,j,k;
    int n = pow(2,currSize);
    for(i=0;i<n;i++)
    		for(j=0;j<n;j++)
    			for(k=0;k<n;k++)
    				Z[zi+i][zj+j] += X[xi+i][xj+k]*Y[yi+k][yj+j];

}


void parRecMM(int zi, int zj, int xi, int xj, int yi, int yj, int currSize, int** Z, int** X, int** Y, int minSize)
{
        if(currSize == minSize){
        		indexMM(zi, zj, xi, xj, yi, yj, currSize, Z, X,  Y);
        }
        else{
            int n = pow(2,currSize-1);

			parRecMM(zi, zj, xi, xj, yi, yj, currSize -1, Z, X, Y, minSize);
			parRecMM(zi, zj+n, xi, xj, yi, yj+n, currSize - 1, Z, X, Y, minSize);
			parRecMM(zi+n, zj, xi+n, xj, yi, yj, currSize - 1, Z, X, Y, minSize);
			parRecMM(zi+n, zj+n, xi+n, xj, yi, yj+n, currSize - 1, Z, X, Y, minSize);

			parRecMM(zi, zj, xi, xj+n, yi+n, yj, currSize - 1, Z, X, Y, minSize);
			parRecMM(zi, zj+n, xi, xj+n, yi+n, yj+n, currSize - 1, Z, X, Y, minSize);
			parRecMM(zi+n, zj, xi+n, xj+n, yi+n, yj, currSize - 1, Z, X, Y, minSize);
			parRecMM(zi+n, zj+n, xi+n, xj+n, yi+n, yj+n, currSize - 1, Z, X, Y, minSize);

//        		cilk_spawn parRecMM(zi, zj, xi, xj, yi, yj, currSize -1, Z, X, Y, minSize);
//        		parRecMM(zi, zj*3/2, xi, xj, yi, yj*3/2, currSize - 1, Z, X, Y, minSize);
//        		parRecMM(zi*3/2, zj, xi*3/2, xj, yi, yj, currSize - 1, Z, X, Y, minSize);
//        		parRecMM(zi*3/2, zj*3/2, xi*3/2, xj, yi, yj*3/2, currSize - 1, Z, X, Y, minSize);
//        		cilk_sync;
//
//        		cilk_spawn parRecMM(zi, zj, xi, xj*3/2, yi*3/2, yj, currSize - 1, Z, X, Y, minSize);
//        		parRecMM(zi, zj*3/2, xi, xj*3/2, yi*3/2, yj*3/2, currSize - 1, Z, X, Y, minSize);
//        		parRecMM(zi*3/2, zj, xi*3/2, xj*3/2, yi*3/2, yj, currSize - 1, Z, X, Y, minSize);
//        		parRecMM(zi*3/2, zj*3/2, xi*3/2, xj*3/2, yi*3/2, yj*3/2, currSize - 1, Z, X, Y, minSize);
//        		cilk_sync;

        }
}


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


int main(int argc,char* argv[])
{

  int i=0;
  int j=0;
  int k=0;
  int r = atoi(argv[1]);
  int n = pow(2,r);

  int **X;
  int **Y;
  int **Zs;
  int **Zp;

  int u = 16*3/2;

  X = createMatrix(n,1);
  Y = createMatrix(n,1);
  Zs = createMatrix(n,0);
  printf("\n-------------Serial---------------\n");
  matMulijk(X,Y,Zs,n);
  printMat(Zs,n);

  Zp = createMatrix(n,0);
  printf("\n-------------Parallel---------------\n");
  parRecMM(0,0,0,0,0,0, r, Zp, X, Y, 2);
  printMat(Zp,n);

  return 0;
}

