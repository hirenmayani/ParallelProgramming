#include<math.h>
#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<time.h>
typedef void(*voidReturnType)(int**,int**,int**,int);

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
printf("\n" );

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
      // printf("multiplication ijk done");


}

void matMulikj(int** X,int** Y,int** Z,int n)
{

      int i,j,k;

      for(i=0;i<n;i++)
      for(k=0;k<n;k++)
      for(j=0;j<n;j++)
        Z[i][j] = Z[i][j]+X[i][k]*Y[k][j];
      // printf("multiplication ikj done");
}



void matMuljik(int** X,int** Y,int** Z,int n)
{
      int i,j,k;

      for(j=0;j<n;j++)
      for(i=0;i<n;i++)
      for(k=0;k<n;k++)
        Z[i][j] = Z[i][j]+X[i][k]*Y[k][j];
      // printf("multiplication jik done");
}

void matMuljki(int** X,int** Y,int** Z,int n)
{
      int i,j,k;

      for(j=0;j<n;j++)
      for(k=0;k<n;k++)
      for(i=0;i<n;i++)
        Z[i][j] = Z[i][j]+X[i][k]*Y[k][j];
      // printf("multiplication jki done");
}

void matMulkij(int** X,int** Y,int** Z,int n)
{
      int i,j,k;

      for(k=0;k<n;k++)
      for(i=0;i<n;i++)
      for(j=0;j<n;j++)
      {

        Z[i][j] = Z[i][j]+X[i][k]*Y[k][j];
      }
      // printf("multiplication kij done");
}
void matMulkji(int** X,int** Y,int** Z,int n)
{
      int i,j,k;

      for(k=0;k<n;k++)
      for(j=0;j<n;j++)
      for(i=0;i<n;i++)
        Z[i][j] = Z[i][j]+X[i][k]*Y[k][j];
      // printf("multiplication kji done");
}

void printMatMulTime(int** X,int** Y,int** Z,int n,void(*matmul_fun_pointer)(int**,int**,int**,int))
{
  int t = clock();
  matmul_fun_pointer(X,Y,Z,n);
  t = clock() - t;
  double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
  printf("\nmatMulijk() took %f seconds to execute \n", time_taken);
  // printMat(Z,n);

}


int main(int argc,char* argv[])
{
  int r = atoi(argv[1]);
  int n = pow(2,r);

  int **X;
  int **Y;
  int **Z;
  double time_taken;
  int t=0;

  X = createMatrix(n,1);
  Y = createMatrix(n,1);
  Z = createMatrix(n,0);
  printMatMulTime(X,Y,Z,n,matMulijk);
  Z = createMatrix(n,0);
  printMatMulTime(X,Y,Z,n,&matMulikj);
  Z = createMatrix(n,0);
  printMatMulTime(X,Y,Z,n,matMuljik);
  Z = createMatrix(n,0);
  printMatMulTime(X,Y,Z,n,matMuljki);
  Z = createMatrix(n,0);
  printMatMulTime(X,Y,Z,n,matMulkij);
  Z = createMatrix(n,0);
  printMatMulTime(X,Y,Z,n,matMulkji);


  return 0;
}

