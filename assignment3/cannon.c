#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>

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

int main(int argc,char *argv[])
{
   int rank,size,row=0,column=0,count=0,i=0,j=0,k=0;
   MPI_Init(NULL,NULL);
   MPI_Comm_size(MPI_COMM_WORLD,&size);
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);
   if(rank==0)
   {  
    
    count =  pow(2,r); 
    if(count!=size) { printf("No of Proc must be equal to %d\nCode terminated",count); exit(0); }
    fseek( fp, 0, SEEK_SET );
    A=(float*)calloc(sizeof(float),row*column);
    B=(float*)calloc(sizeof(float),row*column);
    k=0;
    printf("A matrix:\n");
    for(i=0;i<row;i++) 
    {
       for(j=0;j<column;j++)
       {
          fscanf(fp,"%f",&n);
          A[k]=n;
          printf("%f\t",A[k]);
          k++; 
       } 
       printf("\n"); 
    }
    fclose(fp);
    k=0;
    printf("\nB matrix:\n");
    fp=fopen("B.txt","r");
    for(i=0;i<row;i++) 
    {
       for(j=0;j<column;j++)
       {
          fscanf(fp,"%f",&n);
          B[k]=n; 
          printf("%f\t",B[k]);
          k++; 
       }
       printf("\n");  
    } 
    fclose(fp);
   }
   MPI_Bcast(&row,1,MPI_INT,0,MPI_COMM_WORLD);
   int periods[]={1,1}; //both vertical and horizontal movement; 
   int dims[]={row,row};
   int coords[2]; /* 2 Dimension topology so 2 coordinates */
   int right=0, left=0, down=0, up=0;    // neighbor ranks
   MPI_Comm cart_comm;
   MPI_Cart_create(MPI_COMM_WORLD,2,dims,periods,1,&cart_comm );
   MPI_Scatter(A,1,MPI_FLOAT,&a,1,MPI_FLOAT,0,cart_comm);
   MPI_Scatter(B,1,MPI_FLOAT,&b,1,MPI_FLOAT,0,cart_comm);
   MPI_Comm_rank(cart_comm,&rank);
   MPI_Cart_coords(cart_comm,rank,2,coords);
   MPI_Cart_shift(cart_comm, 1, coords[0], &left,&right);
   MPI_Cart_shift(cart_comm, 0, coords[1], &up,&down);
   MPI_Sendrecv_replace(&a,1,MPI_FLOAT,left,11,right,11,cart_comm,MPI_STATUS_IGNORE);
   MPI_Sendrecv_replace(&b,1,MPI_FLOAT,up,11,down,11,cart_comm,MPI_STATUS_IGNORE);
   c = c + a*b;
   for(i=1;i<row;i++)
   {
     MPI_Cart_shift(cart_comm, 1, 1, &left,&right);
     MPI_Cart_shift(cart_comm, 0, 1, &up,&down);
     MPI_Sendrecv_replace(&a,1,MPI_FLOAT,left,11,right,11,cart_comm,MPI_STATUS_IGNORE);
     MPI_Sendrecv_replace(&b,1,MPI_FLOAT,up,11,down,11,cart_comm,MPI_STATUS_IGNORE);
     c = c + a*b;
   }
   C=(float*)calloc(sizeof(float),row*row);
   MPI_Gather(&c,1,MPI_FLOAT,C,1,MPI_FLOAT,0,cart_comm);
   if(rank==0)
   {
      k=0; 
      printf("\nA * B:\n");
      for(i=0;i<row;i++)
      {
         for(j=0;j<row;j++)
         {
            printf("%f\t",C[k]);
            k++;
         }    
         printf("\n");
      }
   }
   MPI_Finalize();
   return 0;
}
