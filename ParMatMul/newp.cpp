/*
 *  * ParMatMul.c
 *   *
 *    *  Created on: Mar 12, 2018
 *     *      Author: hmayani@cs.stonybrook.edu
 *      */

#include<papi.h>
#include<fstream>
#include <thread>
#include <iostream>
#include <chrono>
#include<math.h>
#include <stdlib.h>
#include<stdio.h>
#include<string.h>
#include <cilk/cilk.h>
using namespace std;

void matMulijk(int* X,int* Y,int* Z,int n)
{

  for(unsigned int i = 0; i < n; ++i){
   for (unsigned int j = 0; j < n; ++j) {
     for (unsigned int k = 0; k < n; ++k) {
       Z[i*n + j] += X[i*n + k] * Y[k*n + j];
      }
    }
 }
}
void i_parallelImplementkij(int zi, int zj, int xi, int xj, int yi, int yj, int currSize, int* Z, int* X, int* Y,int size)
{

		int n = pow(2,currSize);

#pragma cilk grainsize = 5
  for(unsigned int k = 0; k < n; ++k){
   cilk_for (unsigned int i = 0; i < n; ++i) {
     for (unsigned int j = 0; j < n; ++j) {
    	 Z[(zi+i)*size + j+zj] += X[(xi+i)*size + xj + k] * Y[(yi+k)*size + yj + j];
      }
    }
 }
}

void j_parallelImplementkij(int zi, int zj, int xi, int xj, int yi, int yj, int currSize, int* Z, int* X, int* Y,int size)
{

		int n = pow(2,currSize);

#pragma cilk grainsize = 5
  for(unsigned int k = 0; k < n; ++k){
   for (unsigned int i = 0; i < n; ++i) {
     cilk_for (unsigned int j = 0; j < n; ++j) {
    	 Z[(zi+i)*size + j+zj] += X[(xi+i)*size + xj + k] * Y[(yi+k)*size + yj + j];
      }
    }
 }
}

void ij_parallelImplementkij(int zi, int zj, int xi, int xj, int yi, int yj, int currSize, int* Z, int* X, int* Y,int size)
{

		int n = pow(2,currSize);

#pragma cilk grainsize = 5
  for(unsigned int k = 0; k < n; ++k){
   cilk_for (unsigned int i = 0; i < n; ++i) {
     cilk_for (unsigned int j = 0; j < n; ++j) {
    	 Z[(zi+i)*size + j+zj] += X[(xi+i)*size + xj + k] * Y[(yi+k)*size + yj + j];
      }
    }
 }
}
void i_parallelImplementikj(int zi, int zj, int xi, int xj, int yi, int yj, int currSize, int* Z, int* X, int* Y,int size)
{

		int n = pow(2,currSize);


  #pragma cilk grainsize = 5
  cilk_for(unsigned int i = 0; i < n; ++i){
   for (unsigned int k = 0; k < n; ++k) {
     for (unsigned int j = 0; j < n; ++j) {
    	 Z[(zi+i)*size + j+zj] += X[(xi+i)*size + xj + k] * Y[(yi+k)*size + yj + j];
      }
    }
 }

}

void j_parallelImplementikj(int zi, int zj, int xi, int xj, int yi, int yj, int currSize, int* Z, int* X, int* Y,int size)
{

		int n = pow(2,currSize);

  #pragma cilk grainsize = 5
  for(unsigned int i = 0; i < n; ++i){
   for (unsigned int k = 0; k < n; ++k) {
     cilk_for (unsigned int j = 0; j < n; ++j) {
    	 Z[(zi+i)*size + j+zj] += X[(xi+i)*size + xj + k] * Y[(yi+k)*size + yj + j];
      }
    }
 }

}

void ij_parallelImplementikj(int zi, int zj, int xi, int xj, int yi, int yj, int currSize, int* Z, int* X, int* Y,int size)
{

		int n = pow(2,currSize);


  #pragma cilk grainsize = 5
  cilk_for(unsigned int i = 0; i < n; ++i){
   for (unsigned int k = 0; k < n; ++k) {
     cilk_for (unsigned int j = 0; j < n; ++j) {
    	 Z[(zi+i)*size + j+zj] += X[(xi+i)*size + xj + k] * Y[(yi+k)*size + yj + j];
      }
    }
 }

}


void indexMM(int zi, int zj, int xi, int xj, int yi, int yj, int currSize, int* Z, int* X, int* Y,int size)
{

	int n = pow(2,currSize);

	  for(unsigned int k = 0; k < n; ++k){
	   for (unsigned int i = 0; i < n; ++i) {
	     for (unsigned int j = 0; j < n; ++j) {
		     	Z[(zi+i)*size + j+zj] += X[(xi+i)*size + xj + k] * Y[(yi+k)*size + yj + j];
	      }
	    }
	 }
	}




void parRecMM(int zi, int zj, int xi, int xj, int yi, int yj, int currSize, int* Z, int* X, int* Y, int minSize,int Origsize,void(*matmul_base_pointer)(int , int , int , int , int , int , int , int* , int* , int* ,int))
{
		if(currSize <= minSize){
			matmul_base_pointer(zi, zj, xi, xj, yi, yj, currSize, Z, X,  Y,Origsize);
		}
		else{
			int n = pow(2,currSize-1);

						cilk_spawn parRecMM(zi, zj, xi, xj, yi, yj, currSize -1, Z, X, Y, minSize,Origsize,matmul_base_pointer);
						parRecMM(zi, zj+n, xi, xj, yi, yj+n, currSize - 1, Z, X, Y, minSize,Origsize, matmul_base_pointer);
						parRecMM(zi+n, zj, xi+n, xj, yi, yj, currSize - 1, Z, X, Y, minSize,Origsize, matmul_base_pointer);
						parRecMM(zi+n, zj+n, xi+n, xj, yi, yj+n, currSize - 1, Z, X, Y, minSize,Origsize, matmul_base_pointer);
						cilk_sync;


						cilk_spawn parRecMM(zi, zj, xi, xj+n, yi+n, yj, currSize - 1, Z, X, Y, minSize,Origsize, matmul_base_pointer);
						parRecMM(zi, zj+n, xi, xj+n, yi+n, yj+n, currSize - 1, Z, X, Y, minSize,Origsize,matmul_base_pointer);
						parRecMM(zi+n, zj, xi+n, xj+n, yi+n, yj, currSize - 1, Z, X, Y, minSize,Origsize, matmul_base_pointer);
						parRecMM(zi+n, zj+n, xi+n, xj+n, yi+n, yj+n, currSize - 1, Z, X, Y, minSize,Origsize, matmul_base_pointer);
						cilk_sync;
		}
}


void printMat(int* mat,int size)
{
  int i,j;
printf("\n" );

  for(i=0;i<size;i++)
  {
    for(j=0;j<size;j++)
      printf("%d ",mat[i*size+j]);
    printf("\n");
    }
}
int* createMatrix(int size,int init)
{
  int i=0;
  int j=0;


  int *mat;
  mat = (int*)calloc(sizeof(int),size*size);


  if(init == 0)
  {
    for(i=0;i<size;i++)
    {
    for(j=0;j<size;j++)
    {
          mat[i*size+j] = 0;


        }
      }
  }
  else
  {
    for(i=0;i<size;i++)
    {
    for(j=0;j<size;j++)
    {
          mat[i*size+j] = rand()%10;

        }
    }
  }
return mat;

}


int* printMatMulTime(int r, int* X, int* Y, int minSize,int n,void(*matmul_fun_pointer)(int , int , int , int , int , int , int , int*, int*, int*, int ,int,void(*matmul_base_pointer1)(int , int , int , int , int , int , int , int* , int* , int* ,int )), void(*matmul_base_pointer)(int , int , int , int , int , int , int , int* , int* , int* ,int) )
{
	int retval, EventSet = PAPI_NULL;
	long long values[5]={0,0,0,0,0};
	    unsigned counter;
	    unsigned c;
	    unsigned long fact;
	    unsigned stoppoint;


	        /* Initialize the PAPI library */
	        retval = PAPI_library_init(PAPI_VER_CURRENT);

	        if (retval != PAPI_VER_CURRENT) {
	          fprintf(stderr, "PAPI library init error!\n");
	          exit(1);
	        }

		/* Create the Event Set */
	        if (PAPI_create_eventset(&EventSet) != PAPI_OK)
	            printf ("%s:%d\t ERROR\n", __FILE__, __LINE__);

	        /* Add Total Instructions Executed to our EventSet */
	        if (PAPI_add_event(EventSet, PAPI_L1_TCM) != PAPI_OK)
	            printf ("%s:%d\t ERROR\n", __FILE__, __LINE__);

	        /* Add Total Instructions Executed to our EventSet */
	        if (PAPI_add_event(EventSet, PAPI_L2_TCM) != PAPI_OK)
	            printf ("%s:%d\t ERROR\n", __FILE__, __LINE__);

	        /* Add Total Instructions Executed to our EventSet */
	        if (PAPI_add_event(EventSet, PAPI_L3_TCM) != PAPI_OK)
	            printf ("%s:%d\t ERROR\n", __FILE__, __LINE__);
	        int ret = PAPI_add_event(EventSet, PAPI_L1_DCM);
	        if ( ret!= PAPI_OK)
	                    printf ("%s:%d\t ERROR %d\n", __FILE__, __LINE__,ret);
	        ret = PAPI_add_event(EventSet, PAPI_L2_DCM);
	        if (ret != PAPI_OK)
	                    printf ("%s:%d\t ERROR %d\n", __FILE__, __LINE__,ret);


//	   srand ( time(NULL) );
//	   stoppoint = 2;//50+(rand() % 100);
	/* Do some computation here*/

	     	/* Initialize the PAPI library */
	        retval = PAPI_library_init(PAPI_VER_CURRENT);

	        if (retval != PAPI_VER_CURRENT) {
	          fprintf(stderr, "PAPI library init error!\n");
	          exit(1);
	        }


	        /* Start counting */
	        if (PAPI_start(EventSet) != PAPI_OK)

	        	printf ("%s:%d\t ERROR\n", __FILE__, __LINE__);
	        int* Z = createMatrix(n,0);
	        auto start = std::chrono::system_clock::now();

	        matmul_fun_pointer(0,0,0,0,0,0, r, Z, X, Y, minSize,n,matmul_base_pointer);



	        auto end = std::chrono::system_clock::now();
	        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
	        std::cout << "nano seconds = "<<elapsed.count();
	        auto nns = elapsed.count();
	        elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
	        std::cout << ","<<elapsed.count();

//	        printf("Factorial of %d is %lu\n", counter, fact);
	        if (PAPI_read(EventSet, values) != PAPI_OK)
	                    printf ("%s:%d\t ERROR\n", __FILE__, __LINE__);
	                printf ("\nL1_TCM:\t%lld\t L2_TCM:\t%lld\n\n", values[0], values[1]);
	                /* Do some computation here */
	               // printf("\nL1_DCM %lld\nL2_DCM %lld",values[3],values[4]);

	                if (PAPI_stop(EventSet, values) != PAPI_OK)
	                    printf ("%s:%d\t ERROR\n", __FILE__, __LINE__);

	                ofstream myfile ("time.txt",ios::app);
	                // myfile<<nns<<","<<values[0]<<","<<values[1]<<endl;
	        myfile<<nns<<endl;          
		  myfile.close();

	                return Z;
}


int isEqual(int*X,int*Y,int size)
{
  int i=0;
  int j=0;

  for(i=0;i<size;i++)
  {
  for(j=0;j<size;j++)
  {
        if(X[i*size+j] != Y[i*size+j])
        {

                    printf("wrong..returning 0");
                    return 0;
          }
    }
  }
  return 1;
}



int main(int argc,char* argv[])
    {

      int i=0;
      int j=0;
      int k=0;
      int r = atoi(argv[1]);
      int n = pow(2,r);

      int *X;
      int *Y;
      int *Zs;
      int *Zp;

      X = createMatrix(n,1);
      Y = createMatrix(n,1);
      Zs = createMatrix(n,0);
      printf("\n-------------Serial---------------\n");
      matMulijk(X,Y,Zs,n);
//      printMat(Zs,n);

      Zp = createMatrix(n,0);
      printf("\n-------------Parallel---------------\n");

//      parRecMM(0,0,0,0,0,0, r, Zp, X, Y, ,n);
//      printMat(Zp,n);
//      int minSize = 1;
//      for(i=1;i<r-1;i++)
//    	  	  Zp = printMatMulTime(r,  X,  Y, i,n,&parRecMM);
//      	  isEqual(Zs,Zp,n);
      int minSize = 5;
      printMatMulTime(r,  X,  Y, minSize,n,&parRecMM, &i_parallelImplementkij);
      printMatMulTime(r,  X,  Y, minSize,n,&parRecMM, &j_parallelImplementkij);
      printMatMulTime(r,  X,  Y, minSize,n,&parRecMM, &ij_parallelImplementkij);
      printMatMulTime(r,  X,  Y, minSize,n,&parRecMM, &i_parallelImplementikj);
      printMatMulTime(r,  X,  Y, minSize,n,&parRecMM, &j_parallelImplementikj);
      printMatMulTime(r,  X,  Y, minSize,n,&parRecMM, &ij_parallelImplementikj);
/*
//matMulijk(X,Y,Zs,n);
//printMat(Zs,n);
      Zp = printMatMulTime(r,  X,  Y, minSize,n,&parRecMM, &i_parallelImplementkij);
//printMat(Zp,n);      
isEqual(Zs,Zp,n);
      Zp = printMatMulTime(r,  X,  Y, minSize,n,&parRecMM, &j_parallelImplementkij);
      isEqual(Zs,Zp,n);
//printMat(Zp,n);
      Zp = printMatMulTime(r,  X,  Y, minSize,n,&parRecMM, &ij_parallelImplementkij);
      isEqual(Zs,Zp,n);
//printMat(Zp,n);     
 Zp = printMatMulTime(r,  X,  Y, minSize,n,&parRecMM, &i_parallelImplementikj);
      isEqual(Zs,Zp,n);
//printMat(Zp,n);      
Zp = printMatMulTime(r,  X,  Y, minSize,n,&parRecMM, &j_parallelImplementikj);
      isEqual(Zs,Zp,n);
//printMat(Zp,n);      
Zp = printMatMulTime(r,  X,  Y, minSize,n,&parRecMM, &ij_parallelImplementikj);
      isEqual(Zs,Zp,n);
//printMat(Zp,n);
*/
      return 0;
    }

