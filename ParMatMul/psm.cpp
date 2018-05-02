#include <papi.h>
#include <stdio.h>
#include <stdlib.h>
#include<time.h>
#include<iostream>
#include<chrono>
#include<cilk/cilk.h>
#include<math.h>
using namespace std;
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
                    return 0;
                    printf("wrong..returning 0");
          }
    }
  }
  return 1;
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

void matMulikj(int* X,int* Y,int* Z,int n)
{

  for(unsigned int i = 0; i < n; ++i){
   for (unsigned int k = 0; k < n; ++k) {
     for (unsigned int j = 0; j < n; ++j) {
       Z[i*n + j] += X[i*n + k] * Y[k*n + j];
      }
    }
 }
}
void matMuljik(int* X,int* Y,int* Z,int n)
{

  for(unsigned int j = 0; j < n; ++j){
   for (unsigned int i = 0; i < n; ++i) {
     for (unsigned int k = 0; k < n; ++k) {
       Z[i*n + j] += X[i*n + k] * Y[k*n + j];
      }
    }
 }
}
void matMuljki(int* X,int* Y,int* Z,int n)
{

  for(unsigned int j = 0; j < n; ++j){
   for (unsigned int k = 0; k < n; ++k) {
     for (unsigned int i = 0; i < n; ++i) {
       Z[i*n + j] += X[i*n + k] * Y[k*n + j];
      }
    }
 }
}
void matMulkij(int* X,int* Y,int* Z,int n)
{

  for(unsigned int k = 0; k < n; ++k){
   for (unsigned int i = 0; i < n; ++i) {
     for (unsigned int j = 0; j < n; ++j) {
       Z[i*n + j] += X[i*n + k] * Y[k*n + j];
      }
    }
 }
}

void matMulkji(int* X,int* Y,int* Z,int n)
{

  for(unsigned int k = 0; k < n; ++k){
   for (unsigned int j = 0; j < n; ++j) {
     for (unsigned int i = 0; i < n; ++i) {
       Z[i*n + j] += X[i*n + k] * Y[k*n + j];
      }
    }
 }
}

void i_parallelImplementkij(int* X,int* Y,int* Z,int n)
{
#pragma cilk grainsize = 5
  for(unsigned int k = 0; k < n; ++k){
   cilk_for (unsigned int i = 0; i < n; ++i) {
     for (unsigned int j = 0; j < n; ++j) {
       Z[i*n + j] += X[i*n + k] * Y[k*n + j];
      }
    }
 }
}

void j_parallelImplementkij(int* X,int* Y,int* Z,int n)
{
#pragma cilk grainsize = 5
  for(unsigned int k = 0; k < n; ++k){
   for (unsigned int i = 0; i < n; ++i) {
     cilk_for (unsigned int j = 0; j < n; ++j) {
       Z[i*n + j] += X[i*n + k] * Y[k*n + j];
      }
    }
 }
}
 void ik_parallelImplementkij(int* X,int* Y,int* Z,int n)
 {
 #pragma cilk grainsize = 5
   cilk_for(unsigned int k = 0; k < n; ++k){
    cilk_for (unsigned int i = 0; i < n; ++i) {
      for (unsigned int j = 0; j < n; ++j) {
        Z[i*n + j] += X[i*n + k] * Y[k*n + j];
       }
     }
  }
 }

void ij_parallelImplementkij(int* X,int* Y,int* Z,int n)
{
#pragma cilk grainsize = 5
  for(unsigned int k = 0; k < n; ++k){
   cilk_for (unsigned int i = 0; i < n; ++i) {
     cilk_for (unsigned int j = 0; j < n; ++j) {
       Z[i*n + j] += X[i*n + k] * Y[k*n + j];
      }
    }
 }
}
void i_parallelImplementikj(int* X,int* Y,int* Z,int n)
{

  #pragma cilk grainsize = 5
  cilk_for(unsigned int i = 0; i < n; ++i){
   for (unsigned int k = 0; k < n; ++k) {
     for (unsigned int j = 0; j < n; ++j) {
       Z[i*n + j] += X[i*n + k] * Y[k*n + j];
      }
    }
 }

}

void j_parallelImplementikj(int* X,int* Y,int* Z,int n)
{

  #pragma cilk grainsize = 5
  for(unsigned int i = 0; i < n; ++i){
   for (unsigned int k = 0; k < n; ++k) {
     cilk_for (unsigned int j = 0; j < n; ++j) {
       Z[i*n + j] += X[i*n + k] * Y[k*n + j];
      }
    }
 }

}

 void ik_parallelImplementikj(int* X,int* Y,int* Z,int n)
 {

   #pragma cilk grainsize = 5
   cilk_for(unsigned int i = 0; i < n; ++i){
    cilk_for (unsigned int k = 0; k < n; ++k) {
      for (unsigned int j = 0; j < n; ++j) {
        Z[i*n + j] += X[i*n + k] * Y[k*n + j];
       }
     }
  }

 }
void ij_parallelImplementikj(int* X,int* Y,int* Z,int n)
{

  #pragma cilk grainsize = 5
  cilk_for(unsigned int i = 0; i < n; ++i){
   for (unsigned int k = 0; k < n; ++k) {
     cilk_for (unsigned int j = 0; j < n; ++j) {
       Z[i*n + j] += X[i*n + k] * Y[k*n + j];
      }
    }
 }

}



int* printMatMulTime(int* X,int* Y,int n,void(*matmul_fun_pointer)(int*,int*,int*,int))
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


	   srand ( time(NULL) );
	   stoppoint = 2;//50+(rand() % 100);
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
	        matmul_fun_pointer(X,Y,Z,n);
	        auto end = std::chrono::system_clock::now();
	        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
	        std::cout << ","<<elapsed.count();

 elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
                std::cout << "\nseconds"<<elapsed.count();

//	        printf("Factorial of %d is %lu\n", counter, fact);
	        if (PAPI_read(EventSet, values) != PAPI_OK)
	                    printf ("%s:%d\t ERROR\n", __FILE__, __LINE__);
	                printf ("\nL1_TCM:\t%lld\t L2_TCM:\t%lld\n\n", values[0], values[1]);
	                /* Do some computation here */
	               // printf("\nL1_DCM %lld\nL2_DCM %lld",values[3],values[4]);

	                if (PAPI_stop(EventSet, values) != PAPI_OK)
	                    printf ("%s:%d\t ERROR\n", __FILE__, __LINE__);

return Z;
}
int main(int argc,char* argv[])
{


  int r = atoi(argv[1]);
  int n = pow(2,r);
  int *Y;
  int *X;

	    X = createMatrix(n,1);
	    Y = createMatrix(n,1);
//	    int* Z = createMatrix(n,0);
		 printf("\nserialikj,%d,",n);
int* Z =		 printMatMulTime(X,Y,n,&matMulikj);


		 printf("\n\nComparision begins:serial ikj, parallel ikj ,serial kij, parallel kij ,\n\n");
		 printf("\narallel i,kij,%d,",n);
		 int* Z1 = printMatMulTime(X,Y,n,i_parallelImplementkij);
		 if(!isEqual(Z1,Z,n))
		 		   printf("Wrong matrix multiplication by parallel i_kij");

/*		 printf("\narallel j,kij,%d,",n);
		 Z1 = printMatMulTime(X,Y,n,&j_parallelImplementkij);
		 if(!isEqual(Z1,Z,n))
		 		   printf("Wrong matrix multiplication by parallel j+kij");

		 printf("\narallel ij,kij,%d,",n);
		 printMatMulTime(X,Y,n,&ij_parallelImplementkij);
		 if(!isEqual(Z1,Z,n))
		   printf("Wrong matrix multiplication by parallel ij_kij");
*/
		 printf("\narallel i,ikj,%d,",n);
		 printMatMulTime(X,Y,n,&i_parallelImplementikj);
		 if(!isEqual(Z1,Z,n))
		   printf("Wrong matrix multiplication by parallel i_ikj");
/*
		 printf("\narallel j,ikj,%d,",n);
		 printMatMulTime(X,Y,n,&j_parallelImplementikj);
		 if(!isEqual(Z1,Z,n))
		   printf("Wrong matrix multiplication by serial j+ikj");

		 printf("\nparallel ij,kij,%d,",n);
		 printMatMulTime(X,Y,n,&ij_parallelImplementikj);
		 if(!isEqual(Z1,Z,n))
		   printf("Wrong matrix multiplication by parallel ij+ikj");
*/






//	printf("\nijk,%d,",n);
//	    printMatMulTime(X,Y,n,matMulijk);
//	    printf("\nikj,%d,",n);
//	    printMatMulTime(X,Y,n,matMulikj);
//	    printf("\njik,%d,",n);
//	    printMatMulTime(X,Y,n,matMuljik);
//	    printf("\njki,%d,",n);
//	    printMatMulTime(X,Y,n,matMuljki);
//	    printf("\nkij,%d,",n);
//	    printMatMulTime(X,Y,n,matMulkij);
//	    printf("\nkji,%d,",n);
//	    printMatMulTime(X,Y,n,matMulkji);
//
}
        	/* End of computation */









