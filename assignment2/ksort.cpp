//#include "PcwRad.h"
#include<math.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<cilk/cilk.h>
#include<cilk/cilk_api.h>

inline int index(int rows, int cols, int m_width)
{
	 return cols + m_width * rows;
}
int compare (const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}

bool arraySortedOrNot(int* arr, int n)
{
    if (n == 0 || n == 1)
        return true;

    for (int i = 1; i < n; i++)
        if (arr[i-1] > arr[i]){
        		printf("%d %d %d",i, arr[i-1],arr[i]);
            return false;
        }
    return true;
}
void free(int** matrix,int row)
{
	for( int i = 0 ; i < row ; i++ )
	{
	    delete[] matrix[i]; // delete array within matrix
	}
	// delete actual matrix
	delete[] matrix;

}
void free(int* matrix,int row)
{
	delete[] matrix;

}

int* createArr(int size,int init)
{
  int i=0;
  int b = ((int)log2(size));
  b = pow(2,b);

  int *arr;
  arr = (int*)calloc(sizeof(int), size);

  if(init == 0)
  {
    for(i=0;i<size;i++)
    {
          arr[i] = 0;
    }
  }
  else
  {
    for(i=0;i<size;i++)
    {
	  arr[i] = rand()%b;
    }
  }

return arr;
}

int** createArr2d(int rowCount,int colCount)//,int b,int init)
{
	int** a = new int*[rowCount];
	for(int i = 0; i < rowCount; ++i)
	    a[i] = new int[colCount];

	/*  int i=0;
	  int j=0;


	  int *mat;
	  mat = (int*)calloc(sizeof(int),rows*cols);


//	  if(init == 0)
//	  {
	    for(i=0;i<rows;i++)
	    {
	    for(j=0;j<cols;j++)
	    {
	          mat[index(i, j, rows)] = 0;


	        }
	      }
	  }
	  else
	  {
	    for(i=0;i<rows;i++)
	    {
	    for(j=0;j<cols;j++)
	    {
	          mat[index(i, j, rows)] = rand()%b;

	        }
	    }

}*/
	  return a;
}

void printArr(int* arr,int size)
{
  int i;
  printf("\n");
  for(i=0;i<size;i++)
  {
      printf("%d ",arr[i]);
  }
  printf("\n");
}
void PrefixSum(int* arr, int n){
	for(int i=1; i <n;i++){
		arr[i] = arr[i]+arr[i-1];
	}
}
void printMat(int** mat,int m,int n)
{
  int i,j;
printf("\n" );

  for(i=0;i<m;i++)
  {
    for(j=0;j<n;j++)
      printf("%d ",mat[i][j]);
    printf("\n");
    }
}

int* parPrefixSum(int* x,int n){
	int* s;

		s= createArr(n,0);
		if (n==1){
			s[0] +=x[0];
		}
		else{
			int* y;
			y = createArr(n, 0);

			int nb2 = int(n/2);
			//TODO
			cilk_for(int i=0; i <nb2;i++){
//				printf(" %d ",(2*i+1));
				y[i]=x[2*i+1]+x[2*i];
			}
			int* z = parPrefixSum(y, nb2);
//TODO
			cilk_for(int i=0; i <n;i++){
				if(i==0){
					s[0] = x[0];
				}
				else if(i%2==0){
//					printf("%d",int((i-2)/2));
					s[i] = z[int((i-2)/2)] + x[i];
				}
				else{
					s[i] = z[int((i-1)/2)];
				}
			}
		}
		return s;
}


void parCountingRank(int* S,int n,int d, int* r,int p)
{
	/*
	 * S - original unsorted array
	 * n - size of array
	 * d - bits in the largest number 2^d-1
	 * r - sorted array container
	 * p - processing elements
	 * */
	printf("%d received p",p);
	int buckets = pow(2,d)-1;
	int b = floor(log2(n))+1;
	int **f = createArr2d(buckets,p);
	int **r1 = createArr2d(buckets,p);
	int *jstart =  createArr(p,0);
	int *jend = createArr(p,0);
	int *ofset = createArr(p,0);
	int i=0,j=0;

	//TODO cilk_for
	cilk_for(int i=0;i<p;i++)
	{
		for(j=0;j<buckets;j++)
//			f[index(j, i, buckets)] = 0;
			f[j][i] = 0;
		if(i==0)
			jstart[i] = 0;
		else
			jstart[i] = (i)*floor(n/p)+1;
		jend[i] = (i+1<p)?((i+1)*floor(n/p)):n-1;
		for(j=jstart[i];j<=jend[i];j++)
//			f[index(S[j], i, buckets)] = f[index(S[j], i, buckets)] + 1;
			f[S[j]][i] = f[S[j]][i] + 1;
		//TODO sync
//		cilk_sync;
	}
	for(j=0;j<buckets;j++)
	{
//				printf("j=%d  \n",j);
//				printArr(f[j],p);
				f[j] = parPrefixSum(f[j],p);
//				 printArr(f[j],p);
			}
//	printMat(f,buckets,p);
	//TODO cilk
	cilk_for(int i=0;i<p;i++)
	{
		ofset[i] = 0;
		for(j=0;j<buckets;j++)
		{
			r1[j][i] = (i==0)?ofset[i]:(ofset[i] + f[j][i-1]);
			ofset[i] = ofset[i] + f[j][p-1];
 //			printf("\nlast processor prefix sum = %d\n",f[j][p-1]);
//			printf("\nf\n");
//			printMat(f,buckets,p);
//			printf("\nr1\n");
//			printMat(r1,buckets,p);
		}
		for(j=jstart[i];j<=jend[i];j++)
		{
//			printArr(r,n);
			r[j] = r1[S[j]][i];
			r1[S[j]][i] = r1[S[j]][i] + 1 ;
		}

	}
//free(f,buckets);
//free(r1,buckets);
//free(jstart,0);
//free(jend,0);
//free(ofset,0);
}

int extractBitSegment(int value,int left, int right)
{

	int mask = ((1 << (right-left)) - 1) << left;
	int isolatedXbits = (value & mask)>>left;
	return isolatedXbits;
}

void parRadixSort(int* A, int n, int b)
{
	int p = __cilkrts_get_nworkers();
	int *S = createArr(n,0);
	int *r = createArr(n,0);
	int *B = createArr(n,0);
	int d = ceil( log2( n/( p*log2(n) ) ) );
	int bucket_size = ceil((b-1)/d); //number of d-bit segments
	int q = 0;
	printf("\nbs=%d,d=%d,b=%d\n",bucket_size,d,b);
	if(bucket_size<=0)
	{
		parCountingRank(A,n,b, r,p);
			cilk_for(int i=0;i<n;i++)
				B[r[i]] = A[i];
		cilk_for(int i=0;i<n;i++)
				A[i] = B[i];
		return;
	}
	for(int k=0;k<bucket_size;k++)
	{
		q = (k+d<=b)?d:b-k;
	    cilk_for(int i=0;i<n;i++)
			S[i] = extractBitSegment(A[i],k,k+q);

//	printArr(S,n);
		parCountingRank(S,n,q+1, r,p);

//	printArr(r,n);
		cilk_for(int i=0;i<n;i++)
			B[r[i]] = A[i];
		cilk_for(int i=0;i<n;i++)
			A[i] = B[i];

	}
printf("inside function");
printArr(A,n);
//	free(S,0);
//	free(r,0);
//	free(B,0);
}
struct Edges
{
	int u,v,w;
};
void par_PCW_RS(int n, Edges* edges,int noe, int* R)
{
	//noe = number of edges
	int* A = createArr(noe,0);
	int k = ceil(log2(noe)) + 1;
	int u,j;
//TODO cilk for
	for(int i=0;i<noe;i++)
		A[i] = (edges[i].u<<k)+i;
	printArr(A,noe);
	parRadixSort(A,noe,k+ceil(log2(n)));
//TODO cilk for
	for(int i=0;i<noe;i++)
	{
		u = A[i]>>k;
		j = A[i] - (u<<k);

	if(i==1||u!=(A[i-1]>>k))
		R[u] = j;
	}
}
//struct Edges
//{
//	int u,v;
//};
void printEdges(Edges* edges,int size)
{
	for(int i=0;i<size;i++)
		printf("\nu=%d v=%d w=%d\n",edges[i].u,edges[i].v,edges[i].w);

}
int main(int argc,char* argv[])
{
	/*int n = 5;
		int p = __cilkrts_get_nworkers();
		printf("PE=%d\n",p);
		int b = floor(log2(n))+1;
		printf("b=%d",b);
		int *arr = createArr(n,1);
		int *sarr = createArr(n,0);
		printf("random array");
		printArr(arr,n);
		int* sorted = createArr(n,0);

		parCountingRank(arr,n,b,sorted,p);

		for(int i=0;i<n;i++)
	                sarr[sorted[i]] = arr[i];
		printArr(sarr,n);
*/
	int n,noe;
//	n = 3;
//	noe = 2;
	scanf("%d %d",&n,&noe);
	printf("\nnumber of vertices = %d\nnumber of edges%d",n,noe);
	Edges* edges = new Edges[noe];

	for(int i=0;i<noe;i++)
		scanf("%d %d %d",&edges[i].u,&edges[i].v,&edges[i].w);
//	edges[0].u = 1;
//	edges[0].v = 2;
//	edges[0].w = 3;
//	edges[1].u = 2;
//	edges[1].v = 1;
//	edges[1].w = 2;
//
	printEdges(edges,noe);
	int* R = createArr(n,0);
	par_PCW_RS(n,edges,noe,R);
printArr(R,n);
	return 0;
}


