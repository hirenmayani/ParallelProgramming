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

int* createArr(int size,int init)
{
  int i=0;
  int b = ((int)log2(size))+1;

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
	int buckets = pow(2,d-1);
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
				printf("j=%d  \n",j);
//				printArr(f[j],p);  
				f[j] = parPrefixSum(f[j],p);
//				 printArr(f[j],p);
			}	
//	printMat(f,buckets,p);
	//TODO cilk
	for(int i=0;i<p;i++)
	{
		ofset[i] = 1;
		for(j=0;j<buckets;j++)
		{
			r1[j][i] = (i==0)?ofset[i]:(ofset[i] + f[j][i-1]);
			ofset[i] = ofset[i] + f[j][p-1];
 //			printf("\nlast processor prefix sum = %d\n",f[j][p-1]);
//			printf("\nf\n");
//			printMat(f,buckets,p);
//			printf("\nr1\n");
			printMat(r1,buckets,p);
		}
		for(j=jstart[i];j<=jend[i];j++)
		{
			printArr(r,n);
			r[j] = r1[S[j]][i];
			r1[S[j]][i] = r1[S[j]][i] + 1 ;
		}
		
	}
}

int extractBitSegment(int value,int left, int right)
{
	
	int mask = ((1 << (right-left)) - 1) << left;
	int isolatedXbits = value & mask;
	return isolatedXbits;
}

void parRadixSort(int* A, int n, int b,int p)
{
	int *S = createArr(n,0);
	int *r = createArr(n,0);
	int *B = createArr(n,0);
	int d = ceil( log2( n/( p*log2(n) ) ) );
	int bucket_size = (b-1)/d; //number of d-bit segments
	int q = 0;
	for(int k=0;k<bucket_size;k++)
	{
		q = (k+d<=b)?d:b-k
	    cilk_for(int i=0;i<n;i++)
			S[i] = extractBitSegment(A[i],k,k+q-1);
			
		parCountingRank(S,n,q, r,p);
		
		cilk_for(int i=0;i<n;i++)
			B[r[i]] = A[i];
		cilk_for(int i=0;i<n;i++)
			A[i] = B[i];
			
	}
	
	
}

int main(int argc,char* argv[])
{
	
	int n = atoi(argv[1]);
	int p = __cilkrts_get_nworkers();
	printf("PE=%d\n",p);
	int b = floor(log2(n))+1;
	printf("%d",b);
	int *arr = createArr(n,1);
	printf("random array");
	printArr(arr,n);
	int* sorted = createArr(n,0);
	parCountingRank(arr,n,b,sorted,p);
	printArr(sorted,n);
	parRadixSort(arr, n, b, p);
	printArr(arr,n);
	return 0;
}