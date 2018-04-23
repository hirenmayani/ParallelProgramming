//./our 6 1 </work/01905/rezaul/CSE613/HW2/turn-in/roadNet-TX-in.txt
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<cilk/cilk.h>
#include<cilk/cilk_api.h>
#include<iostream>
#include<algorithm>
#include<random>
#include<fstream>
#include<string>
#include<stdint.h>
#include<chrono>
using namespace std;

struct Edges
{
	uint64_t  u,v;
	double w;
};

inline int index(int rows, int cols, int m_width)
{
	 return cols + m_width * rows;
}
int compare (const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}
int compareR (const void * a, const void * b)
{
  return ( ((Edges*)a)->w - ((Edges*)b)->w );
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

bool uarraySortedOrNot(uint64_t * arr, int n)
{
    if (n == 0 || n == 1)
        return true;
    printf("checking");
    for (int i = 1; i < n; i++)
        if (arr[i-1] > arr[i]){
                        printf("%d %d %d",i, arr[i-1],arr[i]);
            return false;
        }
    return true;
}

void free(uint64_t ** matrix,uint64_t  row)
{
	for( uint64_t  i = 0 ; i < row ; i++ )
	{
	    delete[] matrix[i]; // delete array within matrix
	}
	// delete actual matrix
	delete[] matrix;

}
void free(uint64_t * matrix,uint64_t  row)
{
	delete[] matrix;

}

uint64_t * createArr(uint64_t  size,int init)
{
  uint64_t  i=0;
  uint64_t  b = ((uint64_t )log2(size));
  b = pow(2,b);

  uint64_t  *arr;
  arr = (uint64_t *)calloc(sizeof(uint64_t ), size);

  if(init == 0)
  {
    for(i=0;i<size;i++)
    {
          arr[i] = 0;
    }
  }
else if(init == -1)
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

uint64_t ** createArr2d(uint64_t  rowCount,uint64_t  colCount)//,uint64_t  b,uint64_t  init)
{
	uint64_t ** a = new uint64_t *[rowCount];
	for(uint64_t  i = 0; i < rowCount; ++i)
	    a[i] = new uint64_t [colCount];

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

void printArr(uint64_t * arr,uint64_t  size)
{
  uint64_t  i;
  printf("\n");
  for(i=0;i<size;i++)
  {
      printf("%d ",arr[i]);
  }
  printf("\n");
}
void PrefixSum(uint64_t * arr, uint64_t  n){
	for(uint64_t  i=1; i <n;i++){
		arr[i] = arr[i]+arr[i-1];
	}
}
void printMat(uint64_t ** mat,uint64_t  m,uint64_t  n)
{
  uint64_t  i,j;
printf("\n" );

  for(i=0;i<m;i++)
  {
    for(j=0;j<n;j++)
      printf("%d ",mat[i][j]);
    printf("\n");
    }
}

uint64_t * parPrefixSum(uint64_t * x,uint64_t  n){
	uint64_t * s;

		s= createArr(n,0);
		if (n==1){
			s[0] +=x[0];
		}
		else{
			uint64_t * y;
			y = createArr(n, 0);

			uint64_t  nb2 = uint64_t (n/2);
			//TODO
		cilk_for(uint64_t  i=0; i <nb2;i++){
//				printf(" %d ",(2*i+1));
				y[i]=x[2*i+1]+x[2*i];
			}
			uint64_t * z = parPrefixSum(y, nb2);
//TODO
			cilk_for(uint64_t  i=0; i <n;i++){
				if(i==0){
					s[0] = x[0];
				}
				else if(i%2==0){
//					printf("%d",int((i-2)/2));
					s[i] = z[uint64_t ((i-2)/2)] + x[i];
				}
				else{
					s[i] = z[uint64_t ((i-1)/2)];
				}
			}
		}
		return s;
}

void parCW_BS(uint64_t  n, Edges* E,uint64_t  noe, uint64_t * R)
{
	//noe = number of edges
	uint64_t * B = createArr(n,0);
	uint64_t * l = createArr(n,0);
	uint64_t * h = createArr(n,0);
	uint64_t * lo = createArr(n,0);
	uint64_t * hi = createArr(n,0);
	uint64_t * md = createArr(n,0);

	uint64_t  ks = ceil(log2(noe)) + 1;
#pragma cilk grainsize = 1
	cilk_for(uint64_t  u=0;u<n;u++){
		l[u]=0;
		h[u]=noe-1;
	}

	for(uint64_t  k=0;k<ks;k++){
#pragma cilk grainsize = 1
		cilk_for(uint64_t  u=0;u<n;u++){
			B[u]=0;
			lo[u]= l[u];
			hi[u]=h[u];
		}
#pragma cilk grainsize = 1
		cilk_for(uint64_t  i=0;i<noe;i++){
			uint64_t  u = E[i].u;
			md[u] = uint64_t (floor((lo[u]+hi[u])/2));
			if (i>=lo[u] and i<= md[u]){
				B[u]=1;
			}
		}
#pragma cilk grainsize = 1

		cilk_for(uint64_t  i=0;i<noe;i++){
			uint64_t  u = E[i].u;
			md[u] = uint64_t (floor((lo[u]+hi[u])/2));
			if (B[u]== 1 and i>=lo[u] and i<= md[u]){
				h[u] = md[u];
			}else if(B[u]== 0 and i<= hi[u] and i> md[u]){
				l[u] = md[u]+1;
			}
		}

	}
#pragma cilk grainsize = 1
                cilk_for(uint64_t  i=0;i<noe;i++){
                        uint64_t  u = E[i].u;
                        if (i == l[u]){
                                R[u]=i;
                        }
                }

free(B,n);
free(l,n);
free(h,n);
free(lo,n);
free(hi,n);
free(md,n);
}

void parCountingRank(uint64_t * S,uint64_t  n,uint64_t  d, uint64_t * r,uint64_t  p)
{
	/*
	 * S - original unsorted array
	 * n - size of array
	 * d - bits in the largest number 2^d-1
	 * r - sorted array container
	 * p - processing elements
	 * */
//	printf("%d received p",p);
	uint64_t  buckets = pow(2,d);
	uint64_t  b = floor(log2(n))+1;
	uint64_t  **f = createArr2d(buckets,p);
	uint64_t  **r1 = createArr2d(buckets,p);
	uint64_t  *jstart =  createArr(p,0);
	uint64_t  *jend = createArr(p,0);
	uint64_t  *ofset = createArr(p,0);
	uint64_t  i=0;

	//TODO cilk_for
	cilk_for(uint64_t  i=0;i<p;i++)
	{
		for(uint64_t j=0;j<buckets;j++)
			f[j][i] = 0;
		jstart[i] = (i)*floor(n*1.0/p);
		jend[i] = (i+1<p)?(((i+1)*floor(n*1.0/p))-1):(n-1);
		for(uint64_t j=jstart[i];j<=jend[i];j++)
			f[S[j]][i] = f[S[j]][i] + 1;
	}
//uint64_t j =0;
//printArr(jstart,p);
//printArr(jend,p);
//printMat(f,8,p);
	for(uint64_t j=0;j<buckets;j++)
	{
//				printf("j=%d  \n",j);
//				printArr(f[j],p);
				f[j] = parPrefixSum(f[j],p);
//				 printArr(f[j],p);
			}	
//printMat(f,4,p);
	//TODO cilk
	#pragma cilk grainsize = 1
	cilk_for(uint64_t  i=0;i<p;i++)
		{
		ofset[i] = 0;
		for(uint64_t j=0;j<buckets;j++)
		{
			r1[j][i] = (i==0)?(ofset[i]):(ofset[i] + f[j][i-1]);
			ofset[i] = ofset[i] + f[j][p-1];
	//	printf(" %d ",ofset[i]);		
}
		for(uint64_t j=jstart[i];j<=jend[i];j++)
		{
			r[j] = r1[S[j]][i];
			r1[S[j]][i] = r1[S[j]][i] + 1 ;
		}
 //printArr(ofset,p);

	}
//printArr(ofset,p);

free(f,buckets);
free(r1,buckets);
free(jstart,p);
free(jend,p);
free(ofset,p);
}

uint64_t  extractBitSegment(uint64_t  value,uint64_t  left, uint64_t  right)
{

	uint64_t  mask = ((1 << (right-left)) - 1) << left;
	uint64_t  isolatedXbits = (value & mask)>>left;
	return isolatedXbits;
}

void parRadixSort(uint64_t * A, uint64_t  n, uint64_t  b)
{
	uint64_t  p = __cilkrts_get_nworkers();
	uint64_t  *S = createArr(n,0);
	uint64_t  *r = createArr(n,0);
	uint64_t  *B = createArr(n,0);
	uint64_t  d = ceil( log2( n/( p*log2(n) ) ) );
	uint64_t  bucket_size = ceil((b-1)/d); //number of d-bit segments
	uint64_t  q = 0;
//	printf("\nbs=%d,d=%d,b=%d\n",bucket_size,d,b);
	if(bucket_size<=0)
	{
		parCountingRank(A,n,b, r,p);
//printf("recievd ranking causing sigsev");
//printArr(r,n);
			cilk_for(uint64_t  i=0;i<n;i++)
				B[r[i]] = A[i];
		cilk_for(uint64_t  i=0;i<n;i++)
				A[i] = B[i];
		return;
	}
	for(uint64_t  k=0;k<bucket_size;k+=d)
	{
		q = (k+d<=b)?d:b-k;
	    cilk_for(uint64_t  i=0;i<n;i++)
			S[i] = extractBitSegment(A[i],k,k+q);

//	printArr(S,n);
		parCountingRank(S,n,q+1, r,p);

//	printArr(r,n);
		cilk_for(uint64_t  i=0;i<n;i++)
			B[r[i]] = A[i];
		cilk_for(uint64_t  i=0;i<n;i++)
			A[i] = B[i];

	}
//printf("inside function");
//printArr(A,n);
//	free(S,0);
//	free(r,0);
//	free(B,0);
}

void par_PCW_RS(uint64_t  n, Edges* edges,uint64_t  noe, uint64_t * R)
{
	//noe = number of edges
	uint64_t * A = createArr(noe,0);
//changed TODO
	uint64_t  k = ceil(log2(noe)) ;
	uint64_t  u,j;
//printf("k=%dbits for edges\n",k);
//TODO cilk for
	cilk_for(uint64_t  i=0;i<noe;i++)
		A[i] = (edges[i].u<<k)+i;
	//printArr(A,noe);
	parRadixSort(A,noe,1+k+ceil(log2(n)));
//printf("\nafter radix sort");
//printArr(A,noe);
//TODO cilk for
	cilk_for(uint64_t  i=0;i<noe;i++)
	{
		u = A[i]>>k;
		j = A[i] - (u<<k);
		//printf("\nu=%d,j=%d",u,j);
	if(i==0 || u!=(A[i-1]>>k) )
	{
			R[u] = j;
		}
	}
//printf("pcw");
//printArr(R,n);
}

int main(int argc,char* argv[])
{
uint64_t r = atoi(argv[1]);
uint64_t n = pow(2,r);
 if (0!= __cilkrts_set_param("nworkers",argv[2]))
 {
    printf("Failed to set worker count\n");
    return 1;
 }

int p = __cilkrts_get_nworkers();
printf("PE=%d\n",p);

uint64_t b = floor(log2(n));
printf("b=%d",b);
           time_t t;

        srand((unsigned) time(&t));

uint64_t *arr = createArr(n,1);
uint64_t *sarr = createArr(n,0);

auto start = chrono::system_clock::now();
parRadixSort(arr,n,r);
auto end = chrono::system_clock::now();
auto elapsedT = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
auto elapsed = elapsedT.count();


	  ofstream myfile ("1rd.csv",ios::app);
	  	myfile << argv[2] << "," << r << "," <<elapsed<<"\n";
cout <<r<<"," <<elapsed<<"\n";
if (uarraySortedOrNot(sarr, n))
	printf("\n...yey...\n");
else
	printf("\n...ney....\n");

	return 0;
}


