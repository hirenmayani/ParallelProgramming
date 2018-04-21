
//./our 6 1 </work/01905/rezaul/CSE613/HW2/turn-in/roadNet-TX-in.txt
#include<math.h>
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
using namespace std;
struct Edges
{
	int u,v;
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

void parCW_BS(int n, Edges* E,int noe, int* R)
{
	//noe = number of edges
	int* B = createArr(n,0);
	int* l = createArr(n,0);
	int* h = createArr(n,0);
	int* lo = createArr(n,0);
	int* hi = createArr(n,0);
	int* md = createArr(n,0);

	int ks = ceil(log2(noe)) + 1;
#pragma cilk grainsize = 1
	cilk_for(int u=0;u<n;u++){
		l[u]=0;
		h[u]=noe-1;
	}

	for(int k=0;k<ks;k++){
#pragma cilk grainsize = 1
		cilk_for(int u=0;u<n;u++){
			B[u]=0;
			lo[u]= l[u];
			hi[u]=h[u];
		}
#pragma cilk grainsize = 1
		cilk_for(int i=0;i<noe;i++){
			int u = E[i].u;
			md[u] = int(floor((lo[u]+hi[u])/2));
			if (i>=lo[u] and i<= md[u]){
				B[u]=1;
			}
		}
#pragma cilk grainsize = 1

		cilk_for(int i=0;i<noe;i++){
			int u = E[i].u;
			md[u] = int(floor((lo[u]+hi[u])/2));
			if (B[u]== 1 and i>=lo[u] and i<= md[u]){
				h[u] = md[u];
			}else if(B[u]== 0 and i<= hi[u] and i> md[u]){
				l[u] = md[u]+1;
			}
		}

	}
#pragma cilk grainsize = 1
                cilk_for(int i=0;i<noe;i++){
                        int u = E[i].u;
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


void parCountingRank(int* S,int n,int d, int* r,int p)
{
	/*
	 * S - original unsorted array
	 * n - size of array
	 * d - bits in the largest number 2^d-1
	 * r - sorted array container
	 * p - processing elements
	 * */
//	printf("%d received p",p);
	int buckets = pow(2,d);
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
		//if(i==0)
		//	jstart[i] = 0;
		//else
			jstart[i] = (i)*floor(n*1.0/p);
		jend[i] = (i+1<p)?(((i+1)*floor(n*1.0/p))-1):(n-1);
		for(j=jstart[i];j<=jend[i];j++)
//			f[index(S[j], i, buckets)] = f[index(S[j], i, buckets)] + 1;
			f[S[j]][i] = f[S[j]][i] + 1;
		//TODO sync
//		cilk_sync;
	}
printArr(jstart,p);
printArr(jend,p);
	for(j=0;j<buckets;j++)
	{
//				printf("j=%d  \n",j);
//				printArr(f[j],p);
				f[j] = parPrefixSum(f[j],p);
//				 printArr(f[j],p);
			}	
printMat(f,buckets,p);
	//TODO cilk
for(int i=0;i<p;i++)
	{
		ofset[i] = 0;
		for(j=0;j<buckets;j++)
		{
			r1[j][i] = (i==0)?(ofset[i]):(ofset[i] + f[j][i-1]);
			ofset[i] = ofset[i] + f[j][p-1];
		}
		for(j=jstart[i];j<=jend[i];j++)
		{
			r[j] = r1[S[j]][i];
			r1[S[j]][i] = r1[S[j]][i] + 1 ;
		}

	}
printArr(ofset,p);
printMat(f,buckets,p);

free(f,buckets);
free(r1,buckets);
free(jstart,0);
free(jend,0);
free(ofset,0);
}
void parCountingRank1(int* S,int n,int d, int* r,int p)
{
	int buckets = pow(2,d)-1;
//		int b = floor(log2(n))+1;
		int **f = createArr2d(buckets+1,p+1);
		int **r1 = createArr2d(buckets+1,p+1);
		int *jstart =  createArr(p+1,0);
		int *jend = createArr(p+1,0);
		int *ofs = createArr(p+1,0);
		int i=0,j=0;
//#TODO cilk for
		for(int i=1;i<=p;i++)
		{
			for(j=0;j<=buckets;j++)
				f[j][i] = 0;
			jstart[i] = (i-1)*ceil(n/p)+1;
			jend[i] = (i<p)?(i*ceil(n/p)):n;
			for(j=jstart[i];j<=jend[i];j++)
				f[S[j-1]][i] = f[S[j-1]][i] +1;
		}

		for(j=0;j<=buckets;j++)
			f[j] = parPrefixSum(f[j],p);
		//TODO cilkfor
		for(int i=1;i<=p;i++)
		{
			ofs[i] = 1;
			for(j=0;j<=buckets;j++)
			{
				r1[j][i] = (i==1)?ofs[i]:(ofs[i]+f[j][i-1]);
				ofs[i] = ofs[i]+f[j][p];
			}
			for(j=jstart[i];j<=jend[i];j++)
			{
				r[j-1] = r1[S[j-1]][i];
				r1[S[j-1]][i] = r1[S[j-1]][i] + 1;
			}
		}



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
//	printf("\nbs=%d,d=%d,b=%d\n",bucket_size,d,b);
	if(bucket_size<=0)
	{
		parCountingRank(A,n,b, r,p);
//printf("recievd ranking causing sigsev");
//printArr(r,n);
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
//printf("inside function");
//printArr(A,n);
//	free(S,0);
//	free(r,0);
//	free(B,0);
}

void par_PCW_RS(int n, Edges* edges,int noe, int* R)
{
	//noe = number of edges
	int* A = createArr(noe,0);
//changed TODO
	int k = ceil(log2(noe)) ;
	int u,j;
//printf("k=%dbits for edges\n",k);
//TODO cilk for
	cilk_for(int i=0;i<noe;i++)
		A[i] = (edges[i].u<<k)+i;
	//printArr(A,noe);
	parRadixSort(A,noe,1+k+ceil(log2(n)));
//printf("\nafter radix sort");
//printArr(A,noe);
//TODO cilk for
	cilk_for(int i=0;i<noe;i++)
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
//struct Edges
//{
//	int u,v;
//};
void printEdges(Edges* edges,int size)
{
	for(int i=0;i<size;i++)
		printf("\nu=%d v=%d w=%lf\n",edges[i].u,edges[i].v,edges[i].w);

}
int parPartition(int* arr,int q, int r, int x){
	int n = r-q+1;
	if(n==1){
		return q;
	}

	int* b;
	int* lt;
	int* gt;
	b = createArr(n, 0);
	lt = createArr(n, 0);
	gt = createArr(n, 0);

	cilk_for(int i=0;i<n;i++){
		b[i] = arr[q+i];
		if (b[i]<x){
			lt[i] = 1;
		}else{
			lt[i] = 0;
		}

		if (b[i]>x){
			gt[i] = 1;
		}else{
			gt[i] = 0;
		}
	}

	lt = parPrefixSum(lt,n);
	gt = parPrefixSum(gt,n);

	int k = q + lt[n-1];
	arr[k] = x;

	cilk_for(int i=0;i<n;i++){
		if (b[i]<x)
			arr[q+lt[i]-1] = b[i];
		else if(b[i]>x)
			arr[k+gt[i]] = b[i];
	}
	return k;
}

void parQuick(int* arr,int q, int r){
	int n = r-q+1;
	if (n <= 3){
		qsort (arr+q, n, sizeof(int), compare);

	}else{
		int pI = rand()%(r+1-q)+q;
//		printf("pi=%d r=%d q=%d\n",pI, q,r);
		int x = arr[pI];
		int k = parPartition(arr, q, r, x);
		cilk_spawn parQuick(arr, q, k-1);
		parQuick(arr, k+1, r);
		cilk_sync;
	}
}

void mst(int n, Edges* edges, Edges* edgeso,int noe, int* mst)
{
	int* L = createArr(n+1,0);
	int* C = createArr(n,0);
	int* R = createArr(n,0);
	int u1,v1;
int count = 0;
int cp = noe;
//	printEdges(edges,noe);
qsort (edges, noe, sizeof(Edges), compareR);	
//printEdges(edges,noe);
#pragma cilk grainsize = 1
cilk_for(int i=0;i<noe;i++)
{
 edgeso[i].u = edges[i].u;
edgeso[i].v = edges[i].v;
edgeso[i].w = edges[i].w;
}
//printEdges(edges,noe);
#pragma cilk grainsize = 1
	cilk_for(int v=0;v<n;v++)
		L[v] = v;
	bool F = noe>0?true:false;
	std::random_device rd;
std::mt19937 gen(rd());
    	std::bernoulli_distribution dis(0.6);
bool head = true, tail = false;	
bool isEven = true;
while(F)
	{
if(isEven)
{
	std::bernoulli_distribution dis(0.6);
	isEven = false;
}
else
{
std::bernoulli_distribution dis(0.4);
isEven = true;
}
//printf("%d",count);
count+=1;
#pragma cilk grainsize = 1
		cilk_for(int v=0;v<n;v++)
		{
			C[v] = dis(gen);
		}
printf("head tail array");
printArr(C,30);
//		par_PCW_RS(n,edges,noe,R);
	parCW_BS(n,edges,noe,R);
//printArr(R,noe);
//printf("ranking");
#pragma cilk grainsize = 1
		cilk_for(int i=0;i<noe;i++)
		{
			u1 = edges[i].u;
			v1 = edges[i].v;
			//tails - 1
			if( C[u1] == tail && C[v1] == head && R[u1] == i)
			{
				//printf("\ninside if setting u=%d and v=%d",u1,v1);
				L[u1] = v1;
				mst[i] = 1;
			}

		}
//printf("pcw array");
//printArr(L,n);
//printf("selected mst");
//printArr(mst,10);
//printf("\nb4.........");

//printArr(L,n);
#pragma cilk grainsize = 1
		cilk_for(int i=0;i<noe;i++)
		{
			if(edges[i].u == (n)|| edges[i].v == (n))
				continue;
			if(L[edges[i].u]!=L[edges[i].v])
			{
				edges[i].u = L[edges[i].u];
				edges[i].v = L[edges[i].v];
			}
				else
			{
				edges[i].u = n;
                                edges[i].v = n;
				}
		}
//printEdges(edges,noe);
		F = false;
count = 0;
//printf("ht array");
//printArr(C,n);
//printf("after");
//printArr(L,n);
#pragma cilk grainsize = 1
	for(int i=0;i<noe;i++)
		{
			if(edges[i].u!=edges[i].v)
			{
				F = true;
				count += 1;
				}
		}
printf("\nnoe=%d",count);
/*if(count>cp)
{
	printf("number of edges increased");
	printEdges(edges,noe);
	cp = count;
}*/

//printf("\nNumber of remaining edges = %d",count);	
}



}
int main(int argc,char* argv[])
{
///*
int n = atoi(argv[1]);
__cilkrts_set_param("nworkers","64");
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
printArr(sorted,n);
		for(int i=0;i<n;i++)
	                sarr[sorted[i]] = arr[i];
		printArr(sarr,n);
//*/
/*
__cilkrts_set_param("nworkers","64");
printf("please give file number and mode[mode - 0 radix sort mode-1 binary search;]");
	int filen = atoi(argv[1]);
	int mode = atoi(argv[2]);//mode - 0 radix sort mode-1 binary search;
string filenames[] = {"dummy","s-skitter-in.txt",
"com-amazon-in.txt",
"com-friendster-in.txt",
"com-orkut-in.txt",
"roadNet-CA-in.txt",
"roadNet-TX-in.txt",
"ca-AstroPh-in.txt",
"com-dblp-in.txt",
"com-lj-in.txt",
"roadNet-PA-in.txt"};
string fileName = filenames[filen];	
int n,noe;
	scanf("%d %d",&n,&noe);
	printf("\nnumber of vertices = %d\nnumber of edges%d",n,noe);
	Edges* edges = new Edges[noe];
	Edges* edgeso = new Edges[noe];
	int* mstArr = createArr(noe,0);
	for(int i=0;i<noe;i++)
	{	
		scanf("%d %d %lf",&edges[i].u,&edges[i].v,&edges[i].w);
		edges[i].u = edges[i].u-1;
		edges[i].v = edges[i].v-1;
		}	
	//printEdges(edgeso,noe);
//	int* R = createArr(n,0);
//	par_PCW_RS(n,edges,noe,R);
//	printArr(R,n);
//	int* S = createArr(n,1,R,p);
//parCountingRank(S,n,)
	mst(n, edges,edgeso, noe, mstArr);
//printArr(mstArr,noe);
double cost = 0;

cilk_for(int i=0;i<noe;i++)
 if(mstArr[i]==1)
 	cost += edges[i].w;
printf("\n making out file");
ofstream outFile("");
if(mode == 0)
	ofstream outFile (fileName+"-MST-sort-out.txt",ios::out);
else
	ofstream outFile (fileName+"-MST-search-out.txt",ios::out);
outFile<<cost<<endl;
//printEdges(edgeso,noe);
for(int i=0;i<noe;i++)
{
 if(mstArr[i]==1)
 {
 	outFile<<edgeso[i].u<<" "<<edgeso[i].v<<" "<<edgeso[i].w<<endl;
 }
}

outFile.close();
printf("cost=%lf",cost);
printArr(mstArr,noe);
*/
	return 0;
}


