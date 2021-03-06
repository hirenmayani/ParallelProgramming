#include <cilk/cilk.h>
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include <cilk/reducer_opadd.h>
#include<chrono>
#include<fstream>
#include <iostream>
#include <cilk/cilk_api.h>

using namespace std;

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
	  arr[i] = rand();
    }
  }
return arr;

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

void parPrefixSumNW(int* x, int* y,int* z,int* s, int n){
	if (n==1){
		s[0] +=x[0];
	}
	else{
		int nb2 = int(n/2);
		for(int i=0; i <nb2;i++){
			y[i]=x[2*i+1]+x[2*i];
		}
		parPrefixSumNW(x,y,z,s,nb2);
		for(int i=0; i <n;i++){
			if(i==0){
				s[0] += x[0];
			}
			else if(i%2==0){
				s[i] = s[i] + z[int(i/2)];
			}
			else{
				s[i] = s[i] + z[int((i-1)/2)] + x[i];
			}
		}

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
		cilk_for(int i=0; i <nb2;i++){
			y[i]=x[2*i+1]+x[2*i];
		}
		int* z = parPrefixSum(y, nb2);

		cilk_for(int i=0; i <n;i++){
			if(i==0){
				s[0] = x[0];
			}
			else if(i%2==0){
				s[i] = z[int((i-2)/2)] + x[i];
			}
			else{
				s[i] = z[int((i-1)/2)];
			}
		}
		free(y);
		free(z);
	}
	return s;
}

void PrefixSum(int* arr, int n){
	for(int i=1; i <n;i++){
		arr[i] = arr[i]+arr[i-1];
	}
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

	int dup = 0 ;
	cilk::reducer< cilk::op_add<unsigned long> > sum;

	cilk_for(int i=0;i<n;i++){
		if(arr[q+i]==x){
			*sum+=1;
		//dup+=1;	
		}
	}
	dup  = sum.get_value();
	//printf("%d - %d \n",x,dup);
	int k = q + lt[n-1];
	//arr[k] = x;
	cilk_for (int i=0;i<dup;i++){
	arr[k+i] = x;
	}
	cilk_for(int i=0;i<n;i++){
		if (b[i]<x)
			arr[q+lt[i]-1] = b[i];
		else if(b[i]>x)
			arr[dup+k+gt[i]-1] = b[i];
			
	}
	free(lt);
	free(gt);
	free(b);
	return k;
}

void isort(int* arr, int q, int n)
{
int i, key, j;
   for (i = 1; i < n; i++)
   {
       key = arr[q+i];
       j = i-1;
 
       while (j >= 0 && arr[q+j] > key)
       {
           arr[q+j+1] = arr[q+j];
           j = j-1;
       }
       arr[q+j+1] = key;
   }
}

void parQuick(int* arr,int q, int r, int m){
	int n = r-q+1;
	if (n <= m){
//		qsort (arr+q, n, sizeof(int), compare);
		isort(arr,q,n);
	}else{
		int pI = rand()%(r+1-q)+q;
//		printf("pi=%d r=%d q=%d\n",pI, q,r);
		int x = arr[pI];
		int k = parPartition(arr, q, r, x);
		cilk_spawn parQuick(arr, q, k-1, m);
		parQuick(arr, k+1, r, m);
		cilk_sync;
	}
}

bool samearr(int* arr,int* arr1,int n){
	  for(int i=0;i<n;i++)
	  {
		 if(arr[i] != arr1[1]){
			 return false;
		 }
	  }
	  return true;
}

int main(int argc,char* argv[]){
	  int nr = atoi(argv[1]);
	  int mr = atoi(argv[2]);
	  int n = pow(2,nr);
	  int m = pow(2,mr);
	  int *arr;
	  arr = createArr(n,1);
	  /*
	  int *arr1;
	  arr1 = createArr(n,0);
	  int *y;
	  y = createArr(n, 0);
	  int *z;
	  z = createArr(n, 0);
	  int *s;
	  s = createArr(n, 0);
*/
	  //printf("\nbefore sorting\n");
	  //printArr(arr,n);
   	   time_t t;
 if (0!= __cilkrts_set_param("nworkers",argv[3]))
 {
    printf("Failed to set worker count\n");
    return 1;
 }
	  srand((unsigned) time(&t));

	  auto start = chrono::system_clock::now();

	  parQuick(arr,0,n-1, m);

	  auto end = chrono::system_clock::now();
	  auto elapsedT = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	  auto elapsed = elapsedT.count();

	  cout<< argv[3]<<"," <<nr << "," << mr << "," <<elapsed<<"\n";

	  ofstream myfile ("1b.csv",ios::app);
	  	myfile<< argv[3] << "," <<nr << "," << mr << "," <<elapsed<<"\n";

	  //printf("\nafter sorting\n");
	  //printArr(arr,n);

//	printf("end");
/*

		  PrefixSum(arr, n);
	  	  printArr(arr, n);
		  s  = parPrefixSum(arr,n);
		  printArr(s, n);
		  if (samearr(s,arr,n))
*/
		  if (arraySortedOrNot(arr, n))
			  printf("\n..................... yey................... \n");
		  else
			  printf("\n..................... ney................... \n");

// 	  }
	free(arr);
	  return 0;
}
