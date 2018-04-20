#include <cilk/cilk.h>
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include <cilk/reducer_opadd.h>

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

void parQuick(int* arr,int q, int r){
	int n = r-q+1;
	if (n <= 30){
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
	  int n = pow(2,nr);
	//  int n = 8;

//	  for (n=2;n<200;n++){
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

	  parQuick(arr,0,n-1);
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
