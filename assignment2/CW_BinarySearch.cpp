#include "PcwRad.h"
#include<math.h>
#include<stdlib.h>
#include<stdio.h>

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

	cilk_for(int u=0;u<n;u++){
		l[u]=0;
		h[u]=noe-1;
	}

	for(int k=0;k<ks;k++){
		cilk_for(int u=0;u<n;u++){
			B[u]=0;
			lo[u]= l[u];
			hi[u]=h[u];
		}

		cilk_for(int i=0;i<noe;i++){
			int u = E[i].u;
			md[u] = int(floor((lo[u]+hi[u])/2));
			if (i>=lo[u] and i<= md[u]){
				B[u]=1;
			}
		}

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
                cilk_for(int i=0;i<noe;i++){
                        int u = E[i].u;
                        if (i == l[u]){
                                R[u]=i;
                        }
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

	int n,noe;
	scanf("%d %d",&n,&noe);
	printf("\nnumber of vertices = %d\nnumber of edges%d",n,noe);
	Edges* edges = new Edges[noe];

	for(int i=0;i<noe;i++)
		scanf("%d %d %d",&edges[i].u,&edges[i].v,&edges[i].w);
	printEdges(edges,noe);
	int* R = new int[n];
	parCW_BS(n,edges,noe,R);

	printArr(R, n);
	return 0;
}


