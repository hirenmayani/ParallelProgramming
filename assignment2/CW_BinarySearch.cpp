#include "PcwRad.h"
#include<math.h>
#include<stdlib.h>
#include<stdio.h>

void parCW_BS(int n, Edges* edges,int noe, int* R)
{
	//noe = number of edges
	int* A = createArr(noe,1);
	int k = ceil(log2(noe)) + 1;
	int u,j;
//TODO cilk for
	for(int i=0;i<noe;i++)
		A[i] = edges[i].u<<k+i;
	parRadixSort(A,noe,k+ceil(log2(n)));
//TODO cilk for
	for(int i=0;i<noe;i++)
	{
		u = A[i]>>k;
		j = A[i] - u<<k;

	if(i==1||u!=A[i-1]>>k)
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

	int n,noe;
	scanf("%d %d",&n,&noe);
	printf("\nnumber of vertices = %d\nnumber of edges%d",n,noe);
	Edges* edges = new Edges[n];

	for(int i=0;i<noe;i++)
		scanf("%d %d %d",&edges[i].u,&edges[i].v,&edges[i].w);
	printEdges(edges,noe);
	int* R = new int[noe];
	//par_PCW_RS(n,edges,noe,R);

	return 0;
}


