#include<stdio.h>
#include<iostream>
struct Edges
{
	int u,v,w;
};
int compare (const void * a, const void * b)
{
	 Edges *edge1 = (Edges *)a;
	  Edges *edge2 = (Edges *)b;

	  return ( edge1->w - edge2->w );

}
void printEdges(Edges* edges,int size)
{
	for(int i=0;i<size;i++)
		printf("\nu=%d v=%d w=%d\n",edges[i].u,edges[i].v,edges[i].w);

}
int main()
{
	int n,noe;

	scanf("%d %d",&n,&noe);
	printf("\nnumber of vertices = %d\nnumber of edges%d",n,noe);
	Edges* edges = new Edges[noe];
	for(int i=0;i<noe;i++)
		scanf("%d %d %d",&edges[i].u,&edges[i].v,&edges[i].w);
	printEdges(edges,noe);
	//Edges e;
	qsort (edges, noe, sizeof(Edges), compare);
	printEdges(edges,noe);
	return 0;
}


