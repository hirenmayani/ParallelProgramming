#include<stdio.h>
int extractBitSegment(int value,int left, int right)
{

	int mask = ((1 << (right-left)) - 1) << left;
	int isolatedXbits = (value & mask)>>left;
	return isolatedXbits;
}
int main()
{
	printf("\n%d",extractBitSegment(255,5,6));
	return 0;
}
