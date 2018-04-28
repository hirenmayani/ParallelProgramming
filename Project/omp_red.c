#define  INT_MIN 0
int mymax(int r,int n) {
  // r is the already reduced value
  // n is the new value
  int m;
  if (n>r) {
    m = n;
  } else {
    m = r;
  }
  return m;
}
int main()
{
	int ndata = 5;
	int data[] = {1,2,3,4,5};
	int m;
#pragma omp declare reduction \
  (rwz:int:omp_out=mymax(omp_out,omp_in)) \
  initializer(omp_priv=INT_MIN)
  m = INT_MIN;
#pragma omp parallel for reduction(rwz:m)
  for (int idata=0; idata<ndata; idata++)
    m = mymax(m,data[idata]);
return 0;
}
