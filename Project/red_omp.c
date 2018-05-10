#include<unordered_map>
#include<vector>
#include<string>
#include<iostream>
using namespace std;
unordered_map<string,int> mymax(unordered_map<string,int>  l,unordered_map<string,int>  r) {
for(auto i=r.begin();i!=r.end();i++)
	l[i->first]+=i->second;
return l;
}
int main()
{
vector<string> words;
	words.push_back("a");
	words.push_back("b");
	words.push_back("a");
	words.push_back("b");
	words.push_back("b");
unordered_map<string,int> m,map;
 m = map;
#pragma omp declare reduction \
  (rwz:unordered_map<string,int> :omp_out=mymax(omp_out,omp_in)) \
//  initializer(omp_priv=m)
#pragma omp parallel for reduction(rwz:m)
  for (int idata=0; idata<5; idata++)
{
unordered_map<string,int> t;
t[words[idata]] = 1;    
m = mymax(m,t);
}
cout<<m["a"];
cout<<m["b"];
return 0;
}
