#include<iostream>
#include<string>
#include<unordered_map>
#include<vector>
#include <typeinfo>
#include <fstream>
#include <cstdint>
//#include<omp.h>
using namespace std;

/*WORD COUNT REDUCER*/
template <class mr>
 void reduce(mr* left, mr* right)
  {
	  for(typename mr::const_iterator start=right->cbegin(),end = right->cend();start != end;++start)
		  (*left)[start->first] += start->second;
	  right->clear();
  }


///*HISTOGTRAM REDUCER*/
//struct hist_Monoid:cilk::monoid_base<uint64_t[768]>
//{
//	typedef uint64_t value_type[768];
//  static void reduce(value_type* left, value_type* right)
//  {
//	  for(size_t i=0;i<768;i++)
//		  (*left)[i] = (*right)[i];
//  }
//
//  };



template <typename InputIterator,typename Monoid,class Mapper>

 Monoid  map_reduce(InputIterator ibegin,InputIterator iend, Monoid m1,Mapper mapper)
	{

 Monoid result;
 Monoid m;
#pragma omp declare reduction \
  (rwz:Monoid:reduce(omp_out,omp_in)) \
  initializer(omp_priv=m)


#pragma omp parallel for reduction(rwz:m)
  for(InputIterator it=ibegin, ed = iend; it!=ed; ++it)
  {  	  mapper(&m,*it);
  	  	  reduce(&result,&m);
  }




  return result;
	}
/*WORD COUNT MAPPER*/
template <class Monoid,class keys>

class MapFun
{
public:
	void operator()(Monoid *v,keys it) const {
		v->insert({{it,1}});

	}

};
struct pixel
{
	int arr[3];
};
/*HISTOGRAM MAPPER*/
struct histogram_map
{
	void operator()(pixel pix, uint64_t* histogram[768]) const
	{
		histogram[(size_t)pix.arr[0]]++;
		histogram[256+(size_t)pix.arr[1]]++;
		histogram[512+(size_t)pix.arr[2]]++;

	}
};


int main()
{
	vector<string> words;
	words.push_back("a");
	words.push_back("b");
	words.push_back("a");
	words.push_back("b");
	words.push_back("b");
 unordered_map<string,int>  m1;
MapFun<  unordered_map<string,int>,string >   mf;
auto u1 = map_reduce(words.begin(),words.end(),m1,mf);
	cout<<u1["a"]<<endl;
	cout<<u1["b"]<<endl;
/*hist_Monoid m2;
CImg<unsigned char> src("poster.jpg");
int width = src.width();
int height = src.height();
vector<pixel> pixelData;
pixel pix;
for (int r = 0; r < height; r++)
        for (int c = 0; c < width; c++){
        		pix.arr[0] = (int)src(c,r,0,0);
			pix.arr[1] = (int)src(c,r,0,1);
			pix.arr[2] = (int)src(c,r,0,2);
			pixelData.push_back(pix);
        }
cout<<width<<endl;
cout<<height<<endl;
#pragma omp parallel for
for(auto it=pixelData.begin(), ed = pixelData.begin(); it!=ed; ++it)
		{
			cout<<*it;


		}
histogram_map hm;
auto hist = map_reduce(pixelData.begin(),pixelData.end(),m2,hm);

for(size_t i=0;i<768;i++)
		  cout<<hist[i];
*/
}



