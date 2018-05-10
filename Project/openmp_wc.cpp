#include<iostream>
#include<string>
#include<unordered_map>
#include<vector>
#include<cilk/cilk.h>
#include<cilk/reducer.h>
#include <fstream>
#include <cstdint>
#include <chrono>
#include <sstream>
#include <iterator>

//#include<omp.h>
using namespace std;

/*WORD COUNT REDUCER*/
/*template <class mr>
  mr reduce(mr left, mr right)
  {
	  for(auto start=right.begin(),end = right.end();start != end;++start)
		  left[start->first] += start->second;
	  right.clear();
return left;
}*/
template <class mr>
 mr reduce(mr* left, mr* right)
  {
	  for(typename mr::const_iterator start=right->cbegin(),end = right->cend();start != end;++start)
		  (*left)[start->first] += start->second;
	  right->clear();
return *left;  
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
  (rwz:Monoid:omp_out=reduce<Monoid>(&omp_out,&omp_in)) \
//  initializer(omp_priv=m)

#pragma omp parallel for reduction(rwz:m)
  for(InputIterator it=ibegin, ed = iend; it!=ed; ++it)
  {  	  mapper(&m,*it);
  	  	 result= reduce<Monoid>(&result,&m);
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

template<typename Out>
void split(const std::string &s, char delim, Out result) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
std::transform(item.begin(), item.end(), item.begin(), ::tolower);
*(result++) = item;
    }
}
vector<string> readFile(string path)
{
vector<string> words;
	std::ifstream dict_file(path);
	std::string line;
char delim = ' ';
	while(std::getline(dict_file, line))
		{
			split(line,delim, std::back_inserter(words));
		}

return words;
}
int main(int argc,char* argv[])
{
	cout<<"enter file name"<<endl;
string fname = argv[1];
vector<string>words =readFile("corpus");
/*	vector<string> words;
	words.push_back("a");
	words.push_back("b");
	words.push_back("a");
	words.push_back("b");
	words.push_back("b");*/

	unordered_map<string,int>  m1;
	MapFun<  unordered_map<string,int>,string >   mf;
	auto u1 = map_reduce(words.begin(),words.end(),m1,mf);
	cout<<u1["a"]<<endl;
	cout<<u1["b"]<<endl;
}



