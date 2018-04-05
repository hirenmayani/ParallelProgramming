#include<iostream>
#include<string>
#include<unordered_map>
#include<vector>
#include<cilk/cilk.h>
#include<cilk/reducer.h>
#include <typeinfo>
#include <fstream>
#include <cstdint>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
using namespace cv;

///g++ -I/Users/krishnasharma/Downloads/cilkplus-rtl-src-004516/include mr2.cpp
using namespace std;


/*
 * HashMap monoid
 * A note on template classes.
 * Template classes define operations on a generic type t
 Reducer
 * here we make our own template class and we know monoid_base is also a template class.
 * So that our mr system is independent of type
 * mr just a variable name
 * instead of class mr it could have been <int mr>
 * While using our template class mr, the user must specify the type
 *
 * */
/*WORD COUNT REDUCER*/
template <class mr>
//using Value_Type = typedef value_type<mr>::type;

struct Monoid:cilk::monoid_base<mr>
{
  static void reduce(mr* left, mr* right)
  {
	  /*
	   * A short note on typename
	   * typename specifies the same class as T or its subset or child
	   *
	   * note about iterators:
	   * all STL classes(vectors), have an iterator whcih we can use by
	   * vector<int>::iterator myit;
	   * and assign the beginnig of any vector this iterator and iterate until the iterator reaches the end;
	   * for(myIntVectorIterator = myIntVector.begin(); myIntVectorIterator != myIntVector.end();myIntVectorIterator++)
			{
				cout<<*myIntVectorIterator<<" ";
				//Should output 1 4 8
			}
		* a note about typename
		* just like whatever class of mr is will be replaced
		*
		*a note about pointers and ->
			in order to call functions of a class using its pointer use ->
		* What are we doing in reducer
		* for each word in right pointer we are adding its count in the left list where that word is present
		* -> accesses views, leftmost view accummulator
	   */
	  for(typename mr::const_iterator start=right->cbegin(),end = right->cend();start != end;++start)
		  (*left)[start->first] += start->second;
	  right->clear();
  }

  /*Monoid must define identity
   * I am doubtful why do we need it?
   * We are not returning anything although*/
  static void identity(mr *p)
  {
	  new (p) mr();
  }
};
/*HISTOGTRAM REDUCER*/

struct hist_Monoid:cilk::monoid_base<uint64_t[768]>
{
	typedef uint64_t value_type[768];
  static void reduce(value_type* left, value_type* right)
  {
	  for(size_t i=0;i<768;i++)
		  (*left)[i] = (*right)[i];
  }

  };


/* this class is the backend of mr sys
 * it takes a map func, a reducing monoid and an iterator on the input
 * */
/* A note about attribute
 * It's run when a shared library is loaded, typically during program startup.
That's how all GCC attributes are; presumably to distinguish them from function calls.(like macro or sthg)
here flatten is used for optimizing, so that all functions are inline if possible thus use
-O1 optimization while call.
 *
 * We can most probably get rid of value_type, it just talks about the type of template in Monoid
 * */

template <typename InputIterator,typename Monoid,class MapFun>
//void __attribute__((flatten))
 typename Monoid::value_type  map_reduce(InputIterator ibegin,InputIterator iend, Monoid m1,MapFun mapper)
	{
//cilk::reducer<Monoid<unordered_map<string,int>>> redr;
	cilk::reducer<Monoid> redr;
/*Yey..mapping begins
		 * note that iterators are always pointers
		 * refer unordered map example in main()*/
		cilk_for(InputIterator it=ibegin, ed = iend; it!=ed; ++it)
		{
			/*Note about view
			 * Cilk Plus reducers provide a number of useful properties:
				Each strand has a private view of the reducer,
				so we don't need to use mutexes to serialize access to the reducer.
				 The views are combined by the Cilk runtime by calling the reduce() function of the reducer's
				 monoid when views sync.
				*/
	//redr->insert({{*it,1}});
			mapper(*it,&redr.view());


		}
//		std::swap(op,redr.view());
cilk_sync;
		return redr.view();

	}
/*WORD COUNT MAPPER*/
template <class keys,class Monoid>

class MapFun
{
public:
	void operator()(keys it,Monoid* v) const {
		v->insert({{it,1}});
	    }
};
/*HISTOGRAM MAPPER*/
struct histogram_map
{
	void operator()(const char* pix, uint64_t* histogram[768])
	{
		histogram[(size_t)pix[0]]++;
		histogram[256+(size_t)pix[1]]++;
		histogram[512+(size_t)pix[2]]++;

	}
};


int main()
{
	vector<string> words;
	words.push_back("a");
	words.push_back("b");
	words.push_back("a");
	words.push_back("b");
 Monoid<unordered_map<string,int> > m1;	
MapFun<string,unordered_map<string,int>>  mf;
auto u1 = map_reduce(words.begin(),words.end(),m1,mf);
	cout<<u1["a"];
//	cout<<u1["b"]; 
hist_Monoid m2;
Mat byte_array;
byte_array = imread("poster.jpg", CV_LOAD_IMAGE_COLOR);
cout << "Width : " << byte_array.cols << endl;
cout << "Height: " << byte_array.rows << endl;
auto hist = map_reduce(byte_array,byte_array_len/3,m2,histogram_map);


}



