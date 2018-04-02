#include<iostream>
#include<string>
#include<unordered_map>
#include<vector>
#include "/Users/hiren/Downloads/cilkplus-rtl-src-004516/include/cilk/cilk.h"//<cilk/cilk.h>
#include "/Users/hiren/Downloads/cilkplus-rtl-src-004516/include/cilk/reducer.h"//<cilk/reducer.h>

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
	   */
	  for(typename mr::const_iterator start=right->cbegin(),end = right->cend();start != end;++start)
		  (*left)[start->first] += start->second;
	  right->clear();
  }

  /*Monoid must define identity
   * I am doubtful why do we need it?
   * We are not returning anything although*/
  static void identity(mr *p) const
  {
	  new (p) mr();
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
template <class Monoid, class InputIterator, class MapFun>
void __attribute__((flatten))
map_reduce(InputIterator ibegin,InputIterator iend, MapFun mf, typename Monoid::value_type &op)
	{
		cilk::reducer<Monoid> reduce;
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
			mf(*it, reduce.view());

		}
		std::swap(op,reduce.view());

	}

template <class T>
class MapFun
{
public:
	unordered_map<string,int> operator()(typename T::iterator it) const {
		unordered_map<string,int> u1 = (*it,0);
					return u1;

	    }
};

int main()
{
	vector<string> words;
	words.push_back("a");
	words.push_back("b");
	MapFun<string> obj;

	// = {"a","a","b","b","c","abc"};
//	for(auto it=words.begin();it!=words.end();it++)
//		cout<<*it;

//	unordered_map<string,int> u1;
//	unordered_map<string,int> u2;
//	u1["a"] = 1;
//	u2["b"] = 2;
//	u2["a"] = 5;
////	cout<<u1["a"]<<endl;
//	cout<<u1["ab"]<<endl;
//
//	for(auto start=u2.begin();start!=u2.end();start++)
//		u1[start->first] += start->second;
//	cout<<u1["a"];
    Monoid<unordered_map<string,int>> m1;
//    m1.reduce(&u1,&u2);
	map_reduce(words.begin(),words.end(),&obj,m1);
}

