#include<iostream>
#include<string>
#include<unordered_map>
#include<cilk/cilk.h>
#include<cilk/reducer.h>
using namespace std;

/*
 * A note on template classes.
 * Template classes define operations on a generic type t
 * here we make our own template class and we know monoid_base is also a template class.
 * So that our mr system is independent of type
 * mr just a variable name
 * instead of class mr it could have been <int mr>
 * While using our template class mr, the user must specify the type
 *
 * */
template <class mr>
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

		*a note about pointers and ->
			in order to call functions of a class using its pointer use ->
		* What are we doing in reducer
		* for each word in right pointer we are adding its count in the left list where that word is present
	   */
	  for(typename mr::const_iterator start=right->begin();start!=right->end();start++)
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
int main()
{
	unordered_map<string,int> u1;
	unordered_map<string,int> u2;
	u1["a"] = 1;
	u2["b"] = 2;
	u2["a"] = 5;
//	cout<<u1["a"]<<endl;
//	cout<<u1["ab"]<<endl;
//
	//for(auto start=u2.begin();start!=u2.end();start++)
//		u1[start->first] += start->second;
//	cout<<u1["a"];
    Monoid<unordered_map<string,int>> m1;
    m1.reduce(&u1,&u2);
cout<<u1["a"];
}

