#include<iostream>
#include<string>
#include<unordered_map>
#include<vector>
#include<cilk/cilk.h>
#include<cilk/reducer.h>
#include <typeinfo>
#include <fstream>
#include <cstdint>
#include <stdio.h>
#include <strings.h>
#include <stddef.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <ctype.h>


#define IMG_DATA_OFFSET_POS 10
#define BITS_PER_PIXEL_POS 28
#define CHECK_ERROR(a)                                       \
   if (a)                                                    \
   {                                                         \
      perror("Error at line\n\t" #a "\nSystem Msg");         \
      exit(1);                                               \
   }
int swap1;      // to indicate if we need to swap byte order of header information
void test_endianess() {
   unsigned int num = 0x12345678;
   char *low = (char *)(&(num));
   if (*low ==  0x78) {
      //dprintf("No need to swap\n");
      swap1 = 0;
   }
   else if (*low == 0x12) {
      //dprintf("Need to swap\n");
      swap1 = 1;
   }
   else {
      printf("Error: Invalid value found in memory\n");
      exit(1);
   }
}




void swap_bytes(char *bytes, int num_bytes) {
   int i;
   char tmp;

   for (i = 0; i < num_bytes/2; i++) {
//      //dprintf("Swapping %d and %d\n", bytes[i], bytes[num_bytes - i - 1]);
      tmp = bytes[i];
      bytes[i] = bytes[num_bytes - i - 1];
      bytes[num_bytes - i - 1] = tmp;
   }
}


/*#include "CImg.h"
using namespace cimg_library;
*/
///g++ -I/Users/krishnasharma/Downloads/cilkplus-rtl-src-004516/include mr2.cpp
//icpc -o h.out genericWordCount.cpp -O2 -lm -lpthread -I/usr/X11R6/include -L/usr/X11R6/lib -lm -lpthread -lX11 -std=c++11 -nostartfiles
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

struct map_Monoid:cilk::monoid_base<mr>
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
  static void reduce(value_type *left, value_type *right)
  {
	  for(size_t i=0;i<768;i++)
		  (*left)[i] += (*right)[i];
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

template <typename InputIterator,typename Monoid,class Mapper>
//void __attribute__((flatten))
 typename Monoid::value_type  map_reduce(InputIterator ibegin,InputIterator iend, Monoid m1,Mapper mapper)
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
struct pixel
{
	int arr[3];
};
/*HISTOGRAM MAPPER*/
struct histogram_map
{
	void operator()(pixel pix, uint64_t histogram[768]) const
	{
		histogram[(size_t)pix.arr[0]]++;
		histogram[256+(size_t)pix.arr[1]]++;
		histogram[512+(size_t)pix.arr[2]]++;

	}
};


int main(int argc,char*argv[])
{

	vector<string> words;
	words.push_back("a");
	words.push_back("b");
	words.push_back("a");
	words.push_back("b");
 map_Monoid<unordered_map<string,int> > m1;
MapFun<string,unordered_map<string,int>>  mf;
auto u1 = map_reduce(words.begin(),words.end(),m1,mf);
	cout<<u1["a"];

//	cout<<u1["b"]; 

hist_Monoid m2;
int i;
   int fd;
   char *fdata;
   struct stat finfo;
   char * fname;
   int red[256];
   int green[256];
   int blue[256];


   // Make sure a filename is specified
   if (argv[1] == NULL) {
      printf("USAGE: %s <bitmap filename>\n", argv[0]);
      exit(1);
   }

   fname = argv[1];

   // Read in the file
   CHECK_ERROR((fd = open(fname, O_RDONLY)) < 0);
   // Get the file info (for file length)
   CHECK_ERROR(fstat(fd, &finfo) < 0);
   // Memory map the file
   CHECK_ERROR((fdata = (char*)mmap(0, finfo.st_size + 1,  PROT_READ | PROT_WRITE, MAP_PRIVATE, fd, 0)) == NULL);

   if ((fdata[0] != 'B') || (fdata[1] != 'M')) {
      printf("File is not a valid bitmap file. Exiting\n");
      exit(1);
   }

   test_endianess();    // will set the variable "swap"

   unsigned short *bitsperpixel = (unsigned short *)(&(fdata[BITS_PER_PIXEL_POS]));
   if (swap1) {
      swap_bytes((char *)(bitsperpixel), sizeof(*bitsperpixel));
   }
   if (*bitsperpixel != 24) {    // ensure its 3 bytes per pixel
      printf("Error: Invalid bitmap format - ");
      printf("This application only accepts 24-bit pictures. Exiting\n");
      exit(1);
   }

   unsigned short *data_pos = (unsigned short *)(&(fdata[IMG_DATA_OFFSET_POS]));
   if (swap1) {
      swap_bytes((char *)(data_pos), sizeof(*data_pos));
   }

   int imgdata_bytes = (int)finfo.st_size - (int)(*(data_pos));
   printf("This file has %d bytes of image data, %d pixels\n", imgdata_bytes,
                                                            imgdata_bytes / 3);

   vector<pixel> pixelData;
   pixel pix;
   for (i=*data_pos; i < finfo.st_size; i+=3) {
      unsigned char *val = (unsigned char *)&(fdata[i]);
      pix.arr[0] = *val;
      val = (unsigned char *)&(fdata[i+1]);
      pix.arr[1] = *val;
      val = (unsigned char *)&(fdata[i+2]);
      pix.arr[2] = *val;
      pixelData.push_back(pix);
   }

   /*//dprintf("\n\nBlue\n");
   //dprintf("----------\n\n");
   for (i = 0; i < 256; i++) {
      //dprintf("%d - %d\n", i, blue[i]);
   }

   //dprintf("\n\nGreen\n");
   //dprintf("----------\n\n");
   for (i = 0; i < 256; i++) {
      //dprintf("%d - %d\n", i, green[i]);
   }

   //dprintf("\n\nRed\n");
   //dprintf("----------\n\n");
   for (i = 0; i < 256; i++) {
      //dprintf("%d - %d\n", i, red[i]);
   }
*/
//   cilk_for(auto it=pixelData.begin(), ed = pixelData.begin(); it!=ed; ++it)
//		{
//			cout<<*it;
//
//
//		}
histogram_map mapper;
cilk::reducer<hist_Monoid> redr;
cilk_for(auto it=pixelData.begin(), ed = pixelData.end(); it!=ed; ++it)
{
	pixel pix = *it;
	mapper(pix,redr.view());
}

//auto hist = map_reduce(pixelData.begin(),pixelData.end(),m2,hm);

for(size_t i=0;i<768;i++)
		  cout<<hist[i]<<endl;

}



