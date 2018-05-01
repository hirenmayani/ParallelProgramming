#include "CImg.h"
using namespace cimg_library;


int main(){
	CImg<unsigned char> src("a.png");
	int width = src.width();
	int height = src.height();
	//unsigned char* ptr = src.data(10,10); // get pointer to pixel @ 10,10
	//unsigned char pixel = *ptr;
}
