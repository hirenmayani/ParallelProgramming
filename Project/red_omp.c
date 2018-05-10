#include<iostream>
#include<string>
#include<unordered_map>
#include<vector>
#include <fstream>
#include <chrono>
#include <sstream>
#include <iterator>
#include <chrono>
using namespace std;
unordered_map<string,int> mymax(unordered_map<string,int>  l,unordered_map<string,int>  r) {
for(auto i=r.begin();i!=r.end();i++)
	l[i->first]+=i->second;
return l;
}


template<typename Out>
void split(const std::string &s, char delim, Out result) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
//std::transform(item.begin(), item.end(), item.begin(), ::tolower);
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
	cout<<"enter file name and number of processors"<<endl;
	string fname = argv[1];
	int p = atoi(argv[2]);
vector<string>words =readFile(fname);
cout<<words.size();
unordered_map<string,int> m,map;
 m = map;
 auto start = std::chrono::system_clock::now();

#pragma omp declare reduction \
  (rwz:unordered_map<string,int> :omp_out=mymax(omp_out,omp_in)) \
//  initializer(omp_priv=m)
#pragma omp parallel for reduction(rwz:m) //num_threads(24*p)
  for (int idata=0; idata<words.size(); idata++)
{
m[words[idata]] += 1;    
m = mymax(m,map);
}

  auto end = std::chrono::system_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
  std::cout << "nano seconds = "<<elapsed.count();
  auto nns = elapsed.count();
  elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
  std::cout << ","<<elapsed.count();
  ofstream myfile ("time.txt",ios::app);
  myfile<<"OMP"<<","<<fname<<","<<p<<","<<nns<<endl;
  myfile.close();

cout<<m["a"]<<endl;
cout<<m["b"]<<endl;
return 0;
}
