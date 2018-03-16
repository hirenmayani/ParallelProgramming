/*
 * play.cpp
 *
 *  Created on: Mar 13, 2018
 *      Author: hiren
 */




#include <thread>
#include <iostream>
#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include <thread>
#include <iostream>
#include <chrono>
#include<math.h>
#include <stdlib.h>
#include<stdio.h>
#include<string.h>
#include <deque>
#include <functional>
#include <queue>
#include <set>
using namespace std;
using namespace std;

unsigned int nthreads;


struct as{
	int zi;
	int zj;
	int xi;
	int xj;
	int yi;
	int yj;
	int currSize;
	int** Z;
	int** X;
	int** Y;

};

void parRecMM_SC(as a);

struct task_queue
{
    void push_back( as a)
    {
    		cout << "pushing " ;
    		pending.push_back(a) ;
    }


	void execute_all_back(task_queue tqs[])
		{
		cout << tid << "started \n";
		while(true)
		{
			if (!pending.empty()){
				as a  = pending.back();
				pending.pop_back();
				parRecMM_SC(a);
			}
			else{
				break;
			}
		}
		cout << tid << "done \n";

		}

    std::deque<as> pending ;
    int tid;
};

task_queue *tqs;

void parRecMM_SC(as a)
	{
			as na = a;
			a.currSize -=1;
	        if (a.currSize >= 2){
	        		cout << "if\r\n";
	                tqs[0].push_back(a);
	                parRecMM_SC(a);


	        }
	        else{
	        	cout << "else\r\n";
	        }
	}

int** createMatrix(int size,int init)
{
  int i=0;
  int j=0;


  int **mat;
  mat = (int**)calloc(sizeof(int*),size);
  for(int i = 0; i < size; i++)
  {
      mat[i] = (int*)calloc(sizeof(int), size);
  }

  if(init == 0)
  {
    for(i=0;i<size;i++)
    {
    for(j=0;j<size;j++)
    {
          mat[i][j] = 0;

        }
      }
  }
  else
  {
    for(i=0;i<size;i++)
    {
    for(j=0;j<size;j++)
    {
          mat[i][j] = rand()%30;

        }
    }
  }
return mat;

}

int main ()
{
	unsigned int nthreads = 1;
	cout << nthreads;
	printf("%d ",nthreads);

	tqs = new task_queue[nthreads];

	    for(int i = 0; i < nthreads; i++){
	    		tqs[i].tid = i;
	    }
	    as a ;
		a.zi = 0;
		a.zj = 0;
		a.xi = 0;
		a.xj = 0;
		a.yi = 0;
		a.yj = 0;
		a.currSize = 4;
		int** Z;
		int** X;
		int** Y;
		  X = createMatrix(9,1);
		  Y = createMatrix(9,1);
		  Z = createMatrix(9,0);
		 a.X = X;
		 a.Y = Y;
		 a.Z = Z;

		tqs[0].push_back(a) ;

	    for(int i = 0; i < nthreads; i++){
	    		 tqs[i].execute_all_back(tqs) ;
	    }


}
