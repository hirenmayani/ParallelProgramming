/*
 *  * ParMatMul.c
 *   *
 *    *  Created on: Mar 12, 2018
 *     *      Author: hmayani@cs.stonybrook.edu
 *      */


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

unsigned int nthreads;

struct task_queue
{
    template < typename CALLABLE, typename... ARGS >
    void push_back( CALLABLE fn, ARGS&&... args )
    {
    		pending.push_back( std::bind( fn, args... ) ) ;
    }

    template < typename CALLABLE, typename... ARGS >
    void push_front( CALLABLE fn, ARGS&&... args )
    {
    		pending.push_front( std::bind( fn, args... ) ) ;
    }

/*    void execute_all_front(task_queue tqs[])
    	{
		cout << tid;
        while( !pending.empty())
        {
            pending.front()() ;
            pending.pop_front();
        }
		cout << tid << "done \n";
    	}
*/

	void execute_all_back(task_queue tqs[] )
		{
		cout << "---s--- \r\n";
		set <int, greater <int> > emptyQs;
		while(true)
		{
			
			if (!pending.empty()){
				pending.back()();
				pending.pop_back();
			}
			else{
				break;
				cout << "trying to steal \r\n";
				int victim = rand()%nthreads;


				if (!tqs[victim].pending.empty()){
					tqs[victim].pending.front()();
					tqs[victim].pending.pop_front();
				}else{
					emptyQs.insert(victim);
					cout << victim;
					if (emptyQs.size() == nthreads)
					{
						break;

				}

			}
		}
		cout << "---e---\r\n";
		}
	}

    std::deque< std::function< void() > > pending ;
    std::deque< std::function< void() > > p1 ;
    int tid;
};

#include <cilk/cilk.h>
task_queue *tqs;


void pMat(int** mat,int size)
{
  int i,j;


  for(i=0;i<size;i++)
  {
    for(j=0;j<size;j++)
      cout << mat[i][j]<< " ";
    cout << "\r\n";
    }
}

void indexMM(int zi, int zj, int xi, int xj, int yi, int yj, int currSize, int** Z, int** X, int** Y)
{
    int i,j,k;
    int n = pow(2,currSize);
    for(i=0;i<n;i++)
    		for(j=0;j<n;j++)
    			for(k=0;k<n;k++)
    				Z[zi+i][zj+j] += X[xi+i][xj+k]*Y[yi+k][yj+j];

}

void parRecMM_SC(int zi, int zj, int xi, int xj, int yi, int yj, int currSize, int** Z, int** X, int** Y, int minSize)
{
	if (currSize >= 1){
		cout << currSize<< "\r\n";

		tqs[rand()%nthreads].push_back(parRecMM_SC,zi, zj, xi, xj, yi, yj, currSize -1, Z, X, Y, minSize);
		parRecMM_SC(zi, zj, xi, xj, yi, yj, currSize -1, Z, X, Y, minSize);
	

	}
}

void parRecMM_SC1(int zi, int zj, int xi, int xj, int yi, int yj, int currSize, int** Z, int** X, int** Y, int minSize)
{
        if(currSize == 1){

        	indexMM(zi, zj, xi, xj, yi, yj, currSize, Z, X,  Y);
                pMat(Z,4);
		//cout << "executing small\r\n";
        }
        else{
            int n = pow(2,currSize-1);
            tqs[rand()%nthreads].push_back(parRecMM_SC,zi, zj, xi, xj, yi, yj, currSize -1, Z, X, Y, minSize);
            tqs[rand()%nthreads].push_back(parRecMM_SC,zi, zj+n, xi, xj, yi, yj+n, currSize - 1, Z, X, Y, minSize);
            tqs[rand()%nthreads].push_back(parRecMM_SC,zi+n, zj, xi+n, xj, yi, yj, currSize - 1, Z, X, Y, minSize);
			parRecMM_SC(zi+n, zj+n, xi+n, xj, yi, yj+n, currSize - 1, Z, X, Y, minSize);

		tqs[rand()%nthreads].push_back(parRecMM_SC,zi, zj, xi, xj+n, yi+n, yj, currSize - 1, Z, X, Y, minSize);
		tqs[rand()%nthreads].push_back(parRecMM_SC,zi, zj+n, xi, xj+n, yi+n, yj+n, currSize - 1, Z, X, Y, minSize);
		tqs[rand()%nthreads].push_back(parRecMM_SC,zi+n, zj, xi+n, xj+n, yi+n, yj, currSize - 1, Z, X, Y, minSize);
		parRecMM_SC(zi+n, zj+n, xi+n, xj+n, yi+n, yj+n, currSize - 1, Z, X, Y, minSize);

        }
}


void scheduler( int r, int** Z, int** X, int** Y, int minSize){

	//nthreads = std::thread::hardware_concurrency();
	nthreads = 1;
	printf("total processor : %d \n", nthreads);

	tqs = new task_queue[nthreads];

    for(int i = 0; i < nthreads; i++){
    		tqs[i].tid = i;
    }
	tqs[0].push_back(parRecMM_SC, 0,0,0,0,0,0, r, Z, X, Y, 2 ) ;

    for(int i = 0; i < nthreads; i++){
    		cilk_spawn tqs[i].execute_all_back(tqs) ;
    }

    cilk_sync;


}



void parRecMM(int zi, int zj, int xi, int xj, int yi, int yj, int currSize, int** Z, int** X, int** Y, int minSize)
{
        if(currSize == minSize){
        		indexMM(zi, zj, xi, xj, yi, yj, currSize, Z, X,  Y);
        }
        else{
            int n = pow(2,currSize-1);

			parRecMM(zi, zj, xi, xj, yi, yj, currSize -1, Z, X, Y, minSize);
			parRecMM(zi, zj+n, xi, xj, yi, yj+n, currSize - 1, Z, X, Y, minSize);
			parRecMM(zi+n, zj, xi+n, xj, yi, yj, currSize - 1, Z, X, Y, minSize);
			parRecMM(zi+n, zj+n, xi+n, xj, yi, yj+n, currSize - 1, Z, X, Y, minSize);

			parRecMM(zi, zj, xi, xj+n, yi+n, yj, currSize - 1, Z, X, Y, minSize);
			parRecMM(zi, zj+n, xi, xj+n, yi+n, yj+n, currSize - 1, Z, X, Y, minSize);
			parRecMM(zi+n, zj, xi+n, xj+n, yi+n, yj, currSize - 1, Z, X, Y, minSize);
			parRecMM(zi+n, zj+n, xi+n, xj+n, yi+n, yj+n, currSize - 1, Z, X, Y, minSize);


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


void printMat(int** mat,int size)
{
  int i,j;


  for(i=0;i<size;i++)
  {
    for(j=0;j<size;j++)
      printf("%d ",mat[i][j]);
    printf("\n");
    }
}


void matMulijk(int** X,int** Y,int** Z,int n)
{
      int i,j,k;

      for(i=0;i<n;i++)
      for(j=0;j<n;j++)
      for(k=0;k<n;k++)
        Z[i][j] = Z[i][j]+X[i][k]*Y[k][j];

}

int main(int argc,char* argv[])
{

  int r = atoi("2");
  int n = pow(2,r);

  int **X;
  int **Y;
  int **Zp;
  int **Zs;

  X = createMatrix(n,1);
  Y = createMatrix(n,1);

  Zs = createMatrix(n,0);
  Zp = createMatrix(n,0);

  X = createMatrix(n,1);
  Y = createMatrix(n,1);
  Zs = createMatrix(n,0);
  printf("\n-------------Serial---------------\n");
  matMulijk(X,Y,Zs,n);
  printMat(Zs,n);

  printf("\n-------------schedule Parallel---------------\n");
  scheduler(r, Zp, X, Y, 2);
  printMat(Zp,n);

  return 0;
}

int main1(int argc,char* argv[])
{

  int i=0;
  int j=0;
  int k=0;
  int r = atoi(argv[1]);
  int n = pow(2,r);

  int **X;
  int **Y;
  int **Zs;
  int **Zp;

  X = createMatrix(n,1);
  Y = createMatrix(n,1);
  Zs = createMatrix(n,0);
  printf("\n-------------Serial---------------\n");
  matMulijk(X,Y,Zs,n);
  printMat(Zs,n);

  Zp = createMatrix(n,0);
  printf("\n-------------Parallel---------------\n");
  parRecMM(0,0,0,0,0,0, r, Zp, X, Y, 2);
  printMat(Zp,n);

  return 0;
}


