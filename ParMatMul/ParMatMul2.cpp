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
#include <cilk/cilk.h>
#include <tbb/mutex.h>

#include <unistd.h>
using namespace std;

unsigned int nthreads;
struct mmArgs{
        int zi;
        int zj;
        int xi;
        int xj;
        int yi;
        int yj;
        int currSize;
        int minSize;
        int** Z;
        int** X;
        int** Y;

};

void parRecMM_SC(mmArgs a);
std::deque< mmArgs > ddq;

struct task_queue
{
    void push_back( mmArgs mma)
    {
                pending.push_back( mma) ;
    }
    void push_front( mmArgs mma )
    {
                pending.push_front( mma) ;
    }

    void execute_all_front(task_queue tqs[])
        {
                cout <<" not definded \n";
        }


        void DRSteal(task_queue tqs[])
                {
                tbb::mutex m1,m2;

                set <int, greater <int> > emptyQs;
                while(true)
                {
                        if (!pending.empty()){
                        m1.lock();
                                        mmArgs mma = pending.back();
                                        pending.pop_back();
                        m1.unlock();
                                        parRecMM_SC(mma);
                        }
                        else{
                                int victim = rand()%nthreads;
                                if (!tqs[victim].pending.empty()){
                                        m2.lock();
                                                mmArgs mma = tqs[victim].pending.front();
                                                tqs[victim].pending.pop_front();
                                        m2.unlock();
                                                parRecMM_SC(mma);
                                }else{
                                        emptyQs.insert(victim);
                                        if (emptyQs.size() == nthreads)
                                        {
                                                break;
                                        }
                                }

                        }
                }
        }
        void DRShare(task_queue tqs[])
                    {
                    tbb::mutex m1,m2;
                    unsigned int emptyC = 0;
                    set <int, greater <int> > emptyQs;
                    while(true)
                    {

                            if (!pending.empty()){
                            emptyC = 0;
                                    m1.lock();
                                            mmArgs mma = pending.back();
                                            pending.pop_back();
                                    m1.unlock();
                                    parRecMM_SC(mma);
                            }
                            else{
                                    emptyC +=1;
                                    sleep(1);
                                    if (emptyC ==5){
                                            break;
                                    }
                            }
                    }
            }

            void CShare(){
                    while (!ddq.empty()){
                            mmArgs mma = ddq.front();
                            ddq.pop_front();
                            parRecMM_SC(mma);
                    }
            }

        std::deque< mmArgs > pending ;
        int tid;
    };

    task_queue *tqs;


    void indexMM(int zi, int zj, int xi, int xj, int yi, int yj, int currSize, int** Z, int** X, int** Y)
    {
        int i,j,k;
        int n = pow(2,currSize);
        for(i=0;i<n;i++)
                    for(j=0;j<n;j++)
                            for(k=0;k<n;k++)
                                    Z[zi+i][zj+j] += X[xi+i][xj+k]*Y[yi+k][yj+j];

    }


    void parRecMM_SC(mmArgs mma)
    {

            if(mma.currSize == mma.minSize){
                            indexMM(mma.zi, mma.zj, mma.xi, mma.xj, mma.yi, mma.yj, mma.currSize, mma.Z, mma.X, mma.Y);
            }
            else{
                            mma.currSize-=1;
                int n = pow(2,mma.currSize);
                mmArgs mma1 = mma;
                tqs[rand()%nthreads].push_front(mma1);

                mmArgs mma2 = mma;
                mma2.zj+=n;
                mma2.yj+=n;
                tqs[rand()%nthreads].push_front(mma2);

                mmArgs mma3 = mma;
                mma3.zi+=n;
                mma3.xi+=n;
                tqs[rand()%nthreads].push_front(mma3);


                mmArgs mma4 = mma;
                mma4.zi+=n;
                mma4.zj+=n;
                mma4.xi+=n;
                mma4.yj+=n;

                    ddq.push_front(mma1);
                    ddq.push_front(mma2);
                    ddq.push_front(mma3);

                    parRecMM_SC(mma4);

                    mmArgs mma5 = mma;
                    mma5.xj+=n;
                    mma5.yi+=n;
                        tqs[rand()%nthreads].push_front(mma5);


                    mmArgs mma6 = mma;
                    mma6.xj+=n;
                    mma6.yi+=n;
                    mma6.yj+=n;
                    mma6.zj+=n;
                                tqs[rand()%nthreads].push_front(mma6);


                    mmArgs mma7 = mma;
                    mma7.xi+=n;
                    mma7.xj+=n;
                    mma7.yi+=n;
                    mma7.zi+=n;
                                tqs[rand()%nthreads].push_front(mma7);

                    mmArgs mma8 = mma;
                    mma8.xi+=n;
                    mma8.xj+=n;
                    mma8.yi+=n;
                    mma8.yj+=n;
                    mma8.zi+=n;
                    mma8.zj+=n;

                        ddq.push_front(mma5);
                        ddq.push_front(mma6);
                        ddq.push_front(mma7);
                                parRecMM_SC(mma8);

                }
        }


        void scheduler( int r, int** Z, int** X, int** Y, int minSize){
            nthreads = std::thread::hardware_concurrency();
            printf("total processor : %d \n", nthreads);

            tqs = new task_queue[nthreads];

        for(int i = 0; i < nthreads; i++){
                    tqs[i].tid = i;
        }
        mmArgs a ;
            a.zi = 0;
            a.zj = 0;
            a.xi = 0;
            a.xj = 0;
            a.yi = 0;
            a.yj = 0;
            a.currSize = r;
            a.minSize = minSize;
             a.X = X;
             a.Y = Y;
             a.Z = Z;
            //ddq = {0};
            //ddq = tqs[0].pending;
            ddq.push_back(a);
        for(int i = 0; i < nthreads; i++){
                    cilk_spawn tqs[i].CShare() ;
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
                                                                                                                                                                    270,1-8       63%
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

      int r = atoi(argv[1]);
      int n = pow(2,r);

      int **X;
      int **Y;
      int **Zp;
      int **Zs;

      X = createMatrix(n,1);
      Y = createMatrix(n,1);
      Zp = createMatrix(n,0);
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
