
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


using namespace std;
unsigned int nthreads = 4;
unsigned int freeT;
tbb::mutex l1,l2,l3,l4,l5,l6;

struct mmArgs{
        int currSize;
};

struct task_queue
{

	mmArgs cpop_back()
    {
		mmArgs mma;

		l1.lock();
		if (!pending.empty())
		{
			mma = pending.back() ;
			pending.pop_back() ;
		}
		l1.unlock();

		return mma;
    }

	mmArgs cpop_front()
    {
		mmArgs mma;

		l2.lock();
			if (!pending.empty())
			{
				mma = pending.front() ;
				pending.pop_front() ;
			}
		l2.unlock();

		return mma;
    }

	mmArgs back()
    {
                return pending.back() ;
    }

	void pop_back()
    {
		l3.lock();
			pending.pop_back() ;
		l3.unlock();
    }


	mmArgs front()
    {
                return pending.front() ;
    }

	void pop_front()
    {
		l4.lock();
          pending.pop_front() ;
         l4.unlock();
    }

	void push_back( mmArgs mma)
    {
		l5.lock();
			pending.push_back( mma) ;
		l5.unlock();
    }

    void push_front( mmArgs mma )
    {
    		l6.lock();
    			pending.push_front( mma) ;
    		l6.unlock();
    }

    bool empty(){
    		return pending.empty();
    }

	std::deque< mmArgs > pending ;
	int tid;
 };

task_queue *tqs;



void parRecMM_SC(mmArgs mma)
{

		if(mma.currSize == 1){
			cout << mma.currSize;
		}
		else{
						mma.currSize-=1;
			int n = pow(2,mma.currSize);

			mmArgs mma3 = mma;
			tqs[rand()%nthreads].push_front(mma3);

			mmArgs mma8 = mma;
			parRecMM_SC(mma8);

			}
	}
void DRSteal(task_queue tqs[], int myid)
        {

        int atmps = 0;
        while(true)
        {
			mmArgs mma = tqs[myid].cpop_back();
			if (mma.currSize !=  0){
				atmps = 0;
				parRecMM_SC(mma);
			}
			else{
					int victim = rand()%nthreads;
					mmArgs mma = tqs[victim].cpop_front();
					if (mma.currSize == 0){
						atmps = 0;
						parRecMM_SC(mma);
					}else{
						atmps +=1;
						if (atmps == 1000){
							break;
						}
					}

                }
        }
}

int main()
{
	nthreads=4;
	tqs = new task_queue[nthreads];
	mmArgs a ;
	a.currSize = 5;
	tqs[0].push_back(a);
	for(int i = 0; i < nthreads; i++){
		cilk_spawn DRSteal(tqs,i) ;
		/*cilk_spawn CShare() ;
 * 		/*cilk_spawn DRShare(tqs,i) ;
 * 				/cilk_spawn DRSteal(tqs, i) ;*/
	}

cilk_sync;
}

