/*
 *  *  *  * ParMatMul.c
 *   *   *   *
 *    *    *    *  Created on: Mar 12, 2018
 *     *     *     *      Author: hmayani@cs.stonybrook.edu
 *      *      *      */

#include <iostream>
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
#include<fstream>
#include<papi.h>
#include <unistd.h>
#include<chrono>
using namespace std;

unsigned int nthreads = 64;
unsigned int freeT;
tbb::mutex l, l2;
time_t t;

struct mmArgs {
	int zi;
	int zj;
	int xi;
	int xj;
	int yi;
	int yj;
	int currSize = -1;
	int minSize;
	int** Z;
	int** X;
	int** Y;
};

void parRecMM_SC(mmArgs a, int stype, int myid);
int emptyQs = 0;
int *work;

struct task_queue {

	mmArgs cpop_back() {
		mmArgs mma;

		l.lock();
		if (!pending.empty()) {
			mma = pending.back();
			pending.pop_back();
		}
		l.unlock();

		if (mma.currSize > 0) {
			int n = mma.currSize;
			work[tid] -= n*n*n;
		}
		return mma;
	}

	mmArgs cpop_front() {
		mmArgs mma;

		l.lock();
		if (!pending.empty()) {
			mma = pending.front();
			pending.pop_front();
		}
		l.unlock();

		if (mma.currSize > 0) {
			int n = mma.currSize;
			work[tid] -= n*n*n;
		}
		return mma;
	}

	mmArgs back() {
		return pending.back();
	}

	void pop_back() {
		l.lock();
		pending.pop_back();
		l.unlock();
	}

	mmArgs front() {
		return pending.front();
	}

	void pop_front() {
		l.lock();
		pending.pop_front();
		l.unlock();
	}

	void push_back(mmArgs mma) {
		l.lock();
		pending.push_back(mma);
		l.unlock();
		int n = mma.currSize;
		work[tid] += n*n*n;

	}

	void push_front(mmArgs mma) {
		l.lock();
		pending.push_front(mma);
		l.unlock();
		int n = mma.currSize;
		work[tid] += n*n*n;
	}

	bool empty() {
		return pending.empty();
	}

	std::deque<mmArgs> pending;
	int tid;
};

task_queue *tqs;
task_queue ddq;

void DRSteal(task_queue tqs[], int myid) {

	int atmps = 0;
	while (true) {
		mmArgs mma = tqs[myid].cpop_back();
		if (mma.currSize > 0) {
			atmps = 0;
			parRecMM_SC(mma, 1, myid);
		} else {
			int victim = rand() % nthreads;
			mmArgs mma = tqs[victim].cpop_front();
			if (mma.currSize > 0) {
				atmps = 0;
				parRecMM_SC(mma, 1, myid);
			} else {
				atmps += 1;
				if (atmps == 1000) {
					break;
				}
			}

		}
	}
}

void DRShare(task_queue tqs[], int myid) {
	unsigned int emptyc = 0;
	while (true) {
		mmArgs mma = tqs[myid].cpop_back();
		if (mma.currSize > 0) {
			emptyc = 0;
			parRecMM_SC(mma, 2, myid);
		} else {
			emptyc += 1;
			if (emptyc == 9999) {
				break;
			}
		}
	}
}

void CShare() {
        bool prevFree =false;
        unsigned int emptyc = 0;

        while (true) {

                mmArgs mma = ddq.cpop_back();
                if (mma.currSize <= 0) {
                        emptyc += 1;
                        if (emptyc == 9999) {
                                break;
                        }

                }else{
                        parRecMM_SC(mma, 3, -1);
                }
        }
}


void CShare1() {
        bool prevFree =false;

        while (true) {

                if(prevFree){
                        prevFree = false;
                        l2.lock();
                                emptyQs-=1;
                        l2.unlock();
                }
                mmArgs mma = ddq.cpop_back();
                if (mma.currSize <= 0) {
                        if(!prevFree){
                                prevFree = true;
                                l2.lock();
                                        emptyQs+=1;
                                l2.unlock();
                        }
                        if (emptyQs >= nthreads){
                                break;
                        }
                }else{
                        parRecMM_SC(mma, 3, -1);
                }
        }
}


void DRStealMod(task_queue tqs[], int myid) {

	int atmps = 0;
	while (true) {
		mmArgs mma = tqs[myid].cpop_back();
		if (mma.currSize > 0) {
			atmps = 0;
			parRecMM_SC(mma, 1, myid);
		} else {
			int victim1 = rand() % nthreads;
			int victim2 = rand() % nthreads;
			int victim;
			if (work[victim1]  > work[victim2]){
				victim = victim1;
			}else{
				victim = victim2;
			}
			mmArgs mma = tqs[victim].cpop_front();
			if (mma.currSize > 0) {
				atmps = 0;
				parRecMM_SC(mma, 1, myid);
			} else {
				atmps += 1;
				if (atmps == 1000) {
					break;
				}
			}

		}
	}
}

void DRShareMod(task_queue tqs[], int myid) {
	unsigned int emptyc = 0;
	while (true) {
		mmArgs mma = tqs[myid].cpop_back();
		if (mma.currSize > 0) {
			emptyc = 0;
			parRecMM_SC(mma, 2, myid);
		} else {
			emptyc += 1;
			if (emptyc == 9999) {
				break;
			}
		}
	}
}


void indexMM(int zi, int zj, int xi, int xj, int yi, int yj, int currSize,
		int** Z, int** X, int** Y) {
	int i, j, k;
	int n = pow(2, currSize);
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			for (k = 0; k < n; k++)
				Z[zi + i][zj + j] += X[xi + i][xj + k] * Y[yi + k][yj + j];

}

void parRecMM_SC(mmArgs mma, int stype, int myid) {

	if (mma.currSize == mma.minSize) {
		indexMM(mma.zi, mma.zj, mma.xi, mma.xj, mma.yi, mma.yj, mma.currSize,
				mma.Z, mma.X, mma.Y);
	} else {
		mma.currSize -= 1;
		int n = pow(2, mma.currSize);
		mmArgs mma1 = mma;

		mmArgs mma2 = mma;
		mma2.zj += n;
		mma2.yj += n;
		mmArgs mma3 = mma;
		mma3.zi += n;
		mma3.xi += n;

		mmArgs mma4 = mma;
		mma4.zi += n;
		mma4.zj += n;
		mma4.xi += n;
		mma4.yj += n;

		if (stype == 3) {
			ddq.push_front(mma1);
			ddq.push_front(mma2);
			ddq.push_front(mma3);
		} else if(stype == 2) {
			tqs[rand() % nthreads].push_front(mma1);
			tqs[rand() % nthreads].push_front(mma2);
			tqs[rand() % nthreads].push_front(mma3);
		} else if(stype == 1 or stype==4) {
			tqs[myid].push_front(mma1);
			tqs[myid].push_front(mma2);
			tqs[myid].push_front(mma3);
		} else if(stype == 5) {
			int victim1 = rand() % nthreads;
			int victim2 = rand() % nthreads;
			int victim;
			if (work[victim1]  < work[victim2]){
				victim = victim1;
			}else{
				victim = victim2;
			}

			tqs[victim].push_front(mma1);
			tqs[victim].push_front(mma2);
			tqs[victim].push_front(mma3);
		}
		parRecMM_SC(mma4, stype, myid);

		mmArgs mma5 = mma;
		mma5.xj += n;
		mma5.yi += n;

		mmArgs mma6 = mma;
		mma6.xj += n;
		mma6.yi += n;
		mma6.yj += n;
		mma6.zj += n;

		mmArgs mma7 = mma;
		mma7.xi += n;
		mma7.xj += n;
		mma7.yi += n;
		mma7.zi += n;

		if (stype == 3) {
			ddq.push_front(mma4);
			ddq.push_front(mma5);
			ddq.push_front(mma6);
		} else if(stype == 2) {
			tqs[rand() % nthreads].push_front(mma5);
			tqs[rand() % nthreads].push_front(mma6);
			tqs[rand() % nthreads].push_front(mma7);
		} else if(stype == 1 or stype == 4) {
			tqs[myid].push_front(mma5);
			tqs[myid].push_front(mma6);
			tqs[myid].push_front(mma7);
		} else if(stype == 5) {
			int victim1 = rand() % nthreads;
			int victim2 = rand() % nthreads;
			int victim;
			if (work[victim1]  < work[victim2]){
				victim = victim1;
			}else{
				victim = victim2;
			}

			tqs[victim].push_front(mma5);
			tqs[victim].push_front(mma6);
			tqs[victim].push_front(mma7);
		}
		mmArgs mma8 = mma;
		mma8.xi += n;
		mma8.xj += n;
		mma8.yi += n;
		mma8.yj += n;
		mma8.zi += n;
		mma8.zj += n;

		parRecMM_SC(mma8, stype, myid);

	}
	return;
}

void scheduler(int r, int** Z, int** X, int** Y, int minSize, int sType) {
	/*      nthreads = std::thread::hardware_concurrency();*/
	printf("total processor : %d \n", nthreads);

	tqs = new task_queue[nthreads];
	work = new int[nthreads];

	for (int i = 0; i < nthreads; i++) {
		tqs[i].tid = i;
	}

	mmArgs a;
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
	ddq.push_back(a);
	tqs[0].push_back(a);
	int n = pow(2, r);
	work[0] = n*n*n;

    for (int i = 0; i < nthreads; i++) {
            if (sType == 1)
                    cilk_spawn DRSteal(tqs, i);
            else if (sType == 2)
                    cilk_spawn DRShare(tqs,i) ;
            else if (sType == 3)
                    cilk_spawn CShare();
            else if (sType == 4)
                cilk_spawn DRStealMod(tqs, i);
            else if (sType == 5)
                cilk_spawn DRShareMod(tqs,i) ;

	}

	cilk_sync;

}

void parRecMM(int zi, int zj, int xi, int xj, int yi, int yj, int currSize,
		int** Z, int** X, int** Y, int minSize) {
	if (currSize == minSize) {
		indexMM(zi, zj, xi, xj, yi, yj, currSize, Z, X, Y);
	} else {
		int n = pow(2, currSize - 1);

		cilk_spawn parRecMM(zi, zj, xi, xj, yi, yj, currSize - 1, Z, X, Y,
				minSize);
		parRecMM(zi, zj + n, xi, xj, yi, yj + n, currSize - 1, Z, X, Y,
				minSize);
		parRecMM(zi + n, zj, xi + n, xj, yi, yj, currSize - 1, Z, X, Y,
				minSize);
		parRecMM(zi + n, zj + n, xi + n, xj, yi, yj + n, currSize - 1, Z, X, Y,
				minSize);
		cilk_sync;

		cilk_spawn parRecMM(zi, zj, xi, xj + n, yi + n, yj, currSize - 1, Z, X,
				Y, minSize);
		parRecMM(zi, zj + n, xi, xj + n, yi + n, yj + n, currSize - 1, Z, X, Y,
				minSize);
		parRecMM(zi + n, zj, xi + n, xj + n, yi + n, yj, currSize - 1, Z, X, Y,
				minSize);
		parRecMM(zi + n, zj + n, xi + n, xj + n, yi + n, yj + n, currSize - 1,
				Z, X, Y, minSize);
		cilk_sync;
	}
}

int** createMatrix(int size, int init) {
	int i = 0;
	int j = 0;

	int **mat;
	mat = (int**) calloc(sizeof(int*), size);
	for (int i = 0; i < size; i++) {
		mat[i] = (int*) calloc(sizeof(int), size);
	}

	if (init == 0) {
		for (i = 0; i < size; i++) {
			for (j = 0; j < size; j++) {
				mat[i][j] = 0;

			}
		}
	} else {
		for (i = 0; i < size; i++) {
			for (j = 0; j < size; j++) {
				mat[i][j] = rand() % 30;

			}
		}
	}
	return mat;

}

int** printMatMulTime(int** X,int** Y,int n,void(*matmul_fun_pointer)(int, int**, int**, int**, int, int),int r,int sType)

{

int retval, EventSet = PAPI_NULL;

long long values[5]={0,0,0,0,0};

    unsigned counter;

    unsigned c;

    unsigned long fact;

    unsigned stoppoint;





        /* Initialize the PAPI library */

        retval = PAPI_library_init(PAPI_VER_CURRENT);



        if (retval != PAPI_VER_CURRENT) {

          fprintf(stderr, "PAPI library init error!\n");

          exit(1);

        }



/* Create the Event Set */

        if (PAPI_create_eventset(&EventSet) != PAPI_OK)

            printf ("%s:%d\t ERROR\n", __FILE__, __LINE__);



        /* Add Total Instructions Executed to our EventSet */

        if (PAPI_add_event(EventSet, PAPI_L1_TCM) != PAPI_OK)

            printf ("%s:%d\t ERROR\n", __FILE__, __LINE__);



        /* Add Total Instructions Executed to our EventSet */

        if (PAPI_add_event(EventSet, PAPI_L2_TCM) != PAPI_OK)

            printf ("%s:%d\t ERROR\n", __FILE__, __LINE__);



        /* Add Total Instructions Executed to our EventSet */

        if (PAPI_add_event(EventSet, PAPI_L3_TCM) != PAPI_OK)

            printf ("%s:%d\t ERROR\n", __FILE__, __LINE__);

        int ret = PAPI_add_event(EventSet, PAPI_L1_DCM);

        if ( ret!= PAPI_OK)

                    printf ("%s:%d\t ERROR %d\n", __FILE__, __LINE__,ret);

        ret = PAPI_add_event(EventSet, PAPI_L2_DCM);

        if (ret != PAPI_OK)

                    printf ("%s:%d\t ERROR %d\n", __FILE__, __LINE__,ret);





   srand ( time(NULL) );

   stoppoint = 2;//50+(rand() % 100);

/* Do some computation here*/



     	/* Initialize the PAPI library */

        retval = PAPI_library_init(PAPI_VER_CURRENT);



        if (retval != PAPI_VER_CURRENT) {

          fprintf(stderr, "PAPI library init error!\n");

          exit(1);

        }





        /* Start counting */

        if (PAPI_start(EventSet) != PAPI_OK)

        	printf ("%s:%d\t ERROR\n", __FILE__, __LINE__);

//        struct timespec start, finish;



        /* ... */






	auto start = chrono::system_clock::now();
        int** Z = createMatrix(n,0);

 //       clock_gettime(CLOCK_REALTIME, &start);

        matmul_fun_pointer(r, Z, X, Y, 5, sType);
	auto end = chrono::system_clock::now();

      /*  clock_gettime(CLOCK_REALTIME, &finish);

        elapsed = (finish.tv_sec - start.tv_sec);

        elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
*/
        //elapsed = 2*n*n*n/elapsed;
        //elapsed = elapsed/1000000000;


	auto elapsedT = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
        auto elapsed = 2*n*n*n/elapsedT.count();

	if (PAPI_read(EventSet, values) != PAPI_OK)

                    printf ("%s:%d\t ERROR\n", __FILE__, __LINE__);

                printf ("\nL1_TCM:\t%lld\t L2_TCM:\t%lld\n\n", values[0], values[1]);

                /* Do some computation here */

 

                if (PAPI_stop(EventSet, values) != PAPI_OK)

                    printf ("%s:%d\t ERROR\n", __FILE__, __LINE__);



                ofstream myfile ("ans2_ab.csv",ios::app);
myfile<<sType<<","<< n<< ",";
myfile<<elapsedT.count()<<",";
	 myfile <<elapsed<<",";
myfile<<values[0]<<",";
myfile<<values[1]<<endl;
                    myfile.close();

/*

                ofstream myfile ("ans2_cd.csv",ios::app);
 myfile << n<<","<<elapsedT.count()<<endl;
                    myfile.close();

*/


                return Z;

}
void printMat(int** mat, int size) {
	int i, j;

	for (i = 0; i < size; i++) {
		for (j = 0; j < size; j++)
			printf("%d ", mat[i][j]);
		printf("\n");
	}
}

void matMulijk(int** X, int** Y, int** Z, int n) {
	int i, j, k;

	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			for (k = 0; k < n; k++)
				Z[i][j] = Z[i][j] + X[i][k] * Y[k][j];

}

bool areEql(int** m1, int** m2, int size) {
	int i, j;

	for (i = 0; i < size; i++) {
		for (j = 0; j < size; j++)
			if (m1[i][j] != m2[i][j]) {
				return false;
			}
	}
	return true;
}

int main(int argc, char* argv[]) {

	int r = atoi(argv[1]);
	int n = pow(2, r);
	int sType = atoi(argv[2]);
	
	int **X;
	int **Y;
	int **Zp;
	int **Zs;
	nthreads = atoi(argv[3]);
	srand((unsigned) time(&t));
	X = createMatrix(n, 1);
	srand((unsigned) time(&t));
	Y = createMatrix(n, 1);
	Zp = createMatrix(n, 0);
	Zs = createMatrix(n, 0);

	//printf("\n-------------Serial---------------\n");
	matMulijk(X, Y, Zs, n);
	//printMat(Zs,n);
	srand((unsigned) time(&t));

	//printf("\n-------------schedule Parallel---------------\n");
	//scheduler(r, Zp, X, Y, 3, sType);
	//printMat(Zp,n);

	printMatMulTime(X,Y,n,&scheduler,r,sType);

	if (areEql(Zs, Zp, n)) {
		cout << "valid\n";
	} else {
		cout << "invalid\n";
	}
	return 0;
}

int main1(int argc, char* argv[]) {

	int i = 0;
	int j = 0;
	int k = 0;
	int r = atoi(argv[1]);
	int n = pow(2, r);

	int **X;
	int **Y;
	int **Zs;
	int **Zp;

	X = createMatrix(n, 1);
	Y = createMatrix(n, 1);
	Zs = createMatrix(n, 0);
	printf("\n-------------Serial---------------\n");
	matMulijk(X, Y, Zs, n);
	printMat(Zs, n);

	Zp = createMatrix(n, 0);
	printf("\n-------------Parallel---------------\n");
	parRecMM(0, 0, 0, 0, 0, 0, r, Zp, X, Y, 2);
	printMat(Zp, n);

	return 0;
}


