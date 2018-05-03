/*
 * ParMatMul.c
 *
 *  Created on: Mar 12, 2018
 *      Author: hmayani@cs.stonybrook.edu
 */

#include <mpi.h>
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
#include <iostream>
#include <unistd.h>
#include<cilk/cilk.h>

using namespace std;

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

int** createContMatrix(int size, int init) {
	int rows = size;
	int cols = size;
	int *data = (int *) malloc(rows * cols * sizeof(int));
	int **array = (int **) malloc(rows * sizeof(int*));
	int i = 0, j = 0;
	for (i = 0; i < rows; i++)
		array[i] = &(data[cols * i]);

	if (init == 0) {
		for (i = 0; i < size; i++) {
			for (j = 0; j < size; j++) {
				array[i][j] = 0;

			}
		}
	} else if (init == -1) {
		for (i = 0; i < size; i++) {
			for (j = 0; j < size; j++) {
				array[i][j] = rand() % 30;

			}
		}
	} else {
		for (i = 0; i < size; i++) {
			for (j = 0; j < size; j++) {
				array[i][j] = init;

			}
		}
	}
	return array;

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

void matmul(int** Z, int** X, int** Y, int n) {
	/*
	 #pragma cilk grainsize = 5
	 cilk_for(unsigned int i = 0; i < n; ++i){
	 for (unsigned int k = 0; k < n; ++k) {
	 for (unsigned int j = 0; j < n; ++j) {
	 //  Z[i*n + j] += X[i*n + k] * Y[k*n + j];
	 Z[i][j] = Z[i][j] + X[i][k] * Y[k][j];
	 }
	 }
	 }*/

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++)
			for (int k = 0; k < n; k++)
				Z[i][j] = Z[i][j] + X[i][k] * Y[k][j];
	}

}

int main(int argc, char* argv[]) {

	int myrank, n = 0, p = 4;

	int r = atoi("3");
	n = pow(2, r);
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	int rootp = sqrt(p * 1.0);
	MPI_Request sendreq[rootp], recvreq[rootp];

	int i = floor(myrank / rootp);
	int j = myrank % rootp;
	int nbrp = n / rootp;
	int** A;
	int** B;
	int** C;

	A = createContMatrix(nbrp, 0);
	B = createContMatrix(nbrp, 0);
	C = createContMatrix(nbrp, 0);
	Cf = createContMatrix(nbrp, 0);

	/*
	 * divide A and B and send to all proc
	 */
	int ii, jj;
	int AA[n * n];
	int BB[n * n];
	int CC[n * n];
	int Ct[nbrp * nbrp];

	if (myrank == 0) {
		for (ii = 0; ii < n * n; ii++) {
			AA[ii] = ii;
			BB[ii] = ii;
		}
	}

	int a[nbrp * nbrp];
	int b[nbrp * nbrp];
	for (ii = 0; ii < nbrp * nbrp; ii++) {
		a[ii] = 0;
		b[ii] = 0;
	}

	MPI_Datatype blocktype;
	MPI_Datatype blocktype2;

	MPI_Type_vector(nbrp, nbrp, n, MPI_INT, &blocktype2);
	MPI_Type_create_resized(blocktype2, 0, sizeof(int), &blocktype);
	MPI_Type_commit(&blocktype);

	int disps[rootp * rootp];
	int counts[rootp * rootp];

	for (ii = 0; ii < rootp; ii++) {
		for (jj = 0; jj < rootp; jj++) {
			disps[ii * rootp + jj] = ii * n * nbrp + jj * nbrp;
			counts[ii * rootp + jj] = 1;
		}
	}

	MPI_Scatterv(AA, counts, disps, blocktype, a, nbrp * nbrp, MPI_INT, 0,
			MPI_COMM_WORLD);
	MPI_Scatterv(BB, counts, disps, blocktype, b, nbrp * nbrp, MPI_INT, 0,
			MPI_COMM_WORLD);

	for (ii = 0; ii < nbrp; ii++) {
		for (jj = 0; jj < nbrp; jj++) {
			A[ii][jj] = a[ii * nbrp + jj];
			B[ii][jj] = b[ii * nbrp + jj];
		}
		MPI_Barrier (MPI_COMM_WORLD);
	}
	/*
	 * A and B created
	 */

	//printMat(A,nbrp);
	//printMat(B, nbrp);
	int left = (rootp + j - i) % rootp;
	int up = (rootp - j + i) % rootp;
	int destA = i * rootp + left;
	int destB = up * rootp + j;

	int right = (rootp + j + i) % rootp;
	int down = (rootp + j + i) % rootp;
	int srcA = i * rootp + right;
	int srcB = down * rootp + j;

	MPI_Status sstatus[rootp + 1];
	MPI_Status rstatus[rootp + 1];

	MPI_Sendrecv_replace(&(A[0][0]), nbrp * nbrp, MPI_INT, destA, 123, srcA,
			123, MPI_COMM_WORLD, &sstatus[0]);
	MPI_Sendrecv_replace(&(B[0][0]), nbrp * nbrp, MPI_INT, destB, 23, srcB, 23,
			MPI_COMM_WORLD, &rstatus[0]);

	for (int l = 0; l < rootp; l++) {
		int k = (j + i + l - 1) % rootp;

		left = (rootp + j - 1) % rootp;
		up = (rootp + i - 1) % rootp;
		destA = i * rootp + left;
		destB = up * rootp + j;

		right = (rootp + j + 1) % rootp;
		down = (rootp + i + 1) % rootp;
		srcA = i * rootp + right;
		srcB = down * rootp + j;

		matmul(C, A, B, nbrp);
		//	printMat(C, nbrp);
		if (l < rootp) {
			MPI_Sendrecv_replace(&(A[0][0]), nbrp * nbrp, MPI_INT, destA, l,
					srcA, l, MPI_COMM_WORLD, &sstatus[l + 1]);
		}
		if (l < rootp) {
			MPI_Sendrecv_replace(&(B[0][0]), nbrp * nbrp, MPI_INT, destB,
					rootp + l, srcB, rootp + l, MPI_COMM_WORLD,
					&rstatus[l + 1]);
		}
	}
	for (ii = 0; ii < nbrp; ii++) {
		for (jj = 0; jj < nbrp; jj++) {
			Ct[ii * nbrp + jj] =C[ii][jj];
		}
		MPI_Barrier (MPI_COMM_WORLD);
}
	MPI_Gatherv(Ct, nbrp * nbrp, MPI_INT, CC, counts, disps, blocktype, 0,
			MPI_COMM_WORLD);

	if (myrank == 0) {

		for (ii = 0; ii < n * n; ii++) {
			/*if (ii % n == 0) {
				printf("\n");
			}
			printf("%d ", CC[ii]);*/
			Cf[int(ii/n)][int(ii%n)] = CC[ii];
		}
			printMat(Cf,n);

	}
//	printMat(C, nbrp);
	MPI_Finalize();

	return 0;
}

int broadcastAbroadcastB(int argc, char* argv[]) {

	int myrank, n = 0, p = 4;
	int r = atoi("10");
	n = pow(2, r);

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	int rootp = sqrt(p * 1.0);
	MPI_Status sendreq[rootp], recvreq[rootp];

	int i = floor(myrank / rootp);
	int j = myrank % rootp;
	int nbrp = n / rootp;
	int** A;
	int** B;
	int** C;
	int** At;
	int** Bt;
	A = createContMatrix(nbrp, myrank);
	B = createContMatrix(nbrp, myrank);
	C = createContMatrix(nbrp, 0);
	Bt = createContMatrix(nbrp, 0);
	At = createContMatrix(nbrp, 0);

//	cout << myrank << " left: " << left << " right:" << right << ": sending to: " << destA << "  receive from:" << srcA << "\n";

	MPI_Status sstatus[rootp + 1];

	int color = myrank / rootp;
	MPI_Comm row_comm;
	MPI_Comm_split(MPI_COMM_WORLD, color, myrank, &row_comm);

	int ccolor = myrank % rootp;
	MPI_Comm col_comm;
	MPI_Comm_split(MPI_COMM_WORLD, ccolor, myrank, &col_comm);

	int row_rank, row_size;
	MPI_Comm_rank(row_comm, &row_rank);
	MPI_Comm_size(row_comm, &row_size);

	int col_rank, col_size;
	MPI_Comm_rank(col_comm, &col_rank);
	MPI_Comm_size(col_comm, &col_size);
//	printMat(A,nbrp);
//	printf("WORLD RANK/SIZE: %d/%d \t ROW RANK/SIZE: %d/%d\n", myrank, p, row_rank, row_size);
	for (int l = 1; l <= rootp; l++) {
		int k = l - 1;

		if (k == i) {
			// TODO
			memcpy(&(Bt[0][0]), &(B[0][0]), nbrp * nbrp * sizeof(int));
		}
		if (k == j) {
			// TODO
			memcpy(&(At[0][0]), &(A[0][0]), nbrp * nbrp * sizeof(int));
		}
		//printMat(Bt,nbrp);

		int rootB = (l - 1) % rootp;	//TODO proper root
		MPI_Bcast(&(Bt[0][0]), nbrp * nbrp, MPI_INT, rootB, col_comm);
		MPI_Barrier(col_comm);

		int rootA = (l - 1) % rootp;	//TODO proper root
		MPI_Bcast(&(At[0][0]), nbrp * nbrp, MPI_INT, rootA, row_comm);
		MPI_Barrier(row_comm);
		//printf("%d %d %d \n", myrank, At[0][0], Bt[0][0]);
		matmul(C, At, Bt, nbrp);
	}

	//printMat(C, nbrp);
	//cout << myrank << "\n";
	MPI_Comm_free(&row_comm);
	MPI_Comm_free(&col_comm);
	MPI_Finalize();

	return 0;
}

int rotateAbroadcastB(int argc, char* argv[]) {

	int myrank, n = 0, p = 4;
	//int p=4;
	int r = atoi("3");
	n = pow(2, r);
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	int rootp = sqrt(p * 1.0);
	MPI_Status sendreq[rootp], recvreq[rootp];

	int i = floor(myrank / rootp);
	int j = myrank % rootp;
	int nbrp = n / rootp;
	int** A;
	int** B;
	int** C;
	A = createContMatrix(nbrp, myrank);
	B = createContMatrix(nbrp, myrank);
	C = createContMatrix(nbrp, 0);
	int** Bt;
	Bt = createContMatrix(nbrp, 0);

//	cout << myrank << " left: " << left << " right:" << right << ": sending to: " << destA << "  receive from:" << srcA << "\n";

	MPI_Status sstatus[rootp + 1];

	int color = myrank / rootp;
	MPI_Comm row_comm;
	MPI_Comm_split(MPI_COMM_WORLD, color, myrank, &row_comm);

	int ccolor = myrank % rootp;
	MPI_Comm col_comm;
	MPI_Comm_split(MPI_COMM_WORLD, ccolor, myrank, &col_comm);

	int row_rank, row_size;
	MPI_Comm_rank(row_comm, &row_rank);
	MPI_Comm_size(row_comm, &row_size);

	int col_rank, col_size;
	MPI_Comm_rank(col_comm, &col_rank);
	MPI_Comm_size(col_comm, &col_size);

//	printf("WORLD RANK/SIZE: %d/%d \t ROW RANK/SIZE: %d/%d\n", myrank, p, row_rank, row_size);
	printf("\n");
	for (int l = 1; l <= rootp; l++) {
		int k = (j + l - 1) % rootp;

		if (k == i) {
			// TODO
			memcpy(&(Bt[0][0]), &(B[0][0]), nbrp * nbrp * sizeof(int));
//printMat(Bt,nbrp);
//printMat(B,nbrp);
		}
		int root = (l + j - 1) % rootp;	//TODO proper root
		MPI_Bcast(&(Bt[0][0]), nbrp * nbrp, MPI_INT, root, col_comm);
		// Synchronize again before obtaining final time
		MPI_Barrier(col_comm);
//		}
//		else{
		//		printf("WORLD RANK/SIZE: %d/%d \t col RANK/SIZE: %d/%d\n", myrank, p, col_rank, col_size);

		//MPI_Bcast(&(B[0][0]), nbrp * nbrp, MPI_INT, , col_comm);
//		}
//		MPI_Barrier(col_comm);
		//printf("%d %d %d \n", myrank, A[0][0], Bt[0][0]);	
		matmul(C, A, Bt, nbrp);
		//cout <<"in loop" << myrank << "\n";

		int left = (rootp + j - 1) % rootp;
		int destA = i * rootp + left;
		int right = (rootp + j + 1) % rootp;
		int srcA = i * rootp + right;

//		cout << myrank << " "<< destA << " " << srcA << "\n"; 
		if (l < rootp) {
			//	if (myrank != destA)
			//                         MPI_Send(&(A[0][0]), nbrp*nbrp, MPI_INT, destA, l ,MPI_COMM_WORLD);
			//       if (myrank != srcA)
			//                       MPI_Recv(&(A[0][0]), nbrp*nbrp, MPI_INT, srcA, l,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Sendrecv_replace(&(A[0][0]), nbrp * nbrp, MPI_INT, destA, l,
					srcA, l, MPI_COMM_WORLD, &sstatus[l]);
		}
	}

	printMat(C, nbrp);
	//cout << myrank << "\n";
	MPI_Comm_free(&row_comm);
	MPI_Comm_free(&col_comm);
	MPI_Finalize();

	return 0;
}

int rotateBoth(int argc, char* argv[]) {

	int myrank, n = 0, p = 4;
	//int p=4;
	int r = atoi("1");
	n = pow(2, r);

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	int rootp = sqrt(p * 1.0);
	MPI_Request sendreq[rootp], recvreq[rootp];

	int i = floor(myrank / rootp);
	int j = myrank % rootp;
	int nbrp = n / rootp;

	int** A;
	int** B;
	int** C;
	A = createContMatrix(n, myrank);
	B = createContMatrix(n, myrank);
	C = createContMatrix(nbrp, 0);

	int left = (rootp + j - i) % rootp;
	int up = (rootp - j + i) % rootp;
	int destA = i * rootp + left;
	int destB = up * rootp + j;

	int right = (rootp + j + i) % rootp;
	int down = (rootp + j + i) % rootp;
	int srcA = i * rootp + right;
	int srcB = down * rootp + j;

//	cout << myrank << " left: " << left << " right:" << right << ": sending to: " << destA << "  receive from:" << srcA << "\n";

	MPI_Status sstatus[rootp + 1];
	MPI_Status rstatus[rootp + 1];

	MPI_Sendrecv_replace(&(A[0][0]), nbrp * nbrp, MPI_INT, destA, 123, srcA,
			123, MPI_COMM_WORLD, &sstatus[0]);
	MPI_Sendrecv_replace(&(B[0][0]), nbrp * nbrp, MPI_INT, destB, 23, srcB, 23,
			MPI_COMM_WORLD, &rstatus[0]);

	int** Aik;
	int** Bkj;

	for (int l = 0; l < rootp; l++) {
		int k = (j + i + l - 1) % rootp;

		left = (rootp + j - 1) % rootp;
		up = (rootp + i - 1) % rootp;
		destA = i * rootp + left;
		destB = up * rootp + j;

		right = (rootp + j + 1) % rootp;
		down = (rootp + i + 1) % rootp;
		srcA = i * rootp + right;
		srcB = down * rootp + j;

		matmul(C, A, B, nbrp);

		if (l < rootp) {
			MPI_Sendrecv_replace(&(A[0][0]), nbrp * nbrp, MPI_INT, destA, l,
					srcA, l, MPI_COMM_WORLD, &sstatus[l + 1]);
		}
		if (l < rootp) {
			MPI_Sendrecv_replace(&(B[0][0]), nbrp * nbrp, MPI_INT, destB,
					rootp + l, srcB, rootp + l, MPI_COMM_WORLD,
					&rstatus[l + 1]);
		}
	}

	printMat(C, nbrp);
	MPI_Finalize();

	return 0;
}

int kachara(int argc, char* argv[]) {

	int myrank, n = 0, p = 4;
	//int p=4;
	int r = atoi("3");
	n = pow(2, r);
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	int rootp = sqrt(p * 1.0);
	MPI_Request sendreq[rootp], recvreq[rootp];

//	int coords[2]; /* 2 Dimension topology so 2 coordinates */
//	MPI_Cart_coords(cart_comm,rank,2,coords);

	int i = floor(myrank / rootp);
	int j = myrank % rootp;

	int nbrp = n / rootp;

//	cout << myrank << ": rootp " << rootp << " i: " << i << " j: " << j << " nbrt:" << nbrp << "\n";
	int** A;
	int** B;
	int** C;
	A = createContMatrix(nbrp, myrank);
	B = createContMatrix(nbrp, myrank);
	C = createContMatrix(nbrp, 0);

	int left = (rootp + j - i) % rootp;
	int up = (rootp - j + i) % rootp;
	int destA = i * rootp + left;
	int destB = up * rootp + j;

	int right = (rootp + j + i) % rootp;
	int down = (rootp + j + i) % rootp;
	int srcA = i * rootp + right;
	int srcB = down * rootp + j;

//	cout << myrank << " left: " << left << " right:" << right << ": sending to: " << destA << "  receive from:" << srcA << "\n";

	MPI_Status sstatus[rootp + 1];
	MPI_Status rstatus[rootp + 1];

	//int q[2];
	//q[0]= myrank;
	MPI_Sendrecv_replace(&(A[0][0]), nbrp * nbrp, MPI_INT, destA, 123, srcA,
			123, MPI_COMM_WORLD, &sstatus[0]);
	MPI_Sendrecv_replace(&(B[0][0]), nbrp * nbrp, MPI_INT, destB, 23, srcB, 23,
			MPI_COMM_WORLD, &rstatus[0]);
	/*
	 if (myrank != destA)
	 MPI_Send(A, nbrp*nbrp, MPI_INT, destA, 0 ,MPI_COMM_WORLD);
	 if (myrank != srcA)
	 MPI_Recv(A, nbrp*nbrp, MPI_INT, srcA, 0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	 if (myrank != destB)
	 MPI_Send(B, nbrp*nbrp, MPI_INT, destB, 0 ,MPI_COMM_WORLD);

	 if (myrank != srcB)
	 MPI_Recv(B, nbrp*nbrp, MPI_INT, srcB, 0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	 */
//cout << myrank << " " << A[0][0] << "\n ";
//printMat(A,nbrp);
	int** Aik;
	int** Bkj;

	for (int l = 0; l < rootp; l++) {
		int k = (j + i + l - 1) % rootp;

		left = (rootp + j - 1) % rootp;
		up = (rootp + i - 1) % rootp;
		destA = i * rootp + left;
		destB = up * rootp + j;

		right = (rootp + j + 1) % rootp;
		down = (rootp + i + 1) % rootp;
		srcA = i * rootp + right;
		srcB = down * rootp + j;

//		printMat(A,nbrp);

		matmul(C, A, B, nbrp);

		if (l < rootp) {
			MPI_Sendrecv_replace(&(A[0][0]), nbrp * nbrp, MPI_INT, destA, l,
					srcA, l, MPI_COMM_WORLD, &sstatus[l + 1]);
			/*
			 send(Aik, i, left, myrank);
			 if (myrank != destA)
			 MPI_Send(A, nbrp*nbrp, MPI_INT, destA, l ,MPI_COMM_WORLD);
			 if (myrank != srcA)
			 MPI_Recv(A, nbrp*nbrp, MPI_INT, srcA, l,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			 */}
		if (l < rootp) {
			MPI_Sendrecv_replace(&(B[0][0]), nbrp * nbrp, MPI_INT, destB,
					rootp + l, srcB, rootp + l, MPI_COMM_WORLD,
					&rstatus[l + 1]);

			//	cout << myrank << " "<< destB << " "<< srcB << " "<< destA << " "<< srcA  << "\n";
			//send(Bkj, up, j, myrank);
			/*			if (myrank != destB)
			 MPI_Send(B, nbrp*nbrp, MPI_INT, destB, rootp+l ,MPI_COMM_WORLD);

			 if (myrank != srcB)
			 MPI_Recv(B, nbrp*nbrp, MPI_INT, srcB, rootp+l ,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			 */}
	}

	printMat(C, nbrp);

	/*
	 send(B,up,j, myrank);

	 int** Aik;
	 int** Bkj;

	 for (int l=0; l <rootp; l++){
	 int k = (j+i+l-1)%rootp;
	 int left = (rootp+j-1)%rootp;
	 int up = (rootp+i-1)%rootp;

	 matmul(C, Aik, Bkj);
	 if(l < rootp){
	 send(Aik, i, left, myrank);
	 }
	 if(l < rootp){
	 send(Bkj, up, j, myrank);
	 }
	 }
	 */
	/*
	 for (int i=0; i <rootp; i++){//Par
	 for (int j=0; j <rootp; j++){//Par
	 int left = (rootp+j-i)%rootp;
	 int up = (rootp-j+i)%rootp;

	 send(A,i,left, myrank);
	 send(B,up,j, myrank);
	 }
	 }

	 for (int l=0; l <rootp; l++){
	 for (int i=0; i <rootp; i++){ //Par
	 for (int j=0; j <rootp; j++){//Par
	 int k = (j+i+l-1)%rootp;
	 int left = (rootp+j-1)%rootp;
	 int up = (rootp+i-1)%rootp;

	 matmul(Cij, Aik, Bkj);
	 if(l < rootp){
	 send(Aik, i, left, myrank);
	 }
	 if(l < rootp){
	 send(Bkj, up, j, myrank);
	 }
	 }
	 }
	 }
	 */
	/*
	 if(myrank  == 0){
	 if((n/rootp)* (n/rootp)!=n) {
	 printf(" of Proc must be equal to %d\nCode terminated",r);
	 exit(0);
	 }
	 }
	 */

	MPI_Finalize();

	/*
	 printf("\n-------------Serial---------------\n");
	 matMulijk(X,Y,Zs,n);
	 printMat(Zs,n);

	 printf("\n-------------schedule Parallel---------------\n");
	 */
	// MM_rotateA_rotateB(A,B, C, n, p);
	//printMat(Zp,n);
	return 0;
}

