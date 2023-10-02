#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <chrono>
#include <math.h>

std::vector<double> Ax_row( std::vector<double> &A, std::vector<double> &x ){

  size_t N = x.size();
  std::vector<double> y(N);

  //TOTO: Question 2a: Matrix vector multiplication with a row major matrix

  int block = 0;
  int rowEl = 0;
  int vecRow;
  for (int row=0; row < N; row++) {
    vecRow = 0;
    for (rowEl=block; rowEl < N+block; rowEl++) {
      y[row] += A[rowEl] * x[vecRow];
      vecRow++;
    }
    block += N;
  }
  
  // printf("This is y:\n");
  // for (int i=0; i<N; i++) {
  //   printf("y = %.1f\n", y[i]);
  // }

  return y;
}



std::vector<double> Ax_col( std::vector<double> &A, std::vector<double> &x ){
  size_t N = x.size();
  std::vector<double> y(N);

  //TOTO: Question 2a: Matrix vector multiplication with a column major matrix

  int vecRow;
  for (int row = 0; row < N; row++) {
    vecRow = 0;
    for (int rowEl = row; rowEl < N*N; rowEl+=N) {
      y[row] += A[rowEl] * x[vecRow];
      vecRow++;
    }
  }

  // printf("Column major\n");
  // printf("This is y:\n");
  // for (int i=0; i<N; i++) {
  //   printf("y = %.1f\n", y[i]);
  // }

  return y;
}


double benchmark_Ax( std::vector<double> &A, std::vector<double> &x, bool row_major, double Ns){

  double times = 0;

  for( size_t i=0; i<Ns; i++){
    auto t1 = std::chrono::system_clock::now();
    //TOTO: Question 2a: Call the function to be benchmarked
    if( row_major==true ){
      Ax_row(A, x);
    }
    else{
      Ax_col(A, x);
    }
    auto t2 = std::chrono::system_clock::now();
    times += std::chrono::duration<double>(t2-t1).count();
  }
  printf("Done in total %9.4fs  --  average %9.4fs\n", times, times/Ns);

  return times/Ns;
}


int main( int argc, char **argv )
{
  if( argc<3 ){
    printf("Usage: %s [N|matrix dimension] [Ns|number of iterations]\n",argv[0]);
    exit(0);
  }

  size_t N  = atoi(argv[1]) ;
  size_t Ns = atoi(argv[2]) ;
  std::vector<double> A(N*N), B(N*N), x(N);

  // TODO: Question 2a: Initialize matrices and vector
  //       store A as row major and B as column major
  
  // Store the row major order matrix. 
  // Store the col major order matrix, as the elements hold
  // the same values for a square matrix with the initialization of
  // A_i,j = i + j; x_i = i
  int vIdx = 0; // Variable to access indexes of matrices.
  for(int row=0; row < N; row++) {
    for (int col=0; col < N; col++) {
      A.at(vIdx) = row+col;
      B.at(vIdx) = row+col;
      vIdx++;
    }
  }

  // Initiate vector x. 
  for (int i=0; i < N; i++) {
    x.at(i) = i;
  }

  // // Print row major order matrix.
  // printf("Values of A\n"); 
  // for (int i=0; i < N*N; i++) {
  //   printf("%.1f\n", A[i]);
  // }

  // printf("Values of B\n");
  // // Print col major order matrix. 
  // for (int i=0; i < N*N; i++) {
  //   printf("%.1f\n", B[i]);
  // }


  printf("Working with matrix of dimension %zu\n",N);

  printf("A*x (row major).\n");
  double times1 = benchmark_Ax(A,x,true,Ns);

  printf("A*x (column major).\n");
  double times2 = benchmark_Ax(B,x,false,Ns);

  printf("-----------------\n");
  printf("Speedup %.8fs\n", times1/times2);


  return 0;
}
