#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <chrono>
#include <math.h>
#include <string>


std::vector<double> transpose( std::vector<double> A ){
  
  size_t N = sqrt(A.size());
  std::vector<double> AT(N*N);

  // Matrix stored in row major order, loops over the rows. 
  for(int col = 0; col < N; col++) {
    for(int row = 0; row < N; row++) {
      AT.at(col+row*N) = A.at(row+col*N);
    }
  }
  return AT;
}

void AB( std::vector<double> A, std::vector<double> B ){

  size_t N = sqrt(A.size());
  std::vector<double> C(N*N);

  // TODO: Question 2c: Straightforward matrix-matrix multiplication
  
  double r;
  // Loop over the rows of the matrix A. 
  for (int row=0; row < N; row++) {
    // Loop over every column of B. 
    for (int colB=0; colB < N; colB++) {
      // Reset result value. 
      r = 0;
      // Loop over all columns of A and multiply with rows of B.  
      for (int colA=0; colA < N; colA++) {
        r = r + A[colA + row*N] * B[colA*N + colB];
      }
      // Assign computed value to the result matrix C. 
      C[row*N + colB] = r;
    }
  }

  // // Print matrix.
  // printf("Matrix normal multiply\n"); 
  // for (int row = 0; row < N; row++) {
  //   for (int col = 0; col < N; col++) {
  //     printf("%f     ", C[row*N + col]);
  //   }
  //   printf("\n");
  // }
  // printf("\n");

}


void AB_block_row( std::vector<double> A, std::vector<double> B, size_t blockSize ){

  size_t N = sqrt(A.size());
  std::vector<double> C(N*N);

  // TODO: Question 2c: Block matrix-matrix multiplication - B in row major

  // Store results of partial block multiplications. 
  double r = 0;

  // Moving block over rows. 
  for (int rowBlock = 0; rowBlock < N; rowBlock += blockSize) {
    // Moving block over columns. 
    for (int colBlock = 0; colBlock < N; colBlock += blockSize) {
      // Moving inside of the block. 
      // Loop through rows of block.  
      for (int row = 0; row < N; row++) {
        // Loop through cols in a row of the block. 
        for (int col = rowBlock; col < rowBlock + blockSize; col++) {
          // Matrix block multiplication. 
          r = 0;
          for (int colA = colBlock; colA < colBlock + blockSize; colA++) {
            r = r + A[colA + row*N] * B[colA*N + col];
          }
          // Assign new value to the result matrix C. 
          C[row*N + col] += r;
        }
      }

      // // Print matrix.
      // printf("Matrix block multiply steps:\n"); 
      // for (int row = 0; row < N; row++) {
      //   for (int col = 0; col < N; col++) {
      //     printf("%f     ", C[row*N + col]);
      //   }
      //   printf("\n");
      // }
      // printf("----------------------\n");

    }
  }

  // // Print matrix.
  // printf("Matrix block multiply\n"); 
  // for (int row = 0; row < N; row++) {
  //   for (int col = 0; col < N; col++) {
  //     printf("%f     ", C[row*N + col]);
  //   }
  //   printf("\n");
  // }
  // printf("\n");

}


void AB_block_col( std::vector<double> A, std::vector<double> B, size_t blockSize ){
  
  size_t N = sqrt(A.size());
  std::vector<double> C(N*N);

  // TODO: Question 2c: Block matrix-matrix multiplication - B in column major

  // Store results of partial block multiplications. 
  double r = 0;

  // Moving block over rows. 
  for (int rowBlock = 0; rowBlock < N; rowBlock += blockSize) {
    // Moving block over columns. 
    for (int colBlock = 0; colBlock < N; colBlock += blockSize) {
      // Moving inside of the block. 
      // Loop through rows of block.  
      for (int row = 0; row < N; row++) {
        // Loop through cols in a row of the block. 
        for (int col = rowBlock; col < rowBlock + blockSize; col++) {
          // Matrix block multiplication. 
          r = 0;
          for (int colA = colBlock; colA < colBlock + blockSize; colA++) {
            r = r + A[colA + row*N] * B[colA + col*N];
          }
          // Assign new value to the result matrix C. 
          C[row*N + col] += r;
        }
      }

      // // Print matrix.
      // printf("Matrix block multiply steps:\n"); 
      // for (int row = 0; row < N; row++) {
      //   for (int col = 0; col < N; col++) {
      //     printf("%f     ", C[row*N + col]);
      //   }
      //   printf("\n");
      // }
      // printf("----------------------\n");

    }
  }
}


double benchmark_AB( std::vector<double> A, std::vector<double> B, size_t mode, size_t blockSize, size_t Ns ){

  size_t N = sqrt(A.size());
  double times = 0;

  // TODO: Check that the matrix size is divided by the blockSize when mode==2 or 3
  if( (mode==2 or mode==3) &&  N%blockSize!=0 ){
    printf("Error: the size of the matrix (%zu) should be divided by the blockSize variable (%zu).\n",N,blockSize);
    exit(1);
  }

  for( size_t i=0; i<Ns; i++){
    auto t1 = std::chrono::system_clock::now();
    
    // TODO: Question 2c: Call the function to be benchmarked
    
    if( mode==1 ){
      AB(A, B);
    }
    else if( mode==2 ){
      AB_block_row(A, B, blockSize);      
    }
    else if( mode==3 ){
      AB_block_col(A, B, blockSize);
    }
    auto t2 = std::chrono::system_clock::now();
    times += std::chrono::duration<double>(t2-t1).count();
  }
  printf("Done in total %9.4fs  --  average %9.4fs\n", times, times/Ns);

  return times/Ns;

}


int main( )
{

  std::vector<int> matrixSize{ 256, 512}; //, 1024, 2048  };
  size_t M = matrixSize.size();

  std::vector<size_t> blockSize{ 2, 4, 8, 16, 32, 64, 128 };
  size_t Bs = blockSize.size();

  size_t Ns = 2;   // ORIGINAL VALUE IS 5; CHANGED FOR TESTING

  std::vector<double> times1(M);
  std::vector< std::vector<double>> times2(Bs, std::vector<double>(M) );
  std::vector< std::vector<double>> times3(Bs, std::vector<double>(M) );

  std::vector<double> A, B, C;


  for( size_t m=0; m<M; m++){

    printf("Working with matrices of size %d\n",matrixSize[m]);
    printf("---------------------------------------------\n");

    size_t N = matrixSize[m];

    // TODO: Question 2c: Initialize matrices
    //       store A and B as row major and C as column major

    // Initialize the matrices A and B: M_i,j = i + j
    // Store matrix in row major order. 
    for(int row=0; row < N; row++) {
      for (int col=0; col < N; col++) {
        A.push_back(row+col);
        B.push_back(2*row+col);
      }
    }

    // Initialize the matrix C. 
    // Store matrix in column major order.
    // C is B transposed.  
    C = transpose(B);

    printf("----------MMMMMMMMMM------------\n");
    printf("%f\n", C[0]);

    // printf("Size of A = %d\n", A.size());
    // printf("Size of A = %f\n", sqrt(A.size()));
    // printf("Size of B = %f\n", sqrt(B.size()));
    // printf("Size of C = %f\n", sqrt(C.size())); = %f\n", sqrt(A.size()));
    // printf("Size of B = %f\n", sqrt(B.size()));
    // printf("Size of C = %f\n", sqrt(C.size()));


    printf("Start C=A*B (non optimized).\n");
    times1[m] = benchmark_AB( A, B, 1, 0, Ns );

    printf("---------------------------------------------\n");

    for( size_t b=0; b<Bs; b++){
      printf("Start C=A*B (optimized, row major, block size=%zu).\n", blockSize[b] );
      times2[b][m] = benchmark_AB( A, B, 2, blockSize[b], Ns );
    }

    printf("---------------------------------------------\n");

    for( size_t b=0; b<Bs; b++){
      printf("Start C=A*B (optimized, column major, block size=%zu).\n", blockSize[b] );
      times3[b][m] = benchmark_AB( A, C, 3, blockSize[b], Ns );
    }

    printf("==================================================\n");

    // Clear vectors for the next iteration. 
    A.clear();
    B.clear();
    C.clear();
  }



  // Write resuts to file. 

  // FILE *fp=nullptr;
  // fp = fopen("matrix_matrix_times.txt","w");
  // // write header to the file
  // std::string header = " N   unopt ";
  // for(size_t b=0; b<Bs; b++)
  //   header = header + "  br_" + std::to_string(blockSize[b]);
  // for(size_t b=0; b<Bs; b++)
  //   header = header + "  bc" + std::to_string(blockSize[b]);
  // header = header + "\n";
  // fprintf(fp,"%s",header.c_str());

  // for(size_t m=0; m<M; m++){
  //   fprintf(fp,"%d %lf",matrixSize[m],times1[m]);
  //   for(size_t b=0; b<Bs; b++)
  //     fprintf(fp," %lf ",times2[b][m]);
  //   for(size_t b=0; b<Bs; b++)
  //     fprintf(fp," %lf ",times3[b][m]);
  //   fprintf(fp,"\n");
  // }
  // fclose(fp);


  return 0;
}
