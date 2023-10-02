#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <chrono>
#include <math.h>
#include <string>


void transpose( std::vector<double> A ){
  
  size_t N = sqrt(A.size());
  std::vector<double> AT(N*N);
  
  // TODO: Question 2b: Straightforward matrix transposition

  // // Slow.
  // // Matrix stored in row major order, loops over the columns. 
  // for(int row = 0; row < N; row++) {
  //   for(int col = 0; col < N; col++) {
  //     AT.at(col+row*N) = A.at(row+col*N);
  //   }
  // }

  // Fast
  // Matrix stored in row major order, loops over the rows. 
  for(int col = 0; col < N; col++) {
    for(int row = 0; row < N; row++) {
      AT.at(col+row*N) = A.at(row+col*N);
    }
  }

}


void transpose_block( std::vector<double> A, size_t blockSize ){

  size_t N = sqrt(A.size());
  std::vector<double> AT(N*N);

  // Print the matrix. 
  // printf("\nMatrix\n");
  // for (int i=0; i < N*N; i++) {
  //   printf("%.1f    ", A.at(i));
  // }


  // TODO: Question 2b: Block matrix transposition

  // The matrix is divided into square block of blockSize. 
  // The first block that is computed is in the upper left corner.
  // Then the block is shifted rightwards until the block reaches the
  // right edge of the matrix. 
  // After that the block is moved down by one block size and to the 
  // left edge of the matrix again. 
  // This procedure is used over and over again until the whole matrix
  // is computed. 
  // Moving block over rows. 
  for (int rowBlock = 0; rowBlock < N; rowBlock += blockSize) {
    // Moving block over columns. 
    for (int colBlock = 0; colBlock < N; colBlock += blockSize) {
      // Moving inside the block. 
      // Rows. 
      for (int row = rowBlock; row < rowBlock + blockSize; row++) {
          // Cols. 
          for (int col = colBlock; col < colBlock + blockSize; col++) {
              AT[row + col*N] = A[col + row*N];
          }
      }
    }
  }

  // Print the transposed matrix. 
  // printf("\nTransposed matrix\n");
  // for (int i=0; i < N*N; i++) {
  //   printf("%.1f    ", AT.at(i));
  // }
  
}


double benchmark_transpose( std::vector<double> A, size_t mode, size_t blockSize, size_t Ns ){

  size_t N = sqrt(A.size());
  double times = 0;

  // TODO: Check that the matrix size is divided by the blockSize when mode==2
  // I guess this is not a todo, but rather just a flag to avoid any matrices of sizes of 
  // non multiples of the blockSize. 
  if( mode==2 &&  N%blockSize!=0 ){
    printf("Error: the size of the matrix (%zu) should be divided by the blockSize variable (%zu).\n",N,blockSize);
    exit(1);
  }

  for( size_t i=0; i<Ns; i++){
    auto t1 = std::chrono::system_clock::now();
    // TODO: Call the function to be benchmarked
    if( mode==1 ){
      transpose(A);
    }
    else if( mode==2 ){
      transpose_block(A, blockSize);
    }
    auto t2 = std::chrono::system_clock::now();
    times += std::chrono::duration<double>(t2-t1).count();
  }
  printf("Done in total %9.4fs  --  average %9.4fs\n", times, times/Ns);

  return times/Ns;

}


int main( )
{
  std::vector<int> matrixSize{ 1024, 2048, 4096 };
  size_t M = matrixSize.size();

  std::vector<size_t> blockSize{ 2, 4, 8, 16, 32, 64, 128 };
  size_t B = blockSize.size();

  size_t Ns = 2;

  std::vector<double> times1(M);
  std::vector< std::vector<double>> times2(B, std::vector<double>(M) );

  // Matrix to be transposed. 
  std::vector<double> A;

  // loop over matrix sizes
  for( size_t m=0; m<M; m++){

    printf("Working with a matrix of size %d\n",matrixSize[m]);

    size_t N = matrixSize[m];
    
    // TODO:
    
    // Initialize the matrix: A_i,j = i + j; x_i = i
    // Store matrix in row major order. 
    for(int row=0; row < N; row++) {
      for (int col=0; col < N; col++) {
        A.push_back(row+col);
      }
    }

    printf("Start transposing (non optimized).\n");
    times1[m] = benchmark_transpose( A, 1, 0, Ns );

    // loop over block sizes
    for( size_t b=0; b<B; b++){
      printf("Start transposing (optimized, block size=%zu).\n", blockSize[b] );
      times2[b][m] = benchmark_transpose( A, 2, blockSize[b], Ns );
    }

    // Clear vector for the next iteration. 
    A.clear();

    printf("==================================================\n");
  }


  // // write results to a file
  // FILE *fp=nullptr;
  // fp = fopen("transpose_times.txt","w");
  // // write header to the file
  // std::string header = "# N   time_unoptimized ";
  // for(size_t b=0; b<B; b++)
  //   header = header + "  block_" + std::to_string(blockSize[b]);
  // header = header + "\n";
  // fprintf(fp,"%s",header.c_str());
  // for(size_t m=0; m<M; m++){
  //   fprintf(fp,"%d %lf",matrixSize[m],times1[m]);
  //   for(size_t b=0; b<B; b++)
  //     fprintf(fp," %lf ",times2[b][m]);
  //   fprintf(fp,"\n");
  // }
  // fclose(fp);

  return 0;
}
