// Copyright 2020 ETH Zurich. All Rights Reserved.

#include <iostream>
#include <iomanip>
#include <random>
#include <cmath>
#include <mpi.h>

inline long exact(const long N){
    // TODO b): Implement the analytical solution.
    // Gauss sum formula. 
    long r = (N*(N-1))/2;
    return r;
}

void reduce_mpi(const int rank, long& sum){
    // TODO e): Perform the reduction using blocking collectives.
    MPI_Reduce(rank ? &sum : MPI_IN_PLACE, &sum, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
}

// PRE: size is a power of 2 for simplicity
void reduce_manual(int rank, int size, long& sum){
    // TODO f): Implement a tree based reduction using blocking point-to-point communication.
    long recvSum = 0;
    for (int i = size; i > 1; i/=2) {
        // Ranks above I have already been computed. 
        if (rank < i) {
            // Divide ranks in lower and upper half. 
            if (rank >= i/2) {
                // Upper half ranks sends their sums to the lower half ranks
                MPI_Ssend(&sum, 1, MPI_LONG, rank-i/2, 420, MPI_COMM_WORLD); 
            }
            else {
                // Lower half ranks receive the sums of the upper half ranks. 
                MPI_Recv(&recvSum, 1, MPI_LONG, rank+i/2, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                sum += recvSum;
                // Set temporary sum to zero again, results stored must not add up.
                recvSum = 0;
            }
        }
    }
}


int main(int argc, char** argv){
    const long N = 1000000;
    
    // Initialize MPI
    int rank, size;
    // TODO c): Initialize MPI and obtain the rank and the number of processes (size)
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size); 
   
    // -------------------------
    // Perform the local sum:
    // -------------------------
    long sum = 0;
     
    // Determine work load per rank
    long N_per_rank = N / size;
    
    // TODO d): Determine the range of the subsum that should be calculated by this rank.
    long N_start = N_per_rank*rank;
    long N_end = N_start + N_per_rank;
    
    // N_start + (N_start+1) + ... + (N_start+N_per_rank-1)
    for(long i = N_start; i < N_end; ++i){
        sum += i;
    }
    
    // -------------------------
    // Reduction
    // -------------------------
    //reduce_mpi(rank, sum);
    reduce_manual(rank, size, sum);
    
    // -------------------------
    // Print the result
    // -------------------------
    if(rank == 0){
        std::cout << std::left << std::setw(25) << "Final result (exact): " << exact(N) << std::endl;
        std::cout << std::left << std::setw(25) << "Final result (MPI): " << sum << std::endl;
    }
    // Finalize MPI
    // TODO c): Finalize MPI
    MPI_Finalize();
    
    return 0;
}
