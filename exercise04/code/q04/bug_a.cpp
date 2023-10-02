#include <fstream>
#include <mpi.h>

int main(int argc, char *argv[]) 
{
    // Initialize MPI ranks.
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int tag = 420;
    const int N = 20;
    double *result = new double[N];
    double *printBuf; // Buffer for writing data to file.

    // Write data to array. 
    for (int i = 0; i < N; i++)
        result[i] = (double)i*rank;

    // Write results to file. 
    if (rank == 0)
	    printBuf = new double[size*N];
    MPI_Gather(result, N, MPI_DOUBLE, printBuf, N, MPI_DOUBLE, 0, MPI_COMM_WORLD); 

    if (rank == 0) {
        std::ofstream file("bug_a_results.txt");
        int lim = size*N;
        for (int i = 0; i < lim; i++)
            file << printBuf[i] << std::endl;
        file.close();
        delete[] printBuf;
    }

    delete[] result;
    MPI_Finalize();

    return 0;
}
