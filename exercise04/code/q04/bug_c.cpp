#include <mpi.h>

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int bval;
    
    if (rank == 0)
	    val = 10;

    MPI_Bcast(&bval, 1, MPI_INT, 0, MPI_COMM_WORLD);

    std::cout << "[" << rank << "] " << bval << std::endl;
    
    MPI_Finalize();

    return 0;
}
