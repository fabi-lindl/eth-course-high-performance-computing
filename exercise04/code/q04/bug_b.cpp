#include <mpi.h>

int main(int argc, char *argv[])
{
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double important_value;

    if (rank == 0)
	    important_value = 0;
    else
	    important_value = 1.1;

    if (rank == 1)
        MPI_Send(&important_value, 1, MPI_DOUBLE, 0, 123, MPI_COMM_WORLD);
    else
	    MPI_Send(&important_value, 1, MPI_DOUBLE, 1, 123, MPI_COMM_WORLD);

    MPI_Recv(&important_value, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, 
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    printf("Rank %d received the double: %f\n", rank, important_value);

    MPI_Finalize();

    return 0;
}
