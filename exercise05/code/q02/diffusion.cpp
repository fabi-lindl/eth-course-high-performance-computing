#include <fstream>
#include <iostream>
#include <mpi.h>
#include <stdio.h>
#include <vector>
#include <math.h>

struct Diagnostics {
    double time;
    double concentration;
    Diagnostics(double time, double concentration) :
        time(time), concentration(concentration) {}
};

struct Diffusion
{
    double D, L;  //diffusion constant and domain length
    int N;        //grid points per direction (whole grid is NxN)
    int local_N;  //number of rows of this process

    double h,dt;  //grid spacing and timestep
    double aux;   //auxiliary variable

    std::vector<double> c;    //solution vector
    std::vector<double> c_tmp;

    int rank,size; //MPI rank and total number of ranks

    std::vector<Diagnostics> diag;

    Diffusion(double D, double L, int N, int rank, int size) :
        D(D), L(L), N(N), rank(rank), size(size)
    {
        h = L / (N - 1);
        dt = h*h/(4.0*D); // Largest possible timestep (larger values lead to instabilities).

        local_N = N / size;
        if (rank == size - 1) local_N += N % size; // Correction for the last process.

        c.resize((local_N+2)*(N+2), 0.0); // +2 for the ghost cells
        c_tmp.resize((local_N+2)*(N+2), 0.0);

        aux = dt * D / (h*h);
        initialize_density();
    }

    void advance()
    {
        // TODO: Implement Blocking MPI communication to exchange the ghost
        // cells on a periodic domain required to compute the central finite
        // differences below.

        // In case only one rank is used, the send receive must not be 
        // executed, it is going to cause an error. 
        if (size > 1) {
            
            // Define number of request according to the rank number. 
            int nr;
            if (rank > 0 && rank < (size-1))
                nr = 4;
            else
                // Rank 0 and rank (size-1) are treated differently.
                // These two ranks do only send one message. 
                // Rank 0 sends it to rank 1 and rank (size-1) sends it to rank (size-2). 
                nr = 2;
            
                // Define number of send and receive requests. 
            std::vector<MPI_Request> requests(nr, MPI_REQUEST_NULL);	
            // Define the statuses. 
            std::vector<MPI_Status> stats(nr);
            // Define a tag for this ranke. 
            int tag = rank;

            // Send messages.
            int suIdx = (N+2) * local_N; // Send up index. 
            int sdIdx = N+2;             // Send down index. 
            
            // Send down. 
            if (rank > 0) 
                MPI_Isend(&c[sdIdx], N+2, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD, &requests[0]);
            // Send up. 
            if (rank == 0)
                MPI_Isend(&c[suIdx], N+2, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD, &requests[0]);
            else if (rank < (size-1))
                MPI_Isend(&c[suIdx], N+2, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD, &requests[1]);
            // Receive messages. 
            int raIdx = (N+2) * (local_N + 1); // Receive from above index. 
            int rbIdx = 0;                     // Receive from below index. 	

            // Receive from below. 
            if (rank == (size-1))
                MPI_Irecv(&c[rbIdx], N+2, MPI_DOUBLE, rank-1, tag-1, MPI_COMM_WORLD, &requests[1]);
            else if (rank > 0)
                MPI_Irecv(&c[rbIdx], N+2, MPI_DOUBLE, rank-1, tag-1, MPI_COMM_WORLD, &requests[2]);
            // Receive from above
            if (rank == 0)
                MPI_Irecv(&c[raIdx], N+2, MPI_DOUBLE, rank+1, tag+1, MPI_COMM_WORLD, &requests[1]); 
            else if (rank < (size-1))
                MPI_Irecv(&c[raIdx], N+2, MPI_DOUBLE, rank+1, tag+1, MPI_COMM_WORLD, &requests[3]);
            
            // Wait until all requests have completed. 
            MPI_Waitall(nr, requests.data(), stats.data());	
	    }

	    /* Central differences in space, forward Euler in time, Dirichlet BCs */
        for (int i = 1; i <= local_N; ++i)
            for (int j = 1; j <= N; ++j)
                c_tmp[i * (N + 2) + j] =
                    c[i * (N + 2) + j] +
                    aux * (c[i * (N + 2) + (j + 1)] + c[i * (N + 2) + (j - 1)] +
                           c[(i + 1) * (N + 2) + j] + c[(i - 1) * (N + 2) + j] -
                           4 * c[i * (N + 2) + j]);

        // Use swap instead of rho_ = rho_tmp__. This is much more efficient,
        // because it does not copy element by element, just replaces storage
        // pointers. 
        using std::swap;
        swap(c_tmp, c);
    }

    void compute_diagnostics(const double t)
    {
        /*
         * Sum the concentration of all local grids of the different 
         * ranks and sum them up. 
         * Write the results to the diagnostics data structure. 
         */
        // Storage variable for the concentration amount of this iteration. 
        double ammount = 0.0;

        /* Integration to compute total concentration */
        for (int i = 1; i <= local_N; ++i)
            for (int j = 1; j <= N; ++j)
                ammount += c[i * (N + 2) + j] * h * h;

        // TODO: Sum total ammount from all ranks
	    MPI_Reduce(rank ? &ammount : MPI_IN_PLACE, &ammount, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        if (rank == 0) {
            std::cout << "t = " << t << " ammount = " << ammount << '\n';
            diag.push_back(Diagnostics(t, ammount));
        }
    }


    void write_diagnostics(const std::string& filename) const
    {
        /*
        * Write concentration at every grid point at a certain
        * point in time to a CSV file.
        */
        std::ofstream out_file(filename, std::ios::out);
        for (const Diagnostics& d : diag)
            out_file << d.time << ' ' << d.concentration << '\n';
        out_file.close();
    }


    void compute_histogram()
    {
        /*
        * Distributes the concentration of all the small grid
        * squares into a histogram of M bins.
        */
        /* Number of histogram bins */
        const int M = 10;
        std::vector<int> hist(M, 0);

        /* Find the local max and min density values */
        double max_c, min_c, c0;
        max_c = c[1 * (N + 2) + 1];
        min_c = c[1 * (N + 2) + 1];

        for (int i = 1; i <= local_N; ++i)
            for (int j = 1; j <= N; ++j) {
	        c0 = c[i * (N + 2) + j];
                if (c0 > max_c)
                    max_c = c0;
                if (c0 < min_c)
                    min_c = c0;
            }

        // TODO: Compute the global min and max concentration values on this myRank and
        // store the result in min_c and max_c, respectively.
        MPI_Allreduce(rank ? &min_c : MPI_IN_PLACE, &min_c, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(rank ? &max_c : MPI_IN_PLACE, &max_c, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        double epsilon = 1e-8;
        double dh = (max_c - min_c + epsilon) / M;

        for (int i = 1; i <= local_N; ++i)
            for (int j = 1; j <= N; ++j) {
	            // Compute histrogram bin index. 
                int bin = (c[i * (N + 2) + j] - min_c) / dh;
          	    // Update histogram at a particular bin. 
	            hist[bin]++;
            }

        // TODO: Compute the sum of the histogram bins over all ranks and store
        // the result in the array g_hist.  Only myRank 0 must print the result.
        std::vector<int> g_hist(M, 0);
    	MPI_Reduce(&hist[0], &g_hist[0], M, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        for (int i = 0; i < size; i++) {
            if (rank == i) {
            printf("\nRank %d\n", rank);
                for (int l = 0; l < M; l++) 
                    printf("bin[%d] = %d\n", l, hist[l]);
            }
        }
	
	    // Print results (histogram of all the ranks combined). 
        if (rank == 0) {
            printf("=====================================\n");
            printf("Output of compute_histogram():\n");
            int gl = 0;
	        // Print results to console. 
            for (int i = 0; i < M; i++) {
                printf("g_bin[%d] = %d\n", i, g_hist[i]);
                gl += g_hist[i];
            }
	        // Print the total number of elements as a check. 
            printf("Total elements = %d\n", gl);
        }

    }

    void initialize_density()
    {
        /*
        * Density is 1 in the inner square of length 0.25*L.
        * Outside ofthis boundary the density is zero. 
        */
	    // Global index for the local ranks. 
        int gi;
	    // Boundary of the inner square. 
        double bound = 0.25 * L;

        for (int i = 0; i < local_N; ++i) {
	        // Convert local index to global index of the grid. 
            gi = rank * (N / size) + i;
            for (int j = 0; j < N; ++j) {
                if (fabs(gi * h - 0.5*L) < bound && fabs(j * h - 0.5*L) < bound)
                    c[(i+1) * (N+2) + (j+1)] = 1; // Set density to one at this position. 
                else
                    c[(i+1) * (N+2) + (j+1)] = 0; // Set density to zero at this position. 
            }
        }
    }
};

int main(int argc, char* argv[])
{
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " D L N \n";
        return 1;
    }

    // TODO: Start-up the MPI environment and determine this process' myRank ID
    // as well as the total number of processes (=ranks) involved in the
    // communicator
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Variable to check timing t = 0.5 s. 
    int check = 1;

    // Diffusion constant.
    const double D = std::stod(argv[1]);
    // Length of the square in which the substance diffuses.
    const double L = std::stod(argv[2]);
    // Number of grid points (square grid). 
    const int N = std::stoul(argv[3]);

    if (rank == 0)
        printf("Running Diffusion 2D on a %d x %d grid with %d ranks.\n",N,N,size);

    Diffusion system(D, L, N, rank, size);
    system.compute_diagnostics(0);
    for (int step = 0; step < 10000; ++step) {
	    system.advance();
        system.compute_diagnostics(system.dt * step);
        // Print results of time t = 0.5 s to the console. 
        if (check == 1 && (system.dt * step) > 0.5) {
            check = 0;
            system.compute_histogram();
        }
    }
    system.compute_histogram();
    if (rank == 0)
        system.write_diagnostics("diagnostics.dat");
    
    // TODO: Shutdown the MPI environment
    MPI_Finalize();
    
    return 0;
}
