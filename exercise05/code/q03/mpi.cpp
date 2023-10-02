#define _USE_MATH_DEFINES

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
#include <mpi.h>

static double gammaInit(double x)
{
    /*
     * Initialization of the circulation variable gamm. 
     * Init. according to the formula on the assignment sheet. 
    */
    return 4 * x / std::sqrt(1 - 4 * x * x);
}

static void initialConditions(double start, double end, std::vector<double>& x,
                              std::vector<double>& y,
                              std::vector<double>& gamma)
{
    // Number of vector elements, i.e. number of particles.
    const int n = (int)gamma.size();
    // Distance between particles. 
    const double dx = (end - start) / n;

    // Loop over all values of x, y, and gamma and update their values. 
    // y is zero at the start for all particles, this value changes with
    // every step due to the circulation gamma. 
    for (int i = 0; i < n; ++i) {
        // x is initialized to a formula on the assignment sheet. 
        x[i] = start + dx * (i + 0.5);
        y[i] = 0;
        gamma[i] = dx * gammaInit(x[i]);
    }
}

static void computeVelocities(MPI_Comm comm, double epsSq,
                              std::vector<double>& x,
                              std::vector<double>& y,
                              std::vector<double>& gamma,
                              std::vector<double>& u, std::vector<double>& v)
{
    // TODO: perform multi pass to compute interactions and update the local
    // velocities.
    // Test if MPI vars are available. 
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    
    // Every rank needs all the data of all other ranks. 
    // Use all gather so that the data of every rank is available at every rank. 
    // Number of particles worked with. 
    int numParticles = x.size();
    int numgParticles = numParticles*size;
    // Data receive storage for gathered values. 
    std::vector<double> xg(numgParticles);
    std::vector<double> yg(numgParticles);
    std::vector<double> gg(numgParticles);
    
    // MPI_Allgather to pass data to every rank. 
    MPI_Allgather(&x[0], numParticles, MPI_DOUBLE, &xg[0], numParticles, MPI_DOUBLE, comm);
    MPI_Allgather(&y[0], numParticles, MPI_DOUBLE, &yg[0], numParticles, MPI_DOUBLE, comm);
    MPI_Allgather(&gamma[0], numParticles, MPI_DOUBLE, &gg[0], numParticles, MPI_DOUBLE, comm);
    
    // Compute velocities of this rank. 
    // Declare variables. 
    double xsum, ysum, xi, yi, r;
    // Loop over particle values of this rank. 
    for (int i = 0; i < numParticles; i++) {
        // Let the particles of this rank interact with all other particles. 
        xsum = 0.0;
        ysum = 0.0;
        // Define particles that interact with all other particles. 
        xi = x[i];
        yi = y[i];

        for (int k = 0; k < numgParticles; k++) {
            r = gg[k] / (2 * M_PI * (epsSq + (xi-xg[k])*(xi-xg[k]) + (yi-yg[k])*(yi-yg[k])));
            xsum += - (yi-yg[k]) * r;
            ysum += (xi-xg[k]) * r;
        }
        u[i] = xsum;
        v[i] = ysum;
    }
}

static void forwardEuler(double dt, const std::vector<double>& u,
                         const std::vector<double>& v, std::vector<double>& x,
                         std::vector<double>& y)
{
    /*
     * Make one step in time with the forward Euler method. 
    */ 
    for (int i = 0; i < (int)x.size(); ++i) {
        x[i] += dt * u[i];
        y[i] += dt * v[i];
    }
}

static void dumpToCsv(int step, const std::vector<double>& x,
                      const std::vector<double>& y,
                      const std::vector<double>& gamma)
{
    /*
     * Write x, y, and gamma results to a csv file. 
    */
    char fname[128];
    sprintf(fname, "config_%05d.csv", step);

    FILE* f = fopen(fname, "wb");
    fprintf(f, "x,y,gamma\n");
    for (int i = 0; i < (int)x.size(); ++i) {
        fprintf(f, "%g,%g,%g\n", x[i], y[i], gamma[i]);
    }
    fclose(f);
}

static void dumpToCsv(MPI_Comm comm, int step, std::vector<double>& x,
                      std::vector<double>& y,
                      std::vector<double>& gamma)
{
    int rank, nranks;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nranks);

    int numParticles = x.size();                // Particles of one rank. 
    int numgParticles = numParticles * nranks;  // Gathered particles. 

    std::vector<double> xAll(numgParticles), yAll(numgParticles), gammaAll(numgParticles);

    // TODO Gather the data on rank zero before dumping to the csv files.
    MPI_Gather(&x[0], numParticles, MPI_DOUBLE, &xAll[0], numParticles, MPI_DOUBLE, 0, comm);
    MPI_Gather(&y[0], numParticles, MPI_DOUBLE, &yAll[0], numParticles, MPI_DOUBLE, 0, comm);
    MPI_Gather(&gamma[0], numParticles, MPI_DOUBLE, &gammaAll[0], numParticles, MPI_DOUBLE, 0, comm);
    
    if (rank == 0)
        dumpToCsv(step, xAll, yAll, gammaAll);
}

int main(int argc, char** argv)
{
    // Init. MPI.
    MPI_Init(&argc, &argv);

    // Total number of particles must have a modulo 2 of zero.
    if (argc != 2) {
        fprintf(stderr, "usage: %s <total number of particles>\n", argv[0]);
        exit(1);
    }
    
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, nranks;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nranks);

    const int nglobal = std::atoi(argv[1]);

    if (nglobal % nranks != 0) {
        fprintf(stderr,
                "expected n to be a multiple of the number of ranks.\n");
        exit(1);
    }

    // TODO initialize the data for each rank.
    // Each rank needs to work on the same number of particles. 
    const int n = nglobal / nranks;              // TODO
    const double extents = 1.0 / nranks;         // TODO
    const double startX = -0.5 + extents * rank; // TODO
    const double endX = startX + extents;

    // Size of a timestep. 
    const double dt = 1e-4;
    // Small value to prevent a denominator of zero. 
    const double epsSq = 1e-3;
    const double endTime = 2.5;
    // Frequency of how often the results are written to a csv file.
    const double dumpFreq = 0.1;

    // Write results to csv file every 1000 steps. 
    const int dumpEvery = dumpFreq / dt;
    const int numSteps = endTime / dt;

    // state of the simulation
    std::vector<double> x(n);
    std::vector<double> y(n);
    std::vector<double> gamma(n);

    // workspace
    std::vector<double> u(n);
    std::vector<double> v(n);    

    // Init. x, y, and gamma. 
    initialConditions(startX, endX, x, y, gamma);

    for (int step = 0; step < numSteps; ++step) {
        if (step % dumpEvery == 0) {
            // Write results to a csv file (24 files are created in total).
            const int id = step / dumpEvery;
            dumpToCsv(comm, id, x, y, gamma);
        }
        // Compute updated velocites u, v.
        computeVelocities(comm, epsSq, x, y, gamma, u, v);
        // Go ahead one step in time (dt). 
        forwardEuler(dt, u, v, x, y);
    }

    MPI_Finalize();

    return 0;
}
