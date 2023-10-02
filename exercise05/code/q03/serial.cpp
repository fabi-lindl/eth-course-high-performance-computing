#define _USE_MATH_DEFINES

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

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

static void computeVelocities(double epsSq, const std::vector<double>& x,
                              const std::vector<double>& y,
                              const std::vector<double>& gamma,
                              std::vector<double>& u, std::vector<double>& v)
{
    // TODO compute the interactions between the particles and write the result
    // in the velocity vectors u and v. 

    // Every particle must interact with all other particles. 
    
    // Declar and define variables. 
    double xsum, ysum, xi, yi, r;
    int lim = y.size();

    for (int i = 0; i < lim; i++) {
        // Storage variables for the updated coordinates of every particle interaction. 
        xsum = 0.0;
        ysum = 0.0;
        // Set variables that interact with all other particles. 
        xi = x[i];
        yi = y[i];
        // Loop over all particles for the interaction. 
        for (int k = 0; k < lim; k++) {
            // Descritized form of the PDE for diffusion.
            r = gamma[k] / (2 * M_PI * (epsSq + (xi-x[k])*(xi-x[k]) + (yi-y[k])*(yi-y[k])));
            xsum += - (yi-y[k]) * r;
    	    ysum += (xi-x[k]) * r;
	    }
        // Update the velocity values. 
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

int main(int argc, char** argv)
{
    // Total number of particles must have a modulo 2 of zero. 
    if (argc != 2) {
        fprintf(stderr, "usage: %s <total number of particles>\n", argv[0]);
        exit(1);
    }

    const int n = std::atoi(argv[1]);

    // Start and end values of the 1D initialized particles. 
    const double startX = -0.5;
    const double endX = 0.5;

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

    // Initialize x, y, and gamma. 
    initialConditions(startX, endX, x, y, gamma);

    for (int step = 0; step < numSteps; ++step) {
        if (step % dumpEvery == 0) {
            // Write results to a csv file (24 files are created in total). 
            const int id = step / dumpEvery;
            dumpToCsv(id, x, y, gamma);
        }
        // Compute velocities u, v. 
        computeVelocities(epsSq, x, y, gamma, u, v);
        // Go ahead one step in time (dt). 
        forwardEuler(dt, u, v, x, y);
    }

    return 0;
}
