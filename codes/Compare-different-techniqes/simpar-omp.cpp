#include <iostream>
#include <stdlib.h>
#include <omp.h>
#include "funcDef.h"
#include "simpar-omp.h"
#include "simpar-omp-reduction.h"
#include "simpar-omp-atomic.h"

using namespace std;

int main(int argc, char **argv)
{
    long seed;
    unsigned int ncside, ntstep;
    size_t n_part;

    // Function argument assignments
    seed = atol(argv[1]);    //  seed for the random number generator
    ncside = atol(argv[2]);  //  size of the grid (number of cells on the side)
    n_part = atoll(argv[3]); //  number of particles
    ntstep = atol(argv[4]);  //  number of time steps

    if (n_part / (ncside * ncside) > 5)
    {
        if (ncside < 150)
            reduction(seed, ncside, n_part, ntstep);

        else
            manual(seed, ncside, n_part, ntstep);
    }
    else
    {
        atomic(seed, ncside, n_part, ntstep);
    }
        return 0;
}
