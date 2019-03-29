#include <iostream>
#include <stdlib.h>
#include "funcDef.h"
#include "simpar.h"
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
    serial(seed, ncside, n_part, ntstep);   
    return 0;
}
