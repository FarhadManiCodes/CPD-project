#include <iostream>
#include <stdlib.h>
//#include <vector>
#include <ctime>
#include "funcDef.h"


using namespace std;

int main(int argc, char **argv)
{
    long seed; 
    unsigned int ncside, ntstep;
    unsigned long n_part;

    // Function argument assignments
    seed = atol(argv[1]);   //  seed for the random number generator
    ncside = atol(argv[2]); //  size of the grid (number of cells on the side)
    n_part = atoll(argv[3]); //  number of particles
    ntstep = atol(argv[4]); //  number of time steps

    clock_t time_req;
    time_req = clock();

    // Declarations
    cell_t *cell = (cell_t *)calloc(ncside*ncside,sizeof(cell_t)); // Matrix containing cells of the problem
    particle_t *par =  (particle_t *) malloc(n_part * sizeof(particle_t)); // vector containing all particles of the problem

    // Initialize particles and cells
    init_particles(seed, ncside, n_part, par, cell);    
       
    // Loop over time
    for (unsigned int t_step = 0; t_step < ntstep; t_step++)
    {

        cell_t *cell_aux = (cell_t *)calloc(ncside*ncside,sizeof(cell_t)); // Auxilary matrix containing cells of the problem for the next time step 

        // Loop over particles
        for (unsigned long i = 0; i < n_part; i++)
        {
            unsigned int ci = par[i].c_i; //cell index i
            unsigned int cj = par[i].c_j; //cell index j
            double xp = par[i].x;          // x of the particle
            double yp = par[i].y;          // y of the particel
            double m = par[i].m;           // mass of the particle
            double Fx = 0.0, Fy = 0.0;     // Fx,Fy force in (x,y) direction

            // Calculate force components
            calculate_forces(ci, cj, ncside, xp, yp, m, Fx, Fy, cell);

            // Update particle positions
            update_velocities_and_positions(i, Fx, Fy, par);

            // Update particle's cell info
            locate_and_update_cell_info(i, ncside, par, cell_aux);
        }// end of loop over particles
        // Loop trough cells to calculate CoM positions of each cell
        for (unsigned int j = 0; j < ncside; j++)
        {
            for (unsigned int k = 0; k < ncside; k++)
            {
                cell[ncside*j+k].m = cell_aux[ncside*j+k].m;

                if (cell_aux[ncside*j+k].m){
                    
                    cell[ncside*j+k].x = cell_aux[ncside*j+k].x/cell_aux[ncside*j+k].m;
                    cell[ncside*j+k].y = cell_aux[ncside*j+k].y/cell_aux[ncside*j+k].m;
                }
                /*else
                {
                    cell[j][k].x = (j + 0.5) / ncside;
                    cell[j][k].y = (k + 0.5) / ncside;
                }*/
            }
        }// end loop to update cell
        free(cell_aux);
    }

    // Declaration of global mass info
    double total_mass = 0.0, TotalCenter_x = 0.0, TotalCenter_y = 0.0;

    // Update global CoM and total mass
    update_global_quantities(ncside, TotalCenter_x, TotalCenter_y, total_mass, cell);

    // Print required results
    cout << par[0].x << " " << par[0].y << endl;
    cout << TotalCenter_x << " " << TotalCenter_y << endl;
    free(par);
    free(cell);
    time_req = clock()- time_req;
    cout << "It took " << (float)time_req/CLOCKS_PER_SEC << " seconds" << endl;

    return 0;
}
