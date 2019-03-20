#include <iostream>
#include <stdlib.h>
#include <ctime>
#include <chrono>
#include "funcDef.h"


using namespace std;

int main(int argc, char **argv)
{
    long seed; 
    size_t ncside, ntstep;
    size_t n_part;

    // Function argument assignments
    seed = atol(argv[1]);   //  seed for the random number generator
    ncside = atol(argv[2]); //  size of the grid (number of cells on the side)
    n_part = atoll(argv[3]); //  number of particles
    ntstep = atol(argv[4]); //  number of time steps

    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();

    // Declarations
    cell_t *cell = (cell_t *)calloc(ncside*ncside,sizeof(cell_t)); // Matrix containing cells of the problem
    particle_t *par =  (particle_t *) malloc(n_part * sizeof(particle_t)); // vector containing all particles of the problem

    // Initialize particles and cells
    init_particles(seed, ncside, n_part, par, cell);    
       
    // Loop over time
    for (size_t t_step = 0; t_step < ntstep; t_step++)
    {

        cell_t *cell_aux = (cell_t *)calloc(ncside*ncside,sizeof(cell_t)); // Auxilary matrix containing cells of the problem for the next time step 

        // Loop over particles
        for (size_t i = 0; i < n_part; i++)
        {
            double Fx = 0.0, Fy = 0.0;     // Fx,Fy force in (x,y) direction

            // Calculate force components
            calculate_forces(par[i].c_i,par[i].c_j, ncside, par[i].x, par[i].y,par[i].m, Fx, Fy, cell);

            // Update particle positions
            update_velocities_and_positions(i,ncside, Fx, Fy, par[i]);
            
            // Update particle's cell info
            //locate_and_update_cell_info(par[i].x,par[i].y,par[i].m, cell_aux[par[i].c_i*ncside+par[i].c_j]);
             cell_aux[par[i].c_i*ncside+par[i].c_j].m += par[i].m;
             cell_aux[par[i].c_i*ncside+par[i].c_j].x += par[i].m*par[i].x;
             cell_aux[par[i].c_i*ncside+par[i].c_j].y += par[i].m*par[i].y;

        
        }// end of loop over particles
        // Loop trough cells to calculate CoM positions of each cell
        for (size_t j = 0; j < ncside*ncside; j++)
        {          
                cell[j].m = cell_aux[j].m;
                if (cell_aux[j].m){                    
                    cell[j].x = cell_aux[j].x/cell_aux[j].m;
                    cell[j].y = cell_aux[j].y/cell_aux[j].m;
                }
            
        }// end loop to update cell
        free(cell_aux);
    }

    // Declaration of global mass info
    double total_mass = 0.0, TotalCenter_x = 0.0, TotalCenter_y = 0.0;

    // Update global CoM and total mass
    for (size_t j = 0; j < ncside*ncside; j++)
    {
    // Calculate info of each cell
        TotalCenter_x += cell[j].x * cell[j].m;
        TotalCenter_y += cell[j].y * cell[j].m;
        total_mass += cell[j].m;
    }
        // Update positions
    TotalCenter_x /= total_mass;
    TotalCenter_y /= total_mass;

    // Print required results
    cout << par[0].x << " " << par[0].y << endl;
    cout << TotalCenter_x << " " << TotalCenter_y << endl;
    free(par);
    free(cell);
    chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
    chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
    float time_req = time_span.count(); 
    cout << "It took " << time_req << " seconds" << endl;

    return 0;
}
