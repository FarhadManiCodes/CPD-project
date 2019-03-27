#include <iostream>
#include <stdlib.h>
#include <ctime>
#include <omp.h>
#include <chrono>
#include "funcDef.h"

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
    unsigned int CPU_Cache_line_size = 64;
    //-------------------------------------------------------------------------------------
    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();

    // Declarations
    particle_t *par = (particle_t *)aligned_alloc(CPU_Cache_line_size, n_part * sizeof(particle_t)); // vector containing all particles of the problem

    // Initialize particles and cells
    init_particles(seed, ncside, n_part, par);
    if (2 * ncside < 1.41424 / EPSLON)
    {
        size_t i;
        size_t n_cell = ncside * ncside;
        cell_t *cell = (cell_t *)calloc(n_cell, sizeof(cell_t));

        #pragma omp parallel private(i)
        {
            
        #pragma omp for
            for (i = 0; i < n_part; i++)
            {
                //---------------------------------------------------------------
                unsigned int c_i = par[i].x * ncside;
                unsigned int c_j = par[i].y * ncside;
                //----------------------------------------------------------------
                //================================================================
                #pragma omp atomic
                cell[c_i * ncside + c_j].m += par[i].m; // better use of atomic
                #pragma omp atomic
                cell[c_i * ncside + c_j].x += par[i].m * par[i].x;
                #pragma omp atomic
                cell[c_i * ncside + c_j].y += par[i].m * par[i].y;
            }
        #pragma omp  for
        for (i = 0; i < n_cell; i++)
        {
            if (cell[i].m) // Only consider cells with mass greater then eps
            {
                // Update cell center of mass positions using the total mass of the cell
                cell[i].x /= cell[i].m;
                cell[i].y /= cell[i].m;
            }
        }
        } //end of parallel section

        //Initilization done!
//--------------------------------------------------------------------------------------------
        // Loop over time
        for (unsigned int t_step = 0; t_step < ntstep; t_step++)
        {

            cell_t *cell_aux = (cell_t *)calloc(n_cell, sizeof(cell_t)); // Auxilary matrix containing cells of the problem for the next time step
            size_t i;
            #pragma omp parallel private(i)
            {
            // Loop over particles
            #pragma omp  for
            for ( i = 0; i < n_part; i++)
            {
                double ax = 0.0, ay = 0.0; // ax,ay acceleration in (x,y) direction
                unsigned int c_i = par[i].x * ncside;
                unsigned int c_j = par[i].y * ncside;
                // Calculate force components
                calculate_acceleration(c_i, c_j, ncside, par[i].x, par[i].y, par[i].m, ax, ay, cell); // devide the loops and see what happens

                // Update particle positions
                update_velocities_and_positions(ax, ay, par[i]);
                c_i = par[i].x * ncside;
                c_j = par[i].y * ncside;
                // Update particle's cell info
                //locate_and_update_cell_info(par[i].x,par[i].y,par[i].m, cell_aux[c_i*ncside+c_j]);

                #pragma omp atomic
                cell_aux[c_i * ncside + c_j].m += par[i].m; // better use of atomic
                #pragma omp atomic
                cell_aux[c_i * ncside + c_j].x += par[i].m * par[i].x;
                #pragma omp atomic
                cell_aux[c_i * ncside + c_j].y += par[i].m * par[i].y;
                //update_cell(cell_aux[c_i*ncside+c_j],par[i].m,par[i].x,par[i].y);

            } // end of loop over particles
        // Loop trough cells to calculate CoM positions of each cell
            #pragma omp for simd
            for ( i = 0; i < n_cell; i++)
            {
                cell[i].m = cell_aux[i].m;
            }
            #pragma omp for
            for ( i = 0; i < n_cell; i++)
            {
                if (cell_aux[i].m)
                {
                    cell[i].x = cell_aux[i].x / cell_aux[i].m;
                    cell[i].y = cell_aux[i].y / cell_aux[i].m;
                }

            } // end loop to update cell
            }// end of parallel section
            free(cell_aux);
        }

        // Declaration of global mass info

        double total_mass = 0.0, TotalCenter_x = 0.0, TotalCenter_y = 0.0;

// Update global CoM and total mass
#pragma omp parallel for reduction(+ \
                                   : TotalCenter_x, TotalCenter_y, total_mass)
        for (size_t j = 0; j < n_cell; j++)
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
        free(cell);
    }
    else
    {

        cout << "Entered there" << endl;
        for (unsigned int t_step = 0; t_step < ntstep; t_step++)
        {
#pragma omp parallel for if (n_part > 100)
            for (size_t i = 0; i < n_part; i++)
            {
                update_velocities_and_positions(0, 0, par[i]);
            }
        }
        double total_mass = 0.0, TotalCenter_x = 0.0, TotalCenter_y = 0.0;
#pragma omp parallel for if (n_part > 100)
        for (size_t j = 0; j < n_part; j++)
        {
            TotalCenter_x += par[j].x * par[j].m;
            TotalCenter_y += par[j].y * par[j].m;
            total_mass += par[j].m;
        }
        // Update positions
        TotalCenter_x /= total_mass;
        TotalCenter_y /= total_mass;

        // Print required results
        cout << par[0].x << " " << par[0].y << endl;
        cout << TotalCenter_x << " " << TotalCenter_y << endl;
    }
    free(par);
    chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
    chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
    float time_req = time_span.count();
    int numberofthreads = omp_get_max_threads();
    cout << "Number of threads  = " << numberofthreads << endl;
    cout << "It took " << time_req << " seconds" << endl;

    return 0;
}
