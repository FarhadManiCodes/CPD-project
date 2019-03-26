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
    seed = atol(argv[1]);   //  seed for the random number generator
    ncside = atol(argv[2]); //  size of the grid (number of cells on the side)
    n_part = atol(argv[3]); //  number of particles
    ntstep = atol(argv[4]); //  number of time steps
    unsigned int CPU_Cache_line_size = 64;
    //-------------------------------------------------------------------------------------
    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();

    // Declarations
    particle_t *par = (particle_t *)aligned_alloc(CPU_Cache_line_size, n_part * sizeof(particle_t)); // vector containing all particles of the problem

    // Initialize particles and cells
    init_particles(seed, ncside, n_part, par);
    // Loop trough cells to calculate CoM positions of each cell
    bool useopenmp_particle = (n_part > 100 && n_part / ncside > 10);
    if (2 * ncside <  1.41424 / EPSLON)
    {
        size_t i, j;
        size_t n_cell = ncside * ncside;
        cell_t *cell = (cell_t *)aligned_alloc(CPU_Cache_line_size, n_cell * sizeof(cell_t));
        double *cell_x = (double *)calloc(n_cell , sizeof(double));
        double *cell_y = (double *)calloc(n_cell , sizeof(double));
        double *cell_m = (double *)calloc(n_cell , sizeof(double));
#pragma omp parallel private(i, j) if (useopenmp_particle)
        {
#pragma omp for reduction(+                         \
                          : cell_x[:n_cell], cell_y \
                                  [:n_cell], cell_m \
                                  [:n_cell])
            for (i = 0; i < n_part; i++)
            {
                //---------------------------------------------------------------
                unsigned int c_i = par[i].x * ncside;
                unsigned int c_j = par[i].y * ncside;
                //----------------------------------------------------------------
                cell_m[c_i * ncside + c_j] += par[i].m;
                cell_x[c_i * ncside + c_j] += par[i].x * par[i].m;
                cell_y[c_i * ncside + c_j] += par[i].y * par[i].m;
            }
#pragma omp for
            for (j = 0; j < n_cell; j++) // loop to update cell
            {
                cell[j].m = cell_m[j];
                if (cell_m[j])
                {
                    cell[j].x = cell_x[j] / cell_m[j];
                    cell[j].y = cell_y[j] / cell_m[j];
                }
            }

        } //end of parallel section
         //Initilization done!
        free(cell_x);
        free(cell_y);
        free(cell_m);

        //----------------------------------------------------------------------------------------------------------------------
        //------------- starrt loop over time ----------------------------------------------------------------------------------
        //----------------------------------------------------------------------------------------------------------------------
        // start Loop over time
        for (unsigned int t_step = 0; t_step < ntstep; t_step++)
        {
            size_t i, j;
        double *cell_x = (double *)calloc(n_cell , sizeof(double));
        double *cell_y = (double *)calloc(n_cell , sizeof(double));
        double *cell_m = (double *)calloc(n_cell , sizeof(double));
//==========================================================================================================      firstprivate(cell_aux)
#pragma omp parallel private(i, j) if (useopenmp_particle) //if(n_part > 100 && n_part/ncside > 10)
            {
                //============================================================================================================
                //============================================================================================================
                // Loop over particles

#pragma omp for nowait
                for (i = 0; i < n_part; i++)
                {
                    double ax = 0.0, ay = 0.0; // ax,ay acceleration in (x,y) direction
                    unsigned int c_i = par[i].x * ncside;
                    unsigned int c_j = par[i].y * ncside;
                    // Calculate force components
                    calculate_acceleration(c_i, c_j, ncside, par[i].x, par[i].y, par[i].m, ax, ay, cell); // devide the loops and see what happens

                    // Update particle positions
                    update_velocities_and_positions(ax, ay, par[i]);
                }
#pragma omp for reduction(+                         \
                          : cell_x[:n_cell], cell_y \
                                  [:n_cell], cell_m \
                                  [:n_cell])

                for (i = 0; i < n_part; i++)
                {
                    unsigned int c_i = par[i].x * ncside;
                    unsigned int c_j = par[i].y * ncside;
                    // Update particle's cell info
                    //update_cell(cell_aux[c_i * ncside + c_j], par[i].m, par[i].x, par[i].y);
                    cell_m[c_i * ncside + c_j] += par[i].m;
                    cell_x[c_i * ncside + c_j] += par[i].x * par[i].m;
                    cell_y[c_i * ncside + c_j] += par[i].y * par[i].m;
                } // end of loop over particles
                  //==========================================================================================================
                  // Loop trough cells to calculate CoM positions of each cell
#pragma omp for
                for (j = 0; j < n_cell; j++) // loop to update cell
                {
                    cell[j].m = cell_m[j];
                    if (cell_m[j])
                    {
                        cell[j].x = cell_x[j] / cell_m[j];
                        cell[j].y = cell_y[j] / cell_m[j];
                    }
                }

            } //end parallel section
            free(cell_x);
            free(cell_y);
            free(cell_m);
        } // end loop over time
              //==============================================================================================================================
              // Declaration of global mass info
        double total_mass = 0.0, TotalCenter_x = 0.0, TotalCenter_y = 0.0;

// Update global CoM and total mass
#pragma omp parallel for if (n_cell > 100) reduction(+ \
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
        
        cout <<"Entered there"<<endl;
        for (unsigned int t_step = 0; t_step < ntstep; t_step++)
        {
#pragma omp parallel for if (useopenmp_particle)
            for (size_t i = 0; i < n_part; i++)
            {
                update_velocities_and_positions(0, 0, par[i]);
            }
        }
        double total_mass = 0.0, TotalCenter_x = 0.0, TotalCenter_y = 0.0;
#pragma omp parallel for if (useopenmp_particle)
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
