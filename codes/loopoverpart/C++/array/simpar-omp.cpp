#include <iostream>
#include <stdlib.h>
#include <ctime>
#include <omp.h>
#include "funcDef.h"


using namespace std;

int main(int argc, char **argv)
{
    int numberofthreads = omp_get_max_threads();
    if (numberofthreads > 4)
        numberofthreads = 4;
    omp_set_num_threads (numberofthreads);    
    long seed; 
    unsigned int ncside, ntstep;
    unsigned long n_part;

    // Function argument assignments
    seed = atol(argv[1]);   //  seed for the random number generator
    ncside = atol(argv[2]); //  size of the grid (number of cells on the side)
    n_part = atoll(argv[3]); //  number of particles
    ntstep = atol(argv[4]); //  number of time steps
    int CPU_Cache_line_size = 64;
    clock_t time_req;
    time_req = clock();

    // Declarations
    cell_t *cell = (cell_t *)calloc(ncside*ncside,sizeof(cell_t)); // Matrix containing cells of the problem
    particle_t *par = (particle_t *) aligned_alloc(CPU_Cache_line_size,n_part*sizeof(particle_t)); // vector containing all particles of the problem

    // Initialize particles and cells
    init_particles(seed, ncside, n_part, par, cell);    
       
    // Loop over time
    for (unsigned int t_step = 0; t_step < ntstep; t_step++)
    {

        cell_t *cell_aux = (cell_t *)calloc(ncside*ncside,sizeof(cell_t)); // Auxilary matrix containing cells of the problem for the next time step 
        
        // Loop over particles
        #pragma omp parallel for 
        for (unsigned long i = 0; i < n_part; i++)
        {
            double Fx = 0.0, Fy = 0.0;     // Fx,Fy force in (x,y) direction

            // Calculate force components
            calculate_forces(par[i].c_i,par[i].c_j, ncside, par[i].x, par[i].y,par[i].m, Fx, Fy, cell);

            // Update particle positions
            update_velocities_and_positions(i,ncside, Fx, Fy, par[i]);
            
            // Update particle's cell info
            //locate_and_update_cell_info(par[i].x,par[i].y,par[i].m, cell_aux[par[i].c_i*ncside+par[i].c_j]);
            #pragma omp atomic
                cell_aux[par[i].c_i*ncside+par[i].c_j].m += par[i].m;
            #pragma omp atomic 
                cell_aux[par[i].c_i*ncside+par[i].c_j].x += par[i].m*par[i].x;
            #pragma omp atomic 
                cell_aux[par[i].c_i*ncside+par[i].c_j].y += par[i].m*par[i].y;

        
        }// end of loop over particles
        // Loop trough cells to calculate CoM positions of each cell
        //#pragma omp target map(to: cell_aux) map(tofrom: cell)
        //#pragma omp parallel for simd
        #pragma omp simd
        for (unsigned int j = 0; j < ncside*ncside; j++)
        {          
                cell[j].m = cell_aux[j].m;
                if (cell_aux[j].m){                    
                    cell[j].x = cell_aux[j].x/cell_aux[j].m;
                    cell[j].y = cell_aux[j].y/cell_aux[j].m;
                }
            
        }// end loop to update cell
        //#pragma omp end target
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
    cout << "Number of threads  = " << numberofthreads << endl;
    cout << "It took " << (float)time_req/CLOCKS_PER_SEC/numberofthreads << " seconds" << endl;
    

    return 0;
}
