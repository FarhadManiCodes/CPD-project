#include <iostream>
#include <stdlib.h>
#include <omp.h>
#include "funcDef.h"

using namespace std;

void reduction(long seed,unsigned int ncside,size_t n_part, unsigned int ntstep)
{

    unsigned int CPU_Cache_line_size = 64;
    //-------------------------------------------------------------------------------------

    // Declarations
    particle_t *par = (particle_t *)aligned_alloc(CPU_Cache_line_size, n_part * sizeof(particle_t)); // vector containing all particles of the problem
//=================================================================================
//================     Initilization start ========================================
//=================================================================================
    // Initialize particles and cells
    init_particles(seed, ncside, n_part, par);
    if (2 * ncside <  1.41424 / EPSLON)
    {
        size_t n_cell = ncside * ncside;
        cell_t *cell = (cell_t *)aligned_alloc(CPU_Cache_line_size, n_cell * sizeof(cell_t));
        double *cell_x = (double *)calloc(n_cell , sizeof(double));
        double *cell_y = (double *)calloc(n_cell , sizeof(double));
        double *cell_m = (double *)calloc(n_cell , sizeof(double));
        size_t i;
        #pragma omp parallel private(i)
        {
            //update initiall cells            
            #pragma omp for reduction(+:cell_x[:n_cell],cell_y[:n_cell],cell_m[:n_cell]) 
            for (i = 0; i < n_part; i++) //loop over the particles
            {
                
                //---------------------------------------------------------------
                unsigned int c_i = par[i].x * ncside;
                unsigned int c_j = par[i].y * ncside;
                //----------------------------------------------------------------
                //update_cell
                cell_m[c_i * ncside + c_j] += par[i].m;
                cell_x[c_i * ncside + c_j] += par[i].x * par[i].m;
                cell_y[c_i * ncside + c_j] += par[i].y * par[i].m;
            }//end of loop over the particles 
            #pragma omp for //adjust center of mass
            for (i = 0; i < n_cell; i++) // loop to update cell
            {
                cell[i].m = cell_m[i];
                if (cell_m[i])  // Only consider cells with mass
                {
                    cell[i].x = cell_x[i] / cell_m[i];
                    cell[i].y = cell_y[i] / cell_m[i];
                }
            }

        } //end of parallel section
        free(cell_x);
        free(cell_y);
        free(cell_m);
//=================================================================================
//================     Initilization done! ========================================
//=================================================================================


//=================================================================================
//================    start Loop over time ========================================
//================================================================================= 
        
        for (unsigned int t_step = 0; t_step < ntstep; t_step++) //loop over time
        {
            
            double *cell_x = (double *)calloc(n_cell , sizeof(double));
            double *cell_y = (double *)calloc(n_cell , sizeof(double)); // Auxilary matrix containing cells of the problem for the next time step 
            double *cell_m = (double *)calloc(n_cell , sizeof(double));
            size_t i;
            //==========================================================================================================  
            #pragma omp parallel private(i) 
            {
                //----------------------------------------------------------------
                // Loop over particles
                #pragma omp for nowait 
                for (i = 0; i < n_part; i++)
                {//loop over particles to calculate forces + Update position and velocity  
                    double ax = 0.0, ay = 0.0; // ax,ay acceleration in (x,y) direction
                    unsigned int c_i = par[i].x * ncside;
                    unsigned int c_j = par[i].y * ncside;
                    // Calculate force components
                    calculate_acceleration(c_i, c_j, ncside, par[i].x, par[i].y, par[i].m, ax, ay, cell); // devide the loops and see what happens

                    // Update particle positions
                    update_velocities_and_positions(ax, ay, par[i]);
                }//end of loop over particles to calculate forces + Update position and velocity  
                #pragma omp for reduction(+: cell_x[:n_cell], cell_y[:n_cell], cell_m[:n_cell])
                for (i = 0; i < n_part; i++)// update cells
                {
                    unsigned int c_i = par[i].x * ncside;
                    unsigned int c_j = par[i].y * ncside;
                    // Update particle's cell info
                    //update_cell(cell_aux[c_i * ncside + c_j], par[i].m, par[i].x, par[i].y);
                    cell_m[c_i * ncside + c_j] += par[i].m;
                    cell_x[c_i * ncside + c_j] += par[i].x * par[i].m;
                    cell_y[c_i * ncside + c_j] += par[i].y * par[i].m;
                } //end of update cells
                // end of loop over particles
                //==========================================================================================================
                #pragma omp for // Loop trough cells to calculate CoM positions of each cell
                for (i = 0; i < n_cell; i++) // adjust center of mass
                {
                    cell[i].m = cell_m[i];
                    if (cell_m[i])
                    {
                        cell[i].x = cell_x[i] / cell_m[i];
                        cell[i].y = cell_y[i] / cell_m[i];
                    }
                } // end of adjust center of mass
            } //end parallel section
        //==========================================================================================================  
            free(cell_x);
            free(cell_y);
            free(cell_m);
        } // end loop over time
//=================================================================================
//================    end Loop over time ==========================================
//=================================================================================        
    
    // Declaration of overall mass info
    double total_mass = 0.0, TotalCenter_x = 0.0, TotalCenter_y = 0.0;
    
    // Update overall CoM and total mass
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
        printf("%.2f %.2f \n",par[0].x,par[0].y);
        printf("%.2f %.2f \n",TotalCenter_x,TotalCenter_y);
        free(cell);
    }
    else
    {
        for (unsigned int t_step = 0; t_step < ntstep; t_step++)
        {
            #pragma omp parallel for if (n_part > 16)
            for (size_t i = 0; i < n_part; i++)
            {
                update_velocities_and_positions(0, 0, par[i]);
            }
        }
        double total_mass = 0.0, TotalCenter_x = 0.0, TotalCenter_y = 0.0;
        #pragma omp parallel for if (n_part > 16)
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
        printf("%.2f %.2f \n",par[0].x,par[0].y);
        printf("%.2f %.2f \n",TotalCenter_x,TotalCenter_y);
    }
    free(par);
}
