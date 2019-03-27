#include <iostream>
#include <stdlib.h>
#include <omp.h>
#include "funcDef.h"

using namespace std;

void manual(long seed,unsigned int ncside,size_t n_part, unsigned int ntstep)
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
    if (2*ncside < 1.414214 / EPSLON)
    {
        size_t n_cell = ncside * ncside;
        cell_t *cell = (cell_t *)calloc(n_cell, sizeof(cell_t));
        size_t i;
        #pragma omp parallel private(i)
        {
            //update initiall cells
            cell_t *cell_aux = (cell_t *)calloc(n_cell, sizeof(cell_t)); // Auxilary matrix containing cells of the problem for the next time step
            #pragma omp for  //loop over the particles
            for (i = 0; i < n_part; i++)
            {
                //---------------------------------------------------------------
                unsigned int c_i = par[i].x * ncside;
                unsigned int c_j = par[i].y * ncside;
                //----------------------------------------------------------------
                //================================================================

                update_cell(cell_aux[c_i * ncside + c_j], par[i].m, par[i].x, par[i].y);
            } //end of loop over the particles 
#pragma omp critical
            {
                for (i = 0; i < n_cell; i++) //loop over cell for manual reduction
                {
                    cell[i].m += cell_aux[i].m;
                    cell[i].x += cell_aux[i].x;
                    cell[i].y += cell_aux[i].y;
                }
            } //end of critical section 
#pragma omp barrier//end of update initiall cells 
#pragma omp for            
        for (i = 0; i < n_cell; i++) //adjust center of mass
        {
            if (cell[i].m) // Only consider cells with mass
            {
                // Update cell center of mass positions using the total mass of the cell
                cell[i].x /= cell[i].m;
                cell[i].y /= cell[i].m;
            }
        }//end of adjust center of mass
            free(cell_aux);
        } //end of parallel section

//=================================================================================
//================     Initilization done! ========================================
//=================================================================================


//=================================================================================
//================    start Loop over time ========================================
//=================================================================================        

        for (unsigned int t_step = 0; t_step < ntstep; t_step++) //loop over time
        {
            size_t i;
            #pragma omp parallel private(i) 
                        {
            cell_t *cell_aux = (cell_t *)calloc(n_cell, sizeof(cell_t)); // Auxilary matrix containing cells of the problem for the next time step                
            //====================================================================================================================================
            // Loop over particles
                #pragma omp for nowait
                for (i = 0; i < n_part; i++) // Loop over particles to calculate forces + Update position and velocity
                {
                    double ax = 0.0, ay = 0.0; // ax,ay acceleration in (x,y) direction
                    unsigned int c_i = par[i].x * ncside;
                    unsigned int c_j = par[i].y * ncside;
                    // Calculate force components
                    calculate_acceleration(c_i, c_j, ncside, par[i].x, par[i].y, par[i].m, ax, ay, cell); // devide the loops and see what happens

                    // Update particle positions
                    update_velocities_and_positions(ax, ay, par[i]);
                }//end of loop over particles to calculate forces + Update position and velocity  
                // update cells              
                #pragma omp for nowait
                for (i = 0; i < n_part; i++) // update cells
                {
                    unsigned int c_i = par[i].x * ncside;
                    unsigned int c_j = par[i].y * ncside;
                    // Update particle's cell info
                    update_cell(cell_aux[c_i * ncside + c_j], par[i].m, par[i].x, par[i].y);
                } //end of update cells
            // end of loop over particles
            //==========================================================================================================
            //manual reduction
            #pragma omp for 
                for (i = 0; i < n_cell; i++) // reset the cells
                {
                    cell[i].m = 0.0;
                    cell[i].x = 0.0;
                    cell[i].y = 0.0;
                } // end of reseting
            #pragma omp critical
                {
                    for (i = 0; i < n_cell; i++) // loop to update cell
                    {
                        cell[i].m += cell_aux[i].m;
                        cell[i].x += cell_aux[i].x;
                        cell[i].y += cell_aux[i].y;
                    }
                } //end of critical section
#pragma omp barrier  //manual reduction //end of update cells
#pragma omp for 
            for (i = 0; i < n_cell; i++) // adjust center of mass
            { // Loop trough cells to calculate CoM positions of each cell
                if (cell[i].m)
                {
                    cell[i].x /= cell[i].m;
                    cell[i].y /= cell[i].m;
                }
            } //end of adjust center of mass
                free(cell_aux);
            } //end parallel section
        //======================================================================================================
        }// end loop over time
//=================================================================================
//================    end Loop over time ==========================================
//=================================================================================  

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
                update_velocities_and_positions( 0, 0, par[i]);
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
