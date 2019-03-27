#include <iostream>
#include <stdlib.h>
#include "funcDef.h"

using namespace std;

void serial(long seed,unsigned int ncside,size_t n_part, unsigned int ntstep)
{

    // Declarations

    particle_t *par = (particle_t *)malloc(n_part * sizeof(particle_t)); // vector containing all particles of the problem
//=================================================================================
//================     Initilization start ========================================
//=================================================================================    
    // Initialize particles 
    init_particles(seed, ncside, n_part, par);

    if (2 * ncside < 1.41424 / EPSLON)
    {
        size_t n_cell = ncside * ncside;
        cell_t *cell = (cell_t *)calloc(n_cell, sizeof(cell_t));    
        for (size_t i = 0; i < n_part; i++) //loop over the particles to update cells
        {
            //---------------------------------------------------------------
            unsigned int c_i = par[i].x * ncside;
            unsigned int c_j = par[i].y * ncside;
            //----------------------------------------------------------------
            //================================================================

            update_cell(cell[c_i * ncside + c_j], par[i].m, par[i].x, par[i].y);
        }
        for (size_t j = 0; j < n_cell; j++) //adjust center of mass
        {
            if (cell[j].m) // Only consider cells with mass 
            {
                // Update cell center of mass positions using the total mass of the cell
                cell[j].x /= cell[j].m;
                cell[j].y /= cell[j].m;
            }
        }//end of adjust center of mass
//=================================================================================
//================     Initilization done! ========================================
//=================================================================================   
 
//=================================================================================
//================    start Loop over time ========================================
//================================================================================= 
       
        for (unsigned int t_step = 0; t_step < ntstep; t_step++)  // Loop over time
        {

            cell_t *cell_aux = (cell_t *)calloc(n_cell, sizeof(cell_t)); // Auxilary matrix containing cells of the problem for the next time step

            for (size_t i = 0; i < n_part; i++)
            {
                double ax = 0.0, ay = 0.0; // Fx,Fy force in (x,y) direction
                unsigned int c_i = par[i].x * ncside;
                unsigned int c_j = par[i].y * ncside;
                // Calculate force components
                calculate_acceleration(c_i, c_j, ncside, par[i].x, par[i].y, par[i].m, ax, ay, cell); // devide the loops and see what happens

                // Update particle positions
                update_velocities_and_positions( ax, ay, par[i]);
                c_i = par[i].x * ncside;
                c_j = par[i].y * ncside;

                update_cell(cell_aux[c_i * ncside + c_j], par[i].m, par[i].x, par[i].y);

            } // end of loop over particles
            //----------------------------------------------------------------
            // Loop trough cells to calculate CoM positions of each cell
            for (size_t j = 0; j < n_cell; j++) //adjust center of mass
            {
                cell[j].m = cell_aux[j].m;
                if (cell_aux[j].m)
                {
                    cell[j].x = cell_aux[j].x / cell_aux[j].m;
                    cell[j].y = cell_aux[j].y / cell_aux[j].m;
                }

            } // end adjust center of mass
            free(cell_aux);
        }// end loop over time
//=================================================================================
//================    end Loop over time ==========================================
//=================================================================================

        // Declaration of global mass info
        double total_mass = 0.0, TotalCenter_x = 0.0, TotalCenter_y = 0.0;

        // Update global CoM and total mass
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
            for (size_t i = 0; i < n_part; i++)
            {
                update_velocities_and_positions( 0, 0, par[i]);
            }
        }
        double total_mass = 0.0, TotalCenter_x = 0.0, TotalCenter_y = 0.0;
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
