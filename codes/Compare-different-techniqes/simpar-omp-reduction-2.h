#include <iostream>
#include <stdlib.h>
#include <omp.h>
#include "funcDef.h"

using namespace std;

void reduction_array(long seed,unsigned int ncside,size_t n_part, unsigned int ntstep)
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
        double *cell_info = (double *)calloc(n_cell*3 , sizeof(double));

        size_t i;
        #pragma omp parallel private(i)
        {
            //update initiall cells            
            #pragma omp for reduction(+:cell_info[:n_cell*3]) 
            for (i = 0; i < n_part; i++) //loop over the particles
            {
                
                //---------------------------------------------------------------
                unsigned int c_i = par[i].x * ncside;
                unsigned int c_j = par[i].y * ncside;
                //----------------------------------------------------------------
                //update_cell
                cell_info[(c_i * ncside + c_j)*3] += par[i].m;
                cell_info[(c_i * ncside + c_j)*3+1] += par[i].x * par[i].m;
                cell_info[(c_i * ncside + c_j)*3+2] += par[i].y * par[i].m;
            }//end of loop over the particles 
            #pragma omp for //adjust center of mass
            for (i = 0; i < n_cell; i++) // loop to update cell
            {
                if (cell_info[i*3])  // Only consider cells with mass
                {
                    cell_info[i*3+1] = cell_info[i*3+1] / cell_info[i*3];
                    cell_info[i*3+2] = cell_info[i*3+2] / cell_info[i*3];
                }
            }

        } //end of parallel section
//=================================================================================
//================     Initilization done! ========================================
//=================================================================================


//=================================================================================
//================    start Loop over time ========================================
//================================================================================= 
        
        for (unsigned int t_step = 0; t_step < ntstep; t_step++) //loop over time
        {
            // Auxilary matrix containing cells of the problem for the next time step 

            size_t i;
            //==========================================================================================================  
            #pragma omp parallel private(i) 
            {
                //----------------------------------------------------------------
                // Loop over particles
                #pragma omp for  
                for (i = 0; i < n_part; i++)
                {//loop over particles to calculate forces + Update position and velocity  
                    double ax = 0.0, ay = 0.0; // ax,ay acceleration in (x,y) direction
                    unsigned int c_i = par[i].x * ncside;
                    unsigned int c_j = par[i].y * ncside;
                    // Calculate force components
                    calculate_acceleration_array(c_i, c_j, ncside, par[i].x, par[i].y, par[i].m, ax, ay, cell_info); // devide the loops and see what happens

                    // Update particle positions
                    update_velocities_and_positions(ax, ay, par[i]);
                }//end of loop over particles to calculate forces + Update position and velocity  
                #pragma omp for // Loop trough cells to calculate CoM positions of each cell
                for (i = 0; i < 3*n_cell; i++) // reset cell_info 
                {
                    cell_info[i]=0;
                }
                #pragma omp for reduction(+:cell_info[:n_cell*3])
                for (i = 0; i < n_part; i++)// update cells
                {
                    unsigned int c_i = par[i].x * ncside;
                    unsigned int c_j = par[i].y * ncside;
                    // Update particle's cell info
                    //update_cell(cell_aux[c_i * ncside + c_j], par[i].m, par[i].x, par[i].y);
                    cell_info[(c_i * ncside + c_j)*3] += par[i].m;
                    cell_info[(c_i * ncside + c_j)*3+1] += par[i].x * par[i].m;
                    cell_info[(c_i * ncside + c_j)*3+2] += par[i].y * par[i].m;
                } //end of update cells
                // end of loop over particles
                //==========================================================================================================
                #pragma omp for // Loop trough cells to calculate CoM positions of each cell
                for (i = 0; i < n_cell; i++) // adjust center of mass
                {
                    if (cell_info[i*3])
                    {
                        cell_info[i*3+1]= cell_info[i*3+1] / cell_info[i*3];
                        cell_info[i*3+2]= cell_info[i*3+2] / cell_info[i*3];
                    }
                } // end of adjust center of mass
            } //end parallel section
        //==========================================================================================================  
        } // end loop over time
//=================================================================================
//================    end Loop over time ==========================================
//=================================================================================        
    
    // Declaration of overall mass info
    double total_mass = 0.0, TotalCenter_x = 0.0, TotalCenter_y = 0.0;
    
    // Update overall CoM and total mass
    #pragma omp parallel for if (n_cell > 100) reduction(+ \
                                                    : TotalCenter_x, TotalCenter_y, total_mass)
        for (size_t i = 0; i < n_cell; i++)
        {
            // Calculate info of each cell
            TotalCenter_x += cell_info[i*3+1]*cell_info[i*3];
            TotalCenter_y += cell_info[i*3+2] *cell_info[i*3];
            total_mass += cell_info[i*3];
        }
        // Update positions
        TotalCenter_x /= total_mass;
        TotalCenter_y /= total_mass;
        // Print required results
        printf("%.2f %.2f \n",par[0].x,par[0].y);
        printf("%.2f %.2f \n",TotalCenter_x,TotalCenter_y);
        free(cell_info);
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
