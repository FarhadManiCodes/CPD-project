#include <iostream>
#include <stdlib.h>
#include <ctime>
#include <chrono>
#include <omp.h>
#include <fstream>
#include "funcDef.h"
#include <ratio>


using namespace std;

void testserial(long seed, unsigned int ncside,unsigned long n_part,unsigned int ntstep)
{
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
        for (unsigned int j = 0; j < ncside*ncside; j++)
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
    update_global_quantities(ncside, TotalCenter_x, TotalCenter_y, total_mass, cell);

    // Print required results
    //cout << par[0].x << " " << par[0].y << endl;
    //cout << TotalCenter_x << " " << TotalCenter_y << endl;
    free(par);
    free(cell);
}

void test_omp(long seed, unsigned int ncside,unsigned long n_part,unsigned int ntstep)
{
    int CPU_Cache_line_size = 64;

    // Declarations
    cell_t *cell = (cell_t *)calloc(ncside*ncside,sizeof(cell_t)); // Matrix containing cells of the problem
    particle_t *par = (particle_t *) aligned_alloc(CPU_Cache_line_size,n_part*sizeof(particle_t)); // vector containing all particles of the problem

    // Initialize particles and cells
    init_particles(seed, ncside, n_part, par, cell);    
       
    // Loop over time
    bool useopenmp_particle = (n_part > 100 && n_part/ncside > 10);
    for (unsigned int t_step = 0; t_step < ntstep; t_step++)
    {

        cell_t *cell_aux = (cell_t *)calloc(ncside*ncside,sizeof(cell_t)); // Auxilary matrix containing cells of the problem for the next time step 
        
        // Loop over particles
        #pragma omp parallel for //if(useopenmp_particle)
        for (unsigned long i = 0; i < n_part; i++)
        {
            double Fx = 0.0, Fy = 0.0;     // Fx,Fy force in (x,y) direction

            // Calculate force components
            calculate_forces(par[i].c_i,par[i].c_j, ncside, par[i].x, par[i].y,par[i].m, Fx, Fy, cell);  // devide the loops and see what happens

            // Update particle positions
            update_velocities_and_positions(i,ncside, Fx, Fy, par[i]);

            // Update particle's cell info
            //locate_and_update_cell_info(par[i].x,par[i].y,par[i].m, cell_aux[par[i].c_i*ncside+par[i].c_j]);
            #pragma omp atomic
                cell_aux[par[i].c_i*ncside+par[i].c_j].m += par[i].m;  // better use of atomic
            #pragma omp atomic 
                cell_aux[par[i].c_i*ncside+par[i].c_j].x += par[i].m*par[i].x;
            #pragma omp atomic 
                cell_aux[par[i].c_i*ncside+par[i].c_j].y += par[i].m*par[i].y;

        
        }// end of loop over particles
        // Loop trough cells to calculate CoM positions of each cell
        //#pragma omp parallel for 
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
    //cout << par[0].x << " " << par[0].y << endl;
    //cout << TotalCenter_x << " " << TotalCenter_y << endl;
    free(par);
    free(cell);
}

int main(){
    ofstream myfile;
    int totalseed = 5;
    int totalncside = 20;
    int totaln_part = 4;
    int totalntstep = 2;
    myfile.open ("time.txt", ios::trunc);
for (unsigned int k = 1; k  <= totaln_part ;k++){
    for (unsigned int i = 1; i  < totalntstep ;i++){
        for (unsigned int j = 2; j  <= totalncside ;j++){
            
                float av_serial_time = 0.0;
                float av_omp1_time = 0.0;
                float av_omp2_time = 0.0;
                float av_omp4_time = 0.0;
                //---------------------------------------------------------------------------              
                unsigned long n_part = 20000*pow(10,k);
                unsigned int ncside = n_part/(j*10000);
                unsigned int ntstep = 10;
                for (unsigned int s = 0; s < totalseed  ;s++){
                        long seed = s; 
                        /*unsigned int ncside = j*5;
                        unsigned long n_part = 10000*pow(2,k);
                        unsigned int ntstep = i*10; */
                        //test serial
                        chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
                        testserial(seed, ncside,n_part,ntstep);
                        chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
                        chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
                        float time_serial = time_span.count();
 
                        //test omp
                        t1 = chrono::high_resolution_clock::now();
                        //omp_set_num_threads(1);
                        //test_omp(seed, ncside,n_part,ntstep);
                        t2 = chrono::high_resolution_clock::now();
                        time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
                        float time_omp1 = time_span.count();
                        //---------------------------------------------------------------------------
                        //test omp
                        t1 = chrono::high_resolution_clock::now();
                        omp_set_num_threads(2);
                        //test_omp(seed, ncside,n_part,ntstep);
                        t2 = chrono::high_resolution_clock::now();
                        time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
                        float time_omp2 = time_span.count();
                        //---------------------------------------------------------------------------
                        //test omp
                        t1 = chrono::high_resolution_clock::now();
                        omp_set_num_threads(4);
                        test_omp(seed, ncside,n_part,ntstep);
                        t2 = chrono::high_resolution_clock::now();
                        time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
                        float time_omp4 = time_span.count();
                        //---------------------------------------------------------------------------
                        //--------------------------------------------------
                            av_serial_time +=time_serial;
                            av_omp1_time+=time_omp1;
                            av_omp2_time+=time_omp2;
                            av_omp4_time+=time_omp4;

                        
                       // myfile << seed << "," << ntstep << "," << ncside << "," << n_part << "," << time_serial << ","<<time_omp1 << ","<<time_omp2 << ","<< time_omp4 << endl;
                        
                        //cout << seed << "," << ncside << "," << n_part << "," << ntstep << "," << time_serial << ","<<time_omp1 << ","<<time_omp2 << ","<< time_omp4 << "\n"<< endl;
                        

                }
                av_serial_time/=totalseed;
                av_omp1_time/=totalseed;
                av_omp2_time/=totalseed;
                av_omp4_time/=totalseed;
                myfile << ntstep << "," << ncside << "," << n_part << "," << av_serial_time << ","<<av_omp1_time << ","<<av_omp2_time << ","<< av_omp4_time << endl;
                
            }   

       }

    }
    myfile.close();
}

