#include <iostream>
#include <stdlib.h>
#include <ctime>
#include <chrono>
#include <omp.h>
#include <fstream>
#include "funcDef.h"
#include "simpar.h"
#include "simpar-omp.h"
#include "simpar-omp-reduction.h"
#include "simpar-omp-atomic.h"
#include <ratio>


using namespace std;

int main(){
    ofstream myfile;
    int totalseed = 2;
    int totalncside = 7;
    int totaln_part = 8;
    int totalntstep = 1;
    myfile.open ("time.txt", ios::trunc);
for (unsigned int k = 2; k  < totaln_part ;k++){
    for (unsigned int i = 0; i  < totalntstep ;i++){
        for (unsigned int j = 2; j  <= totalncside ;j++){
            
                float av_serial_time = 0.0;
                float av_omp1_time = 0.0;
                float av_omp2_time = 0.0;
                float av_omp4_time = 0.0;
                //---------------------------------------------------------------------------              
                unsigned long n_part = 1*pow(10,k+1);
                unsigned int ncside = pow(2,j);
                unsigned int ntstep = 10;
                for (unsigned int s = 0; s < totalseed  ;s++){
                        long seed = s; 
                        /*unsigned int ncside = j*5;
                        unsigned long n_part = 10000*pow(2,k);
                        unsigned int ntstep = i*10; */
                        //test serial
                        chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
                        serial(seed, ncside,n_part,ntstep);
                        chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
                        chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
                        float time_serial = time_span.count();
 
                        //test omp
                        t1 = chrono::high_resolution_clock::now();
                        omp_set_num_threads(4);
                        manual(seed, ncside,n_part,ntstep);
                        t2 = chrono::high_resolution_clock::now();
                        time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
                        float time_omp1 = time_span.count();
                        //---------------------------------------------------------------------------
                        //test omp
                        t1 = chrono::high_resolution_clock::now();
                        omp_set_num_threads(4);
                        reduction(seed, ncside,n_part,ntstep);
                        t2 = chrono::high_resolution_clock::now();
                        time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
                        float time_omp2 = time_span.count();
                        //---------------------------------------------------------------------------
                        //test omp
                        t1 = chrono::high_resolution_clock::now();
                        omp_set_num_threads(4);
                        atomic(seed, ncside,n_part,ntstep);
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
    return 0;
}

