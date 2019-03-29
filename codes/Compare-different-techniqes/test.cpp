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

int main()
{
        ofstream myfile;
        int totalseed = 5;
        int totalncside = 8;
        int totaln_part = 24;
        int totalntstep = 1;
        myfile.open("time.txt", ios::trunc);
        printf("Test started\n");
        for (unsigned int k = 0; k <= totaln_part; k++)
        {
                printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
                for (unsigned int i = 0; i < totalntstep; i++)
                {
                        printf("========================================================================\n");
                        for (unsigned int j = 0; j <= totalncside; j++)
                        {

                                float av_serial_time = 0.0;
                                float av_omp1_time = 0.0;
                                float av_omp2_time = 0.0;
                                float av_omp4_time = 0.0;
                                //---------------------------------------------------------------------------
                                size_t n_part = 10 * pow(2, k);
                                unsigned int ncside = 4 * pow(2, j);
                                unsigned int ntstep = 20;
                                printf("----------------------------------------------------------\n");
                                for (unsigned int s = 0; s < totalseed; s++)
                                {
                                        long seed = s;
                                        printf("seed = %ld, ncside = %u, n_part = %lu\n", seed, ncside, n_part);
                                        //test serial
                                        printf("serial\n");
                                        chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
                                        serial(seed, ncside, n_part, ntstep);
                                        chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
                                        chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
                                        float time_serial = time_span.count();

                                        //test omp
                                        //printf("manual\n");
                                        omp_set_num_threads(4);
                                        t1 = chrono::high_resolution_clock::now();
                                        //manual(seed, ncside, n_part, ntstep);
                                        t2 = chrono::high_resolution_clock::now();
                                        time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
                                        float time_omp1 = time_span.count();
                                        //---------------------------------------------------------------------------
                                        //test omp
                                        //printf("reduction\n");
                                        t1 = chrono::high_resolution_clock::now();
                                        //omp_set_num_threads(4);
                                        //reduction(seed, ncside, n_part, ntstep);
                                        t2 = chrono::high_resolution_clock::now();
                                        time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
                                        float time_omp2 = time_span.count();
                                        //---------------------------------------------------------------------------
                                        //test omp
                                        printf("test omp\n");
                                        t1 = chrono::high_resolution_clock::now();
                                        if (n_part / (ncside * ncside) > 5)
                                        {
                                                if (ncside < 150)
                                                        reduction(seed, ncside, n_part, ntstep);

                                                else
                                                        manual(seed, ncside, n_part, ntstep);
                                        }
                                        else
                                        {
                                                atomic(seed, ncside, n_part, ntstep);
                                        }
                                        t2 = chrono::high_resolution_clock::now();
                                        time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
                                        float time_omp4 = time_span.count();
                                        //---------------------------------------------------------------------------
                                        //--------------------------------------------------
                                        av_serial_time += time_serial;
                                        av_omp1_time += time_omp1;
                                        av_omp2_time += time_omp2;
                                        av_omp4_time += time_omp4;
                                }
                                av_serial_time /= totalseed;
                                av_omp1_time /= totalseed;
                                av_omp2_time /= totalseed;
                                av_omp4_time /= totalseed;
                                myfile << ntstep << "," << ncside << "," << n_part << "," << av_serial_time << "," << av_omp1_time << "," << av_omp2_time << "," << av_omp4_time << endl;
                        }
                }
        }
        myfile.close();
        return 0;
}
