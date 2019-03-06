# include <stdio.h>
# include <stdlib.h>
# include <stddef.h>
#include <math.h>
//-------------------------------------------------------------------
#define RND0_1 ((double) random() / ((long long)1<<31))
#define G 6.67408e-11
#define EPSLON 0.01
#define EPSLON2 0.0001
#define EPS 1e-15
//-----------------particle sturucture definnition -------------------------------
typedef struct particle_t particle_t;
struct particle_t {
  double x;     // x possition of the particle [0,1]
  double y;     // y possition of the particle [0,1]
  double vx;    // velocity of the particle in x direction
  double vy;    // velocity of the particle in y direction
  double m;     // mass of the particle
  unsigned int c_i; // particle's cell index 1 (row index)
  unsigned int c_j; // particle's cell index 2 (column index)
};

/*--------------------------------------------------------------------
 * Function:    init_particles
 * Purpose:     initialize the problem by randomly produce particles,
 * In arg:      seed for the random number generator
 *              size of the grid (number of cells on the side)
 *              number of particles
 *              pointer to struct particle
 *              Cx,Cy = matrix of the center of mass
 *              total_mass_C totall mass of each cell
 */
double init_particles(long seed, long ncside, long long n_part, particle_t *par, double Cx[][ncside], double Cy[][ncside],double mass_C[][ncside],double *CTX,double *CTY){

    long long i;

    srandom(seed);
    double MC_local[ncside][ncside];
    MC_local[ncside][ncside] = 0.0;
    double CX_local[ncside][ncside];
    CX_local[ncside][ncside] = 0.0;
    double CY_local[ncside][ncside],CTX_local=0.0,CTY_local=0.0,total_mass=0;
    CY_local[ncside][ncside] = 0.0;
    for(i = 0; i < n_part; i++) {
        par[i].x = RND0_1;
        par[i].y = RND0_1;
        par[i].vx = RND0_1 / ncside / 10.0;
        par[i].vy = RND0_1 / ncside / 10.0;
        //--------------------------------------------------------------------
        par[i].m = RND0_1 * ncside / (G * 1e6 * n_part);
        //======================== added code ===================================

        unsigned int ci = par[i].x*ncside;
        unsigned int cj = par[i].y*ncside;
        //printf("particle = %lld cell = %u,%u\tx = %0.5f, y = %0.5f\n",i,ci,cj,par[i].x,par[i].y);
        //--------------------------------------------------------------------

        MC_local[ci][cj]+=par[i].m;
        //--------------------------------------------------------------------
        par[i].c_i = ci;
        par[i].c_j = cj;
        CX_local[ci][cj]+=(par[i].x*par[i].m);
        CY_local[ci][cj]+=(par[i].y*par[i].m);
        total_mass+=par[i].m;
    }

    size_t j,k;
    for (j = 0; j < ncside; j++){
        for (k = 0; k < ncside; k++){
            mass_C[j][k] = MC_local[j][k];
            if (MC_local[j][k] > EPS){
                Cx[j][k] = CX_local[j][k]/MC_local[j][k];
                Cy[j][k] = CY_local[j][k]/MC_local[j][k];
            }
            else{
                Cx[j][k] = (j+0.5)/ncside;
                Cy[j][k] = (k+0.5)/ncside;

            }
            //total_mass += MC_local[j][k];
            //printf("\nrow = %lu, col = %lu Cx = %.4f ,Cy = %.4f, CX2 = %.3f , CY2 = %.3f, \n",j,k,Cx[j][k],Cy[j][k],CX_local[j][k],CY_local[j][k]);
            CTX_local += CX_local[j][k];
            CTY_local += CY_local[j][k];
        }
    }
    *CTX = CTX_local/total_mass;
    *CTY = CTY_local/total_mass;
return total_mass;
}
//--------------------------------------------------------------------------
int main ( int argc , char* argv [ argc +1]) {
    long seed, ncside,ntstep;
    long long n_part;
    //-----------------------------------------------------------------------
    seed = atol(argv[1]);   //  seed for the random number generator
    ncside = atol(argv[2]); //  size of the grid (number of cells on the side)
    n_part = atol(argv[3]); //  number of particles
    ntstep = atol(argv[4]); //  number of time steps
    //-------------------------------------------------------------------------
    double Cx[ncside][ncside];
    double Cy[ncside][ncside];
    double mass_C[ncside][ncside];
    double total_mass=0.0,CTX = 0.0,CTY =0.0;
    particle_t par[n_part];
    //---------------------------------------------------------------------------
    total_mass=init_particles(seed, ncside, n_part, par,Cx, Cy,mass_C,&CTX,&CTY);
    //---------------------------------------------------------------------------
    size_t i,k,l;
    long long j;
    //printf("Fx = %f , Fy = %f\n",Fx,Fy);
    for(i = 0; i < ntstep; i++) { // loop over time
    double d2,d,dx,dy,Fx = 0,Fy=0,wl=0,wr=0,wu=0,wb=0,ax,ay;
    unsigned int ci,cj,cip1,cim1,cjm1,cjp1;
    double MC_local[ncside][ncside];
    double CX_local[ncside][ncside];
    double CY_local[ncside][ncside],CTX_local=0.0,CTY_local=0.0;
    for (l = 0; l < ncside; l++){
        for (k = 0; k < ncside; k++){
            MC_local[l][k] = 0;
            CX_local[l][k] = 0;
            CY_local[l][k] = 0;
        }
    }
    for(j = 0; j < n_part; j++) { // loop over particles
            ci=par[j].c_i;
            cj=par[j].c_j;
            //-----------------------------------------------------------------------
            //---------------- calculating index of naibour cell --------------------
            //-----------------------------------------------------------------------
             if (ci && ci < ncside-1 ){
                cip1 = ci+1;
                cim1 = ci-1;
            }
            else if(ci > ncside-2) {//
                cip1 = 0;
                cim1 = ci-1;
                wr = 1.0;
            }
            else    // ci == 0
            {
                cip1 = ci+1;
                cim1 = ncside-1;
                wl = 1.0;

            }
            //------------------------------------------------------------------------------
            if (cj && cj < ncside-1 ){
                cjp1 = cj+1;
                cjm1 = cj-1;
            }
            else if(cj > ncside-2) {//
                cjp1 = 0;
                cjm1 = cj-1;
                wu = 1.0;
            }
            else    // cj == 0
            {
                cjp1 = cj+1;
                cim1 = ncside-1;
                wb = 1.0;

            }
            //--------------------------------------------------------------------------------------------
            dx = (Cx[ci][cj]-par[j].x);
            dy = (Cy[ci][cj]-par[j].y);
            //printf("\nx = %f ,ci = %u, X = %f, m = %f , M = %f\n",par[j].x,ci,Cx[ci][cj],par[j].m,mass_C[ci][cj]);
            //printf("y = %f ,cj = %u, Y = %f\n",par[j].y,cj,Cy[ci][cj]);
            d2 = dx*dx+dy*dy;
            if (d2 >= EPSLON2){
                double m = par[j].m;
                d = sqrt(d2);
                //--------------------------------------------------------------
                Fx += ((G*m*mass_C[ci][cj])/d2)*(dx/d);  // 0
                Fy += ((G*m*mass_C[ci][cj]))/d2*(dy/d);  //
                //---------------------------------------------------------------
                dx = (Cx[cip1][cj]-par[j].x+wr);      //1
                dy = (Cy[cip1][cj]-par[j].y);
                d2 = dx*dx+dy*dy;
                d = sqrt(d2);
                Fx += ((G*m*mass_C[cip1][cj])/d2)*(dx/d);  // 1
                Fy += ((G*m*mass_C[cip1][cj]))/d2*(dy/d);  //
                //---------------------------------------------------------------
                dx = (Cx[cip1][cjm1]-par[j].x+wr);      //2
                dy = (Cy[cip1][cjm1]-par[j].y)-wb;
                d2 = dx*dx+dy*dy;
                d = sqrt(d2);
                Fx += ((G*m*mass_C[cip1][cjm1])/d2)*(dx/d);  // 2
                Fy += ((G*m*mass_C[cip1][cjm1]))/d2*(dy/d);  //
                //---------------------------------------------------------------
                dx = (Cx[cip1][cj]-par[j].x);      //3
                dy = (Cy[ci][cjm1]-par[j].y)-wb;
                d2 = dx*dx+dy*dy;
                d = sqrt(d2);
                Fx += ((G*m*mass_C[ci][cj])/d2)*(dx/d);  //
                Fy += ((G*m*mass_C[ci][cj]))/d2*(dy/d);  //
                //---------------------------------------------------------------
                dx = (Cx[cip1][cj]-par[j].x)-wl;      //4
                dy = (Cy[ci][cjm1]-par[j].y)-wb;
                d2 = dx*dx+dy*dy;
                d = sqrt(d2);
                Fx += ((G*m*mass_C[ci][cj])/d2)*(dx/d);  //
                Fy += ((G*m*mass_C[ci][cj]))/d2*(dy/d);  //
                //---------------------------------------------------------------
                dx = (Cx[cip1][cj]-par[j].x)-wl;      //5
                dy = (Cy[ci][cjm1]-par[j].y);
                d2 = dx*dx+dy*dy;
                d = sqrt(d2);
                Fx += ((G*m*mass_C[ci][cj])/d2)*(dx/d);  //
                Fy += ((G*m*mass_C[ci][cj]))/d2*(dy/d);  //
                //---------------------------------------------------------------
                dx = (Cx[cip1][cj]-par[j].x)-wl;      //6
                dy = (Cy[ci][cjm1]-par[j].y)+wu;
                d2 = dx*dx+dy*dy;
                d = sqrt(d2);
                Fx += ((G*m*mass_C[ci][cj])/d2)*(dx/d);  //
                Fy += ((G*m*mass_C[ci][cj]))/d2*(dy/d);  //
                //---------------------------------------------------------------
                dx = (Cx[cip1][cj]-par[j].x);      //7
                dy = (Cy[ci][cjm1]-par[j].y)+wu;
                d2 = dx*dx+dy*dy;
                d = sqrt(d2);
                Fx += ((G*m*mass_C[ci][cj])/d2)*(dx/d);  //
                Fy += ((G*m*mass_C[ci][cj]))/d2*(dy/d);  //
                //---------------------------------------------------------------
                dx = (Cx[cip1][cj]-par[j].x)+wr;      //8
                dy = (Cy[ci][cjm1]-par[j].y)+wu;
                d2 = dx*dx+dy*dy;
                d = sqrt(d2);
                Fx += ((G*m*mass_C[ci][cj])/d2)*(dx/d);  //
                Fy += ((G*m*mass_C[ci][cj]))/d2*(dy/d);  //
                //---------------------------------------------------------------
                //======== define new position and velocity ===================
                //============================================================
                ax = Fx/m;
                ay = Fy/m;
                par[j].vx+=ax;
                par[j].x+=par[j].vx+0.5*ax;
                if (par[j].x > 1){
                    par[j].x+= -1;
                }
                else if (par[j].x < 0){
                        par[j].x+=1;
                }
                //----------------------------------------------------------------------
                par[j].vy+=ay;
                par[j].y+=par[j].vy+0.5*ay;
                if (par[j].y > 1){
                    par[j].y+= -1;
                }
                else if (par[j].y < 0){
                        par[j].y+=1;
                }
            }
            //----------------------------------------------------------------------
            // calculate new cell and center of the cell
            ci = par[j].x*ncside;
            cj = par[j].y*ncside;
            //printf("\nparticle = %lu cell = %u,%u\tx = %0.5f, y = %0.5f\n",j,ci,cj,par[j].x,par[j].y);
            //--------------------------------------------------------------------
            MC_local[ci][cj]+=par[j].m;
            //--------------------------------------------------------------------
            par[j].c_i = ci;
            par[j].c_j = cj;
            CX_local[ci][cj]+=(par[j].x*par[j].m);
            CY_local[ci][cj]+=(par[j].y*par[j].m);
            //printf("CX = %0.5f CY = %0.5f\n", CX_local[ci][cj],CY_local[ci][cj]);
    }
    double totalmass2=0;
    for (l = 0; l < ncside; l++){
        for (k = 0; k < ncside; k++){
            mass_C[l][k] = MC_local[l][k];
            if (MC_local[l][k] > EPS){
                Cx[l][k] = CX_local[l][k]/MC_local[l][k];
                Cy[l][k] = CY_local[l][k]/MC_local[l][k];
            }
            else{
                Cx[l][k] = (l+0.5)/ncside;
                Cy[l][k] = (k+0.5)/ncside;

            }
            //printf("\nrow = %lu, col = %lu Cx = %.4f ,Cy = %.4f, CX2 = %.3f , CY2 = %.3f, M= %.3f\n",l,k,Cx[l][k],Cy[l][k],CX_local[l][k],CY_local[l][k],mass_C[l][k]);
            CTX_local += CX_local[l][k];
            CTY_local += CY_local[l][k];
        }
    }
    CTX = CTX_local/total_mass;
    CTY = CTY_local/total_mass;
    }
    printf("%.2f %.2f \n",par[0].x,par[0].y);
    printf("%.2f %.2f \n",CTX,CTY);
    return 0;
}  /* main */
