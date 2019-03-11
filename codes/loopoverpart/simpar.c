# include <stdio.h>
# include <stdlib.h>
# include <stddef.h>
#include <math.h>
#include <time.h>
//------------------------------------------------------------------------------
#define RND0_1 ((double) random() / ((long long)1<<31))
#define G 6.67408e-11
#define EPSLON 0.01
#define EPSLON2 0.0001
//-----------------particle sturucture definnition -----------------------------
typedef struct particle_t particle_t;
struct particle_t {
  double x;     // x possition of the particle [0,1]
  double y;     // y possition of the particle [0,1]
  double vx;    // velocity of the particle in x direction
  double vy;    // velocity of the particle in y direction
  double m;     // mass of the particle
  unsigned long c_i; // particle's cell index 1 (row index)
  unsigned long c_j; // particle's cell index 2 (column index)
};
//-----------------cell sturucture definnition ---------------------------------
typedef struct cell_t cell_t;
struct cell_t {
  double m;     // mass of the cell
  double cx;    // y possition of center of the mass[0,1]
  double cy;    // x psoition of the center of the mass
  //double np;  // numper of partilces in each cell
  //unsigned long long part[]; array of number of particles inside the cell
};

/*------------------------------------------------------------------------------
 * Function:    init_particles
 * Purpose:     initialize the problem by randomly produce particles,
 *                and calculate the cell initial cell structure
 *
 * In arg:      seed for the random number generator
 *              ncside: size of the grid (number of cells on the side)
 *              number of particles
 *              pointer to struct particle
 *              pointer to struct cell
 */
void init_particles(long seed, long ncside, long long n_part, particle_t *par, cell_t cell[][ncside]){
    unsigned long long i;
    srandom(seed);
    double Mass_Cell_local[ncside][ncside];
    Mass_Cell_local[ncside][ncside] = 0.0;
    double CX_local[ncside][ncside];
    CX_local[ncside][ncside] = 0.0;
    double CY_local[ncside][ncside];
    CY_local[ncside][ncside] = 0.0;
    for(i = 0; i < n_part; i++) {
        par[i].x = RND0_1;
        par[i].y = RND0_1;
        par[i].vx = RND0_1 / ncside / 10.0;
        par[i].vy = RND0_1 / ncside / 10.0;
        //--------------------------------------------------------------------
        par[i].m = RND0_1 * ncside / (G * 1e6 * n_part);
        //======================== added code ===================================
        unsigned long ci = par[i].x*ncside;
        unsigned long cj = par[i].y*ncside;
        //--------------------------------------------------------------------
        Mass_Cell_local[ci][cj]+=par[i].m;
        //--------------------------------------------------------------------
        par[i].c_i = ci;
        par[i].c_j = cj;
        CX_local[ci][cj]+=(par[i].x*par[i].m);
        CY_local[ci][cj]+=(par[i].y*par[i].m);
        //total_mass+=par[i].m;
    }
    size_t j,k;
    for (j = 0; j < ncside; j++){
        for (k = 0; k < ncside; k++){
            cell[j][k].m = Mass_Cell_local[j][k];
            if (Mass_Cell_local[j][k] > __DBL_EPSILON__){
                cell[j][k].cx = CX_local[j][k]/Mass_Cell_local[j][k];
                cell[j][k].cy = CY_local[j][k]/Mass_Cell_local[j][k];
            }
            else{
                cell[j][k].cx = (j+0.5)/ncside;
                cell[j][k].cy = (k+0.5)/ncside;
            }
        }
    }
}

/*------------------------------------------------------------------------------
 * Function:    calcforce
 * Purpose:     Calculate the force on each particle from its own and neighbour
 *                cells
 *
 * In arg:      ci,cj: cell indexes associated to the particles
 *              ncside: size of the grid (number of cells on the side)
 *              xp,yp: coordinates of the particle
 *              mass of the partice
 *              Fx,Fy forces acting on the particle
 *              pointer to struct cell
 */
void calcforce(unsigned long ci,unsigned long cj,unsigned long ncside,double xp,double yp,double m,double *Fx,double *Fy,cell_t cell[][ncside]){
    double Fx_local = 0.0,Fy_local=0.0;
        /*
        Fx_local,Fy_local force in (x,y) direction
        */
    double dx = (cell[ci][cj].cx-xp);
    double dy = (cell[ci][cj].cy-yp);
    double d2 = dx*dx+dy*dy;  // calculating the square of distance
        /*
        d2 is squre of the distance of particle with each center of the mass of each cell
        dx,dy distance in (x,y) direction
        */
    if (d2 >= EPSLON2){ // if the distance is bigger than 0.01 otherwise F = 0
        //----------------------------------------------------------------------
        //------------- calculating indexes of neighbour cells -----------------
        //----------------------------------------------------------------------
        unsigned long cip1,cim1,cjp1,cjm1;
        /*
        ci,cj, : cell index of associated to current particle
        cip1: ci+1
        cim1: ci-1
        cjp1,cjm1 = (cj+1,cj-1)
        */
        double wl=0.0,wr=0.0,wu=0.0,wb=0.0;
        /*
        wl,wr,wu,wb went to the (left,right,up,bottom) of domain
        */
        if (ci && ci < ncside-1 ){
        // if the particle is not in the right or left edge cells
            cip1 = ci+1;
            cim1 = ci-1;
        }
        else if(ci > ncside-2) {// if the particle is in the right edge
            cip1 = 0;
            cim1 = ci-1;
            wr = 1.0;
        }
        else    // ci == 0 // if the particle is in the left edge
        {
            cip1 = ci+1;
            cim1 = ncside-1;
            wl = 1.0;
        }
        //----------------------------------------------------------------------
        if (cj && cj < ncside-1 ){
        // if the particle is not in the top or bottom edge cells
            cjp1 = cj+1;
            cjm1 = cj-1;
        }
        else if(cj > ncside-2) {//if the particle is in the top edge
            cjp1 = 0;
            cjm1 = cj-1;
            wu = 1.0;
        }
        else{    // cj == 0 //if the particle is in the bottom edge
            cjp1 = cj+1;
            cjm1 = ncside-1;
            wb = 1.0;
            }
        //----------------------------------------------------------------------
        double d = sqrt(d2);                           //
        Fx_local += ((G*m*cell[ci][cj].m)/d2)*(dx/d);  // 0 its own cell
        Fy_local += ((G*m*cell[ci][cj].m))/d2*(dy/d);  //
        //----------------------------------------------------------------------
        dx = (cell[cip1][cj].cx-xp+wr);      //1 (I+1,J)
        dy = (cell[cip1][cj].cy-yp);
        d2 = dx*dx+dy*dy;
        d = sqrt(d2);
        Fx_local += ((G*m*cell[cip1][cj].m)/d2)*(dx/d);
        Fy_local += ((G*m*cell[cip1][cj].m))/d2*(dy/d);  //1
        //----------------------------------------------------------------------
        dx = (cell[cip1][cjm1].cx-xp+wr);      //2 (I+1,J-1)
        dy = (cell[cip1][cjm1].cy-yp)-wb;
        d2 = dx*dx+dy*dy;
        d = sqrt(d2);
        Fx_local += ((G*m*cell[cip1][cjm1].m)/d2)*(dx/d);
        Fy_local += ((G*m*cell[cip1][cjm1].m))/d2*(dy/d);  //2
        //----------------------------------------------------------------------
        dx = (cell[ci][cjm1].cx-xp);      //3 (I,J-1)
        dy = (cell[ci][cjm1].cy-yp)-wb;
        d2 = dx*dx+dy*dy;
        d = sqrt(d2);
        Fx_local += ((G*m*cell[ci][cjm1].m)/d2)*(dx/d);
        Fy_local += ((G*m*cell[ci][cjm1].m))/d2*(dy/d);  //3
        //----------------------------------------------------------------------
        dx = (cell[cim1][cjm1].cx-xp)-wl;      //4 (I-1,J-1)
        dy = (cell[cim1][cjm1].cy-yp)-wb;
        d2 = dx*dx+dy*dy;
        d = sqrt(d2);
        Fx_local += ((G*m*cell[cim1][cjm1].m)/d2)*(dx/d);
        Fy_local += ((G*m*cell[cim1][cjm1].m))/d2*(dy/d);  //4
        //----------------------------------------------------------------------
        dx = (cell[cim1][cj].cx-xp)-wl;      //5 (I-1,J)
        dy = (cell[cim1][cj].cy-yp);
        d2 = dx*dx+dy*dy;
        d = sqrt(d2);
        Fx_local += ((G*m*cell[cim1][cj].m)/d2)*(dx/d);
        Fy_local += ((G*m*cell[cim1][cj].m))/d2*(dy/d);  //5
        //----------------------------------------------------------------------
        dx = (cell[cim1][cjp1].cx-xp)-wl;     //6 (I-1,J+1)
        dy = (cell[cim1][cjp1].cy-yp)+wu;
        d2 = dx*dx+dy*dy;
        d = sqrt(d2);
        Fx_local += ((G*m*cell[cim1][cjp1].m)/d2)*(dx/d);
        Fy_local += ((G*m*cell[cim1][cjp1].m))/d2*(dy/d);  //6
        //----------------------------------------------------------------------
        dx = (cell[ci][cjp1].cx-xp);         //7 (I,J+1)
        dy = (cell[ci][cjp1].cy-yp)+wu;
        d2 = dx*dx+dy*dy;
        d = sqrt(d2);
        Fx_local += ((G*m*cell[ci][cjp1].m)/d2)*(dx/d);
        Fy_local += ((G*m*cell[ci][cjp1].m))/d2*(dy/d);  //7
        //----------------------------------------------------------------------
        dx = (cell[cip1][cjp1].cx-xp)+wr;      //8 (I+1,J+1)
        dy = (cell[cip1][cjp1].cy-yp)+wu;
        d2 = dx*dx+dy*dy;
        d = sqrt(d2);
        Fx_local += ((G*m*cell[cip1][cjp1].m)/d2)*(dx/d);
        Fy_local += ((G*m*cell[cip1][cjp1].m))/d2*(dy/d);  //8
        //----------------------------------------------------------------------
    }
        *Fx = Fx_local;
        *Fy = Fy_local;
}
//------------------------------------------------------------------------------
int main ( int argc , char* argv [ argc +1]) {
    long seed, ncside,ntstep;
    long long n_part;
    //--------------------------------------------------------------------------
    seed = atol(argv[1]);   //  seed for the random number generator
    ncside = atol(argv[2]); //  size of the grid (number of cells on the side)
    n_part = atol(argv[3]); //  number of particles
    ntstep = atol(argv[4]); //  number of time steps
    //--------------------------------------------------------------------------
    clock_t start, end;  // to calculate the time
    start = clock();
    //--------------------------------------------------------------------------
    particle_t par[n_part];
    cell_t cell[ncside][ncside];
    //--------------------------------------------------------------------------
    init_particles(seed, ncside, n_part, par,cell); // initilization
    //--------------------------------------------------------------------------
    size_t t_step,j,k;
    unsigned long long i;
    for(t_step = 0; t_step < ntstep; t_step++) { // loop over time
        double Mass_Cell_local[ncside][ncside];
        double CX_local[ncside][ncside];
        double CY_local[ncside][ncside];
        for ( j = 0; j < ncside; j++){ // initilize the mass and center of the cells for each time step
            for (k = 0; k < ncside; k++){
                Mass_Cell_local[j][k] = 0;
                CX_local[j][k] = 0;
                CY_local[j][k] = 0;
            }
        }// End of initilize the mass of the cells for each time step
        for(i = 0; i < n_part; i++) { // loop over particles
            unsigned long ci=par[i].c_i; //cell index i
            unsigned long cj=par[i].c_j; //cell index j
            double xp = par[i].x; // x of the particle
            double yp = par[i].y; // y of the particel
            double m = par[i].m;  // mass of the particle
            //------------------------------------------------------------------
            double Fx = 0.0,Fy=0.0; // Fx,Fy force in (x,y) direction
            //------------------------------------------------------------------
            // calculate the Forces acting on each particle
            calcforce(ci,cj,ncside,xp,yp,m,&Fx,&Fy,cell);
            //------------------------------------------------------------------
            //=============================================================
            //======== define new position and velocity ===================
            //=============================================================
            double ax = Fx/m;  // calculate the acceleration in x direction
            double ay = Fy/m;  // calculate the acceleration in y direction
            //------------- x direction ----------------------------------------
            par[i].vx+=ax;     // update the velocity in x direction
            par[i].x+=par[i].vx+0.5*ax; // update the position in x direction
            if (par[i].x > 1){
                par[i].x-= (int)par[i].x; // correct the position( if it pass the right edge)
            }
            else if (par[i].x < 0){
                par[i].x+= (int)par[i].x+1;// correct the position( if it pass the left edge)
            }
            //------------------------------------------------------------------
            //--------------------- y direction -----------------------------
            par[i].vy+=ay;   // update the velocity in y direction
            par[i].y+=par[i].vy+0.5*ay; // update the position in y direction
            if (par[i].y > 1){
                par[i].y-= (int)par[i].x; // correct the position( if it pass the top edge)
            }
            else if (par[i].y < 0){
                par[i].y+= (int)par[i].y+1;//correct the position( if it pass the bottom edge)
            }
            //============================================================================
            //========  updating new cell association and center and mass of the cell ====
            //============================================================================
            //----------------------------------------------------------------------
            //
            ci = par[i].x*ncside;
            cj = par[i].y*ncside;
            //--------------------------------------------------------------------
            Mass_Cell_local[ci][cj]+=par[i].m;
            //--------------------------------------------------------------------
            par[i].c_i = ci;
            par[i].c_j = cj;
            CX_local[ci][cj]+=(par[i].x*par[i].m);
            CY_local[ci][cj]+=(par[i].y*par[i].m);
            //============================================================================
        } // end of loop over particles
        for (j = 0; j < ncside; j++){ //loop to update cell
            for (k = 0; k < ncside; k++){
                cell[j][k].m = Mass_Cell_local[j][k];
                if (Mass_Cell_local[j][k] > __DBL_EPSILON__){
                    cell[j][k].cx = CX_local[j][k]/Mass_Cell_local[j][k];
                    cell[j][k].cy = CY_local[j][k]/Mass_Cell_local[j][k];
                }
                else{
                    cell[j][k].cx = (j+0.5)/ncside;
                    cell[j][k].cy = (k+0.5)/ncside;
                }
            }
        } // end loop to update cell
    } // end of loop over time
    double total_mass=0.0,TotalCenter_x = 0.0,TotalCenter_y =0.0;
    for (j = 0; j < ncside; j++){ //loop to gain total center of mass
        for (k = 0; k < ncside; k++){
            TotalCenter_x += cell[j][k].cx*cell[j][k].m;
            TotalCenter_y += cell[j][k].cy*cell[j][k].m;
            total_mass +=cell[j][k].m;
        }
    } //end of loop to gain total center of mass
    TotalCenter_x/= total_mass;
    TotalCenter_y/= total_mass;
    printf("%.2f %.2f \n",par[0].x,par[0].y);
    printf("%.2f %.2f \n",TotalCenter_x,TotalCenter_y);
    //----------------------------------------------------------------------------------
    end = clock();
    double time_spent = ((double) (end - start)) / CLOCKS_PER_SEC;
    //---------------------------------------------------------------------------------
    printf("It took %.4f seconds\n",time_spent);
    return EXIT_SUCCESS;
}  /* main */
