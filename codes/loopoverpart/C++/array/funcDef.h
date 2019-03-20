#include <iostream>
#include <stdlib.h>
#include <cmath>
#define RND0_1 ((double)random() / ((long long)1 << 31))
// Global consts
//#define RND0_1 ((double)random() / ((long long)1 << 31))
#define G 6.67408e-11
#define EPSLON2 0.0001
//#define EPSLON 0.01

using namespace std;

// Particle class
class particle_t
{
public:
  double x;         // x possition of the particle [0,1]
  double y;         // y possition of the particle [0,1]
  double vx;        // velocity of the particle in x direction
  double vy;        // velocity of the particle in y direction
  double m;         // mass of the particle
  size_t c_i; // particle's cell index 1 (row index)
  size_t c_j; // particle's cell index 2 (column index)
};

// Cell class
class cell_t
{
public:
  double x; // x possition of the cell CoM
  double y; // y possition of the cell CoM
  double m; // total mass of cell
};


/*------------------------------------------------------------------------------
 * Function:    init_particles
 * Purpose:     initialize the problem by randomly produce particles,
 *                and calculate the cell initial cell structure
 *
 * In arg:      seed for the random number generator
 *              ncside: size of the grid (number of cells3 on the side)
 *              number of particles
 *              pointer to struct particle
 *              pointer to struct cell
 */
void init_particles(long seed, size_t ncside, size_t n_part,particle_t *par,cell_t *cell)
{

  // Declarations
  size_t i;
  srandom(seed);

  // Loop over particles
  
  for (i = 0; i < n_part; i++)
  {
    // Initiliaze random particle variables
    par[i].x = RND0_1;
    par[i].y = RND0_1;
    par[i].vx = RND0_1 / ncside / 10.0;
    par[i].vy = RND0_1 / ncside / 10.0;
    par[i].m = RND0_1 * ncside / (G * 1e6 * n_part);
    //---------------------------------------------------------------
    par[i].c_i = par[i].x * ncside;
    par[i].c_j = par[i].y * ncside;
   //----------------------------------------------------------------
   //================================================================
   cell[par[i].c_i*ncside+par[i].c_j].m += par[i].m;
   cell[par[i].c_i*ncside+par[i].c_j].x += par[i].m*par[i].x;
   cell[par[i].c_i*ncside+par[i].c_j].y += par[i].m*par[i].y;
  }
  //================================================================
  // Loop trough cells to calculate CoM positions of each cell
  //#pragma omp parallel for 
  //#pragma omp simd
  for (size_t j = 0; j < ncside*ncside; j++)
  {
      if (cell[j].m) // Only consider cells with mass greater then eps
      {
        // Update cell center of mass positions using the total mass of the cell
        cell[j].x /= cell[j].m;
        cell[j].y /= cell[j].m;
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
 *              cells 
 */
inline void calculate_forces(size_t ci, size_t cj,size_t ncside, double xp, double yp, double m, double &Fx, double &Fy, const cell_t *cell)
{
  // Calculate distance from cell to point positions
  double dx = (cell[ci*ncside+cj].x - xp);
  double dy = (cell[ci*ncside+cj].y - yp);

  // Calcuate square ditance
  double d2 = dx * dx + dy * dy;

  // Check if the distance is bigger then epslon
  if (d2 >= EPSLON2) {
    // Indexes associated with neighboring cells
    size_t cip1, cim1, cjp1, cjm1;          // cip1: ci+1; cim1: ci-1; cjp1,cjm1 = (cj+1,cj-1)
    double wl = 0.0, wr = 0.0, wu = 0.0, wb = 0.0; // wl,wr,wu,wb went to the (left,right,up,bottom) of domain

    // Check location of neighboring particles (left and right)
    if (ci > 0 && ci < ncside - 1) // Particle is not in the right or left edge cells
    {
      cip1 = ci + 1;
      cim1 = ci - 1;
    }
    else if (ci == ncside - 1) // Particle is in the right edge
    {
      cip1 = 0;
      cim1 = ci - 1;
      wr = 1.0;
    }
    else // Particle is in the left edge, ci == 0
    {
      cip1 = ci + 1;
      cim1 = ncside - 1;
      wl = 1.0;
    }

    // Check location of neighboring particles (botton and top)
    if (cj > 0 && cj < ncside - 1) // Particle is not in the top or bottom edge cells
    {
      cjp1 = cj + 1;
      cjm1 = cj - 1;
    }
    else if (cj == ncside - 1) // Particle is in the top edge
    {
      cjp1 = 0;
      cjm1 = cj - 1;
      wu = 1.0;
    }
    else // Particle is in the bottom edge, cj == 0
    {
      cjp1 = cj + 1;
      cjm1 = ncside - 1;
      wb = 1.0;
    }

    double W = G*m;
    
    // Calculating forces
    // Intercation with own cell : 0
    double fd = sqrt(d2);
    Fx += ((W * cell[ci*ncside+cj].m)) / (d2 *fd)* dx ;
    Fy += ((W * cell[ci*ncside+cj].m)) / (d2 *fd) *dy;

    // Intercation with right cell (I+1,J) : 1
    dx = (cell[cip1*ncside+cj].x - xp + wr);
    dy = (cell[cip1*ncside+cj].y - yp);
    fd = pow(dx * dx + dy * dy,-1.5);
    Fx += W * cell[cip1*ncside+cj].m * fd *dx;
    Fy += W * cell[cip1*ncside+cj].m * fd *dy;

    // Interaction with botton right cell (I+1,J-1) : 2
    dx = cell[cip1*ncside+cjm1].x - xp + wr;
    dy = cell[cip1*ncside+cjm1].y - yp - wb;
    fd = pow(dx * dx + dy * dy,-1.5);
    Fx += W * cell[cip1*ncside+cjm1].m * fd *dx;
    Fy += W * cell[cip1*ncside+cjm1].m * fd *dy; //2

    // Interaction with botton cell (I,J-1) : 3
    dx = cell[ci*ncside+cjm1].x - xp;
    dy = cell[ci*ncside+cjm1].y - yp - wb;
    fd = pow(dx * dx + dy * dy,-1.5);
    Fx += W * cell[ci*ncside+cjm1].m * fd * dx;
    Fy += W * cell[ci*ncside+cjm1].m * fd * dy;

    // Interaction with left botton cell (I-1,J-1) : 4
    dx = (cell[cim1*ncside+cjm1].x - xp) - wl;
    dy = (cell[cim1*ncside+cjm1].y - yp) - wb;
    fd = pow(dx * dx + dy * dy,-1.5);
    Fx += W * cell[cim1*ncside+cjm1].m * fd * dx;
    Fy += W * cell[cim1*ncside+cjm1].m* fd * dy;

    // Interaction with left botton cell (I-1,J) : 5
    dx = (cell[cim1*ncside+cj].x - xp) - wl;
    dy = (cell[cim1*ncside+cj].y - yp);
    d2 = dx * dx + dy * dy;
    fd = pow(dx * dx + dy * dy,-1.5);
    Fx += W * cell[cim1*ncside+cj].m * fd * dx;
    Fy += W * cell[cim1*ncside+cj].m * fd * dy; //5

    // Interaction with left botton cell (I-1,J+1) : 6
    dx = (cell[cim1*ncside+cjp1].x - xp) - wl;
    dy = (cell[cim1*ncside+cjp1].y - yp) + wu;
    fd = pow(dx * dx + dy * dy,-1.5);
    Fx += W * cell[cim1*ncside+cjp1].m * fd * dx;
    Fy += W * cell[cim1*ncside+cjp1].m * fd * dy; //6

    // Interaction with left botton cell (I,J+1) : 7
    dx = (cell[ci*ncside+cjp1].x - xp); //7 (I,J+1)
    dy = (cell[ci*ncside+cjp1].y - yp) + wu;
    d2 = dx * dx + dy * dy;
    fd = pow(dx * dx + dy * dy,-1.5);
    Fx += W * cell[ci*ncside+cjp1].m * fd *dx;
    Fy += W * cell[ci*ncside+cjp1].m * fd *dy; //7

    // Interaction with left botton cell (I+1,J+1) : 8
    dx = (cell[cip1*ncside+cjp1].x - xp) + wr;
    dy = (cell[cip1*ncside+cjp1].y - yp) + wu;
    d2 = dx * dx + dy * dy;
    fd = pow(dx * dx + dy * dy,-1.5);
    Fx += W * cell[cip1*ncside+cjp1].m * fd * dx;
    Fy += W * cell[cip1*ncside+cjp1].m * fd * dy; //8
  }
}

/*------------------------------------------------------------------------------
 * Function:    update_velocties_and_positions
 * Purpose:     update the velocities and position of the ith particle 
 *
 * In arg:      i: index corresponding to particle
 *              Fx, Fy: forces in x and y direction
 *              par: vector containing all particles
 */
void update_velocities_and_positions(size_t i,size_t ncside, double Fx, double Fy, particle_t &par)
{
  double ax = Fx / par.m; // calculate the acceleration in x direction
  double ay = Fy / par.m; // calculate the acceleration in y direction

  // Update velocities and positions in "x" direction
  par.vx += ax;
  par.x += par.vx + 0.5 * ax;

  // Correct the "x" position
  // If it pass the right edge
  if (par.x > 1)
  {
    par.x -= (int)par.x;
  }
  // If it pass the left edge
  else if (par.x < 0)
  {
    par.x += (int)par.x + 1; // correct the position( if it pass the left edge)
  }

  // Update velocities and positions in "y" direction
  par.vy += ay;
  par.y += par.vy + 0.5 * ay;

  // Correct the "y" position
  // If it pass the top edge
  if (par.y > 1)
  {
    par.y -= (int)par.y; // is there better option?
  }
  // If it pass the botton edge
  else if (par.y < 0)
  {
    par.y += (int)par.y + 1;
  }
  par.c_i = par.x * ncside;
  par.c_j = par.y * ncside;
}