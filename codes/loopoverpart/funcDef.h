#include <iostream>
#include <stdlib.h>
//#include <stddef.h>
#include <cmath>
#include <vector>

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
  unsigned int c_i; // particle's cell index 1 (row index)
  unsigned int c_j; // particle's cell index 2 (column index)
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
 * Function:    locate_and_update_cell_info
 * Purpose:     locate cell belonging to ith particle and update its information 
 *
 * In arg:      i: index corresponding to particle
 *              ncside: number of cells in each side
 *              par: vector containing all particles
 *              cell: matrix containing
 */
void locate_and_update_cell_info(unsigned long i, unsigned int ncside, vector<particle_t> &par, vector<vector<cell_t>> &cell)
{
  // Determine cell where the particle is located
  unsigned int ci = par[i].x * ncside;
  unsigned int cj = par[i].y * ncside;

  // Add particles mass to cell mass
  cell[ci][cj].m += par[i].m;

  // Add (x,y)*m to cell position ** this will have to be devided by the total mass of cell later
  cell[ci][cj].x += (par[i].x * par[i].m);
  cell[ci][cj].y += (par[i].y * par[i].m);

  // Assign cell ids to particle info
  par[i].c_i = ci;
  par[i].c_j = cj;
}

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
void init_particles(long seed, unsigned int ncside, unsigned long n_part, vector<particle_t> &par, vector<vector<cell_t>> &cell)
{
  #define RND0_1 ((double)random() / ((long long)1 << 31))
  // Declarations
  unsigned long i;
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

    locate_and_update_cell_info(i, ncside, par, cell);
  }

  // Loop trough cells to calculate CoM positions of each cell
  for (unsigned int j = 0; j < ncside; j++)
  {
    for (unsigned int k = 0; k < ncside; k++)
    {
      // ****************** CHECK THIS ******************
      // What is DBL_EPSILON?
      if (cell[j][k].m > __DBL_EPSILON__) // Only consider cells with mass greater then eps
      {
        // Update cell center of mass positions using the total mass of the cell
        cell[j][k].x /= cell[j][k].m;
        cell[j][k].y /= cell[j][k].m;
      }
      // If we multiple by 0 afterwards why do we need to set a mass center?
      /* else // If the mass of the cell is too small put the center in the middle of the cell (This is not nescessary)
      {
        cell[j][k].x = (j + 0.5) / ncside;
        cell[j][k].y = (k + 0.5) / ncside;
      } */
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
void calculate_forces(unsigned int ci, unsigned int cj,unsigned int ncside, double xp, double yp, double m, double &Fx, double &Fy, const vector<vector<cell_t>> &cell)
{
  // Calculate distance from cell to point positions
  double dx = (cell[ci][cj].x - xp);
  double dy = (cell[ci][cj].y - yp);

  // Calcuate square ditance
  double d2 = dx * dx + dy * dy;

  // Check if the distance is bigger then epslon
  if (d2 >= EPSLON2) {
    // Indexes associated with neighboring cells
    unsigned int cip1, cim1, cjp1, cjm1;          // cip1: ci+1; cim1: ci-1; cjp1,cjm1 = (cj+1,cj-1)
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

    // Calculating forces
    // Intercation with own cell : 0
    double d = sqrt(d2);
    Fx += ((G * m * cell[ci][cj].m) / d2) * (dx / d);
    Fy += ((G * m * cell[ci][cj].m)) / d2 * (dy / d);

    // Intercation with right cell (I+1,J) : 1
    dx = (cell[cip1][cj].x - xp + wr);
    dy = (cell[cip1][cj].y - yp);
    d2 = dx * dx + dy * dy;
    d = sqrt(d2);
    Fx += ((G * m * cell[cip1][cj].m) / d2) * (dx / d);
    Fy += ((G * m * cell[cip1][cj].m)) / d2 * (dy / d);

    // Interaction with botton right cell (I+1,J-1) : 2
    dx = (cell[cip1][cjm1].x - xp + wr);
    dy = (cell[cip1][cjm1].y - yp) - wb;
    d2 = dx * dx + dy * dy;
    d = sqrt(d2);
    Fx += ((G * m * cell[cip1][cjm1].m) / d2) * (dx / d);
    Fy += ((G * m * cell[cip1][cjm1].m)) / d2 * (dy / d); //2

    // Interaction with botton cell (I,J-1) : 3
    dx = (cell[ci][cjm1].x - xp);
    dy = (cell[ci][cjm1].y - yp) - wb;
    d2 = dx * dx + dy * dy;
    d = sqrt(d2);
    Fx += ((G * m * cell[ci][cjm1].m) / d2) * (dx / d);
    Fy += ((G * m * cell[ci][cjm1].m)) / d2 * (dy / d);

    // Interaction with left botton cell (I-1,J-1) : 4
    dx = (cell[cim1][cjm1].x - xp) - wl;
    dy = (cell[cim1][cjm1].y - yp) - wb;
    d2 = dx * dx + dy * dy;
    d = sqrt(d2);
    Fx += ((G * m * cell[cim1][cjm1].m) / d2) * (dx / d);
    Fy += ((G * m * cell[cim1][cjm1].m)) / d2 * (dy / d);

    // Interaction with left botton cell (I-1,J) : 5
    dx = (cell[cim1][cj].x - xp) - wl;
    dy = (cell[cim1][cj].y - yp);
    d2 = dx * dx + dy * dy;
    d = sqrt(d2);
    Fx += ((G * m * cell[cim1][cj].m) / d2) * (dx / d);
    Fy += ((G * m * cell[cim1][cj].m)) / d2 * (dy / d); //5

    // Interaction with left botton cell (I-1,J+1) : 6
    dx = (cell[cim1][cjp1].x - xp) - wl;
    dy = (cell[cim1][cjp1].y - yp) + wu;
    d2 = dx * dx + dy * dy;
    d = sqrt(d2);
    Fx += ((G * m * cell[cim1][cjp1].m) / d2) * (dx / d);
    Fy += ((G * m * cell[cim1][cjp1].m)) / d2 * (dy / d); //6

    // Interaction with left botton cell (I,J+1) : 7
    dx = (cell[ci][cjp1].x - xp); //7 (I,J+1)
    dy = (cell[ci][cjp1].y - yp) + wu;
    d2 = dx * dx + dy * dy;
    d = sqrt(d2);
    Fx += ((G * m * cell[ci][cjp1].m) / d2) * (dx / d);
    Fy += ((G * m * cell[ci][cjp1].m)) / d2 * (dy / d); //7

    // Interaction with left botton cell (I+1,J+1) : 8
    dx = (cell[cip1][cjp1].x - xp) + wr;
    dy = (cell[cip1][cjp1].y - yp) + wu;
    d2 = dx * dx + dy * dy;
    d = sqrt(d2);
    Fx += ((G * m * cell[cip1][cjp1].m) / d2) * (dx / d);
    Fy += ((G * m * cell[cip1][cjp1].m)) / d2 * (dy / d);
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
void update_velocities_and_positions(unsigned long i, double Fx, double Fy, vector<particle_t> &par)
{
  double ax = Fx / par[i].m; // calculate the acceleration in x direction
  double ay = Fy / par[i].m; // calculate the acceleration in y direction

  // Update velocities and positions in "x" direction
  par[i].vx += Fx / par[i].m;
  par[i].x += par[i].vx + 0.5 * ax;

  // Correct the "x" position
  // If it pass the right edge
  if (par[i].x > 1)
  {
    par[i].x -= (int)par[i].x;
  }
  // If it pass the left edge
  else if (par[i].x < 0)
  {
    par[i].x += (int)par[i].x + 1; // correct the position( if it pass the left edge)
  }

  // Update velocities and positions in "y" direction
  par[i].vy += ay;
  par[i].y += par[i].vy + 0.5 * ay;

  // Correct the "y" position
  // If it pass the top edge
  if (par[i].y > 1)
  {
    par[i].y -= (int)par[i].y;
  }
  // If it pass the botton edge
  else if (par[i].y < 0)
  {
    par[i].y += (int)par[i].y + 1;
  }
}

/*------------------------------------------------------------------------------
 * Function:    update_global_quantities
 * Purpose:     update global mass and CoM  
 *
 * In arg:      i: index corresponding to particle
 *              TotalCenter_:  position of center of mass
 *              total: total mass of the problem
 */
void update_global_quantities(unsigned int ncside, double &TotalCenter_x, double &TotalCenter_y, double &total_mass, const vector<particle_t> &par, const vector<vector<cell_t>> &cell)
{
  // Loop trough cells
  for (unsigned int j = 0; j < ncside; j++)
  {
    for (unsigned int k = 0; k < ncside; k++)
    {
      // Calculate info of each cell
      TotalCenter_x += cell[j][k].x * cell[j][k].m;
      TotalCenter_y += cell[j][k].y * cell[j][k].m;
      total_mass += cell[j][k].m;
    }
  }

  // Update positions
  TotalCenter_x /= total_mass;
  TotalCenter_y /= total_mass;
}
/* 
void zero_cells(long ncside, vector<vector<cell_t>> &cell_aux)
{
  for (size_t j = 0; j < ncside; j++)
  {
    for (size_t k = 0; k < ncside; k++)
    {
      // Calculate info of each cell
      cell_aux[j][k].x = 0;
      cell_aux[j][k].y = 0;
      cell_aux[j][k].m = 0;
    }
  }
} */