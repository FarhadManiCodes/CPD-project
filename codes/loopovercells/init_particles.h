#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <ctime>

#define RND0_1 ((double) random() / ((long long)1<<31))
#define G 6.67408e-11
#define EPSLON 0.01  

using namespace std;

class particle_t {
      
  public:
    //! Constructor.
    particle_t() {};
    
    //! Destructor.
    ~particle_t(){};
    
    double x;                                               /* In */
    double y;                                               /* In */
    double vx;                                              /* In */
    double vy;                                              /* In */
    double m;                                               /* In */
    
    double fx;                                               /* In force */
    double fy;
    double ax;                                               /* In acceleration*/
    double ay;                                               /* In acceleration*/
    
  private:
    
};


class cell_t {

  public:
    
    //! Constructor.
    cell_t(){};
    
//     cell_t() : neighbourCells(9) {};
    
    //! Destructor.
    ~cell_t(){};
    
    vector<int> particlesInCell;
    vector<int> neighbourCells;
    
    void set_particle(int particleID){
      particlesInCell.push_back(particleID);
    };
        
//     void set_neighbourCell(int index, int cellID){
//       neighbourCells[index] = cellID;
// //     void set_neighbourCell(int cellID){
// //       neighbourCells.push_back(cellID);
//     };
    
    void add_neighbourCells(vector<int> &cells){
      neighbourCells = cells;
    };
    
    double mass;
    double tempX;
    double tempY;
    
};

class overall_t {
  public:
    
    //! Constructor.
    overall_t(){};
    
    //! Destructor.
    ~overall_t(){};
    
    double mass;
    double centerOfMass_X;
    double centerOfMass_Y;
    
};

void init_particles(long seed, long ncside, long long n_part, particle_t *par)
{
    long long i;
    srandom(seed);

    for(i = 0; i < n_part; i++)
    {
        par[i].x = RND0_1;
        par[i].y = RND0_1;
        par[i].vx = RND0_1 / ncside / 10.0;
        par[i].vy = RND0_1 / ncside / 10.0;
        par[i].m = RND0_1 * ncside / (G * 1e6 * n_part);
    }
};

void computeOverallCenterOfMassPositionAndMass(int n_part, particle_t *par, overall_t *overall) {

    double overallMass = 0.0;
    double tempX = 0.0;
    double tempY = 0.0;
    
    for(int i = 0; i < n_part; i++){
      
      overallMass += par[i].m;
      tempX += par[i].m * par[i].x;
      tempY += par[i].m * par[i].y;
      
    }
        
    overall->mass = overallMass;
    overall->centerOfMass_X = tempX/overallMass;
    overall->centerOfMass_Y = tempY/overallMass;
      
};


class simpar_t {
  public:
    
    //! Constructor.
    simpar_t(){};
    
    //! Destructor.
    ~simpar_t(){};
    
    void init_NeighbourCells(int ncside, cell_t *cell) {
      
      vector< vector<int> > blockMatrix(ncside+2);
      
      //Initialize matrix to zero
      for(int i=0 ; i<ncside+2 ; ++i)
        blockMatrix[i].resize(ncside+2, 0);
      
      
      //Set cell numbering for inside matrix 
      int cellIndex = 0;
      for(int i=1 ; i<ncside+1 ; ++i){
        for(int j=1 ; j<ncside+1 ; ++j){
            blockMatrix[i][j] = cellIndex++;
        }
      }
      
      //Set cell numbering for matrix Edges   
      for(int j=1 ; j<ncside+1 ; ++j){
        blockMatrix[0][j] = blockMatrix[ncside][j];   //Bottom == Top
        blockMatrix[ncside+1][j] = blockMatrix[1][j];  //Top == Bottom
        blockMatrix[j][0] = blockMatrix[j][ncside];    //Left == Rigth
        blockMatrix[j][ncside+1] = blockMatrix[j][1];   //Rigth == Left
      }
      
      //Set cell numbering for matrix Corners  
      blockMatrix[0][0] = blockMatrix[ncside][ncside];      //Bottom Left  == Up Right 
      blockMatrix[ncside+1][ncside+1] = blockMatrix[1][1];  //Up Right == Bottom Left
      blockMatrix[0][ncside+1] = blockMatrix[ncside][1];    //Bottom Right  == Up Left
      blockMatrix[ncside+1][0] = blockMatrix[1][ncside];    //Up Left == Bottom Right
      
      
      //Assign Neighbours to each CELL
      cellIndex = 0;
      for(int i=1 ; i<ncside+1 ; ++i){
        for(int j=1 ; j<ncside+1 ; ++j){
              
          vector<int> neighbourCells;      
          neighbourCells.push_back(blockMatrix[i-1][j-1]);
          neighbourCells.push_back(blockMatrix[i-1][j]);
          neighbourCells.push_back(blockMatrix[i-1][j+1]);
          neighbourCells.push_back(blockMatrix[i][j-1]);
          neighbourCells.push_back(blockMatrix[i][j]);
          neighbourCells.push_back(blockMatrix[i][j+1]);
          neighbourCells.push_back(blockMatrix[i+1][j-1]);
          neighbourCells.push_back(blockMatrix[i+1][j]);
          neighbourCells.push_back(blockMatrix[i+1][j+1]);
          
          cell[cellIndex].add_neighbourCells(neighbourCells);
          
          
          cellIndex +=1;
        }
      }
    };
    
    void computeParticlesAndMassCenterInEachCell(int ncside, int n_part, particle_t *par, cell_t *cell) {
      
      double x, y, mass;
      int posX, posY, position;
        
      for(int i = 0; i < n_part; i++){
        x = par[i].x;
        y = par[i].y;
        mass = par[i].m;
        
        posX = x*ncside;
        posY = y*ncside;
        position = posY*ncside + posX;
        cell[position].set_particle(i);
        cell[position].mass += mass;
        cell[position].tempX += mass*x;
        cell[position].tempY += mass*y;
      
      }
    };

    
    void computeForceAndPositions(int ncside, particle_t *par, cell_t *cell) {

      int numberofCells = ncside*ncside;
      double particleMass, particleX, particleY;
      double massOfCell, centerOfMassX, centerOfMassY;
      double dX, dY, f;
      double currentDistanceSquared, currentDistance;
      
      for(int i = 0; i < numberofCells; i++){

        for (std::vector<int>::const_iterator particleIt = cell[i].particlesInCell.begin(); particleIt != cell[i].particlesInCell.end(); ++particleIt){
          particleMass = par[*particleIt].m;
          particleX = par[*particleIt].x;
          particleY = par[*particleIt].y;
          
          par[*particleIt].fx = 0;
          par[*particleIt].fy = 0;
          
          
          for (std::vector<int>::const_iterator neighbourCellIt = cell[i].neighbourCells.begin(); neighbourCellIt != cell[i].neighbourCells.end(); ++neighbourCellIt){
            
            if (cell[*neighbourCellIt].particlesInCell.size() != 0){
            
              massOfCell = cell[*neighbourCellIt].mass;
              centerOfMassX = cell[*neighbourCellIt].tempX / massOfCell ;
              centerOfMassY = cell[*neighbourCellIt].tempY / massOfCell ;
              
              dX = centerOfMassX-particleX;
              dY = centerOfMassY-particleY;
              currentDistanceSquared = (dX)*(dX) + (dY)*(dY);
              currentDistance = sqrt(currentDistanceSquared);
            
              
              if (currentDistanceSquared < EPSLON)
                continue;
              else
                f = G*particleMass*massOfCell/currentDistanceSquared;
                
              
              par[*particleIt].fx += f * dX / currentDistance;
              par[*particleIt].fy += f * dY / currentDistance;
              
            }
            
          }
          
          par[*particleIt].ax = par[*particleIt].fx / par[*particleIt].m;
          par[*particleIt].ay = par[*particleIt].fy / par[*particleIt].m ;
          
          par[*particleIt].vx = par[*particleIt].vx + par[*particleIt].ax * 1;
          par[*particleIt].vy = par[*particleIt].vy + par[*particleIt].ay * 1;
          
          par[*particleIt].x = par[*particleIt].x + (par[*particleIt].vx*1) + (0.5*par[*particleIt].ax*1*1);
          par[*particleIt].y = par[*particleIt].y + (par[*particleIt].vy*1) + (0.5*par[*particleIt].ay*1*1);
          
          if (par[*particleIt].x > 1.0)
            par[*particleIt].x = par[*particleIt].x - 1.0;
          else if (par[*particleIt].x < 0)
            par[*particleIt].x = par[*particleIt].x + 1.0;
          
          if (par[*particleIt].y > 1.0)
            par[*particleIt].y = par[*particleIt].y - 1.0;
          else if (par[*particleIt].y < 0)
            par[*particleIt].y = par[*particleIt].y + 1.0;
          
          
        }
      }
    };
    
};


// extern void init_particles(long seed, long ncside, long long n_part, particle_t *par);
// extern void computeParticlesAndMassCenterInEachCell(int ncside, int n_part, particle_t *par, cell_t *cell);
// extern void init_NeighbourCells(int ncside, cell_t *cell);
// extern void computeCenterOfMassPositionAndMassPerEachCell(int ncside, particle_t *par, cell_t *cell);
// extern void computeOverallCenterOfMassPositionAndMass(int n_part, particle_t *par, overall_t *overall);
// extern void computeForceAndPositions(int ncside, particle_t *par, cell_t *cell);
// extern void computeAccelerationAndVelocityAndDisplacement(int ncside, particle_t *par, cell_t *cell, double dt);


// /*LBLAS OPERATORS */
// 
// extern "C" //This is important to get the function from the -lblas library you will use when compiling
// {
//     double daxpy_(int *n, double *a, double *A, int *incA, double *B, int *incB);
// //The daxpy fortran function shown above multiplies a first matrix 'A' by a constant 'a'
// //and adds the result to a second matrix 'B.' Both matrices are of size 'n.'
// };
// 
// void daxpy(int n, double a, double *A, int incA, double *B, int incB)
// {
//     daxpy_(&n, &a, A, &incA, B, &incB); //Once again, note the call notation. Important!
// }
  
  

