#include <iostream>
#include <iomanip>
#include <sstream>
// #include <stdio.h>
// #include <stdlib.h>
#include <stdexcept>
#include "init_particles.h"


using namespace std;


int main(int argc, char *argv[]) {

//   for(int i = 0; i < argc; i++){
// //     std::cout << "Argument "<< i << " = " << argv[i] << std::endl;
//   }
  
  
  if (argc != 5) {
    throw std::invalid_argument( "Incomplete arguments" );
  }
  
  
  unsigned int randomNumber, numberOfSideCells, numberOfParticles, numberOfTimeSteps;
  for(int i = 1; i < argc; i++){
    stringstream ss(argv[i]);
    if (i == 1)
      ss >> randomNumber;
    else if (i == 2)
      ss >> numberOfSideCells;
    else if (i == 3)
      ss >> numberOfParticles;
    else if (i == 4)
      ss >> numberOfTimeSteps;
  }
  
  clock_t time_req;
  time_req = clock();
  
  //Initialize Particle Members
  vector<particle_t> myInstance(numberOfParticles);
  particle_t *particle_ptr = &myInstance[0];
  init_particles(randomNumber, numberOfSideCells, numberOfParticles, particle_ptr);
  

  //Initialize Overall Member
  overall_t overall;
  overall_t *overall_ptr = &overall;
  
  
  //Initialize Cells Objects 
  vector<cell_t> myCells(numberOfSideCells*numberOfSideCells);
  cell_t *cell_ptr = &myCells[0];
  
  
  //Initialize simpar object
  simpar_t simpar;
  
  //Create cell Neighbourhood for each cell
  simpar.init_NeighbourCells(numberOfSideCells, cell_ptr);
  
  bool lastStep = false;
  for(int step=0; step<numberOfTimeSteps; step++){
    
    if (step == numberOfTimeSteps-1)
      lastStep = true;
    
         
    //Assign particles to each cell and compute Cell Mass
    simpar.computeParticlesAndMassCenterInEachCell(numberOfSideCells, numberOfParticles, particle_ptr, cell_ptr, overall_ptr, lastStep);
  
    //Compute Force, Velocities and New Positions
    simpar.computeForceAndPositions(numberOfSideCells, particle_ptr, cell_ptr);
    
  
    if(lastStep){
      
      cout << fixed;
      cout << setprecision(5);
      cout << particle_ptr[0].x  << " " << particle_ptr[0].y << std::endl; 
    
      cout << fixed;
      cout << setprecision(5);
      cout << overall_ptr[0].tempX / overall_ptr[0].mass  << " " << overall_ptr[0].tempY / overall_ptr[0].mass << std::endl;
      
    }
    
    //Clean Data
    for(int i = 0; i < numberOfSideCells*numberOfSideCells; i++){
      cell_ptr[i].particlesInCell.clear();
      cell_ptr[i].mass = 0;
      cell_ptr[i].tempX = 0;
      cell_ptr[i].tempY = 0;
    }
    

  }
  cout << fixed;
  cout << setprecision(3);
  time_req = clock()- time_req;
  cout << "It took " << (float)time_req/CLOCKS_PER_SEC << " seconds" << endl;
  
  
  
  
  return 0;
}
