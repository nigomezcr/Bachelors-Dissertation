#include "functions.hpp"
#include <string>
#include <fstream>

int main(void)
{
  Particle train;

  // start system
  initial_conditions(train);
  compute_force(train);
  //print(train, 0.0);
  start_integration(train, DT);

  // evolve
  for(int istep = 0; istep < NSTEPS; ++istep) {
    time_integration(train, DT);
    compute_force(train);
    print(train, istep*DT);
  }
  
  return 0;
}
