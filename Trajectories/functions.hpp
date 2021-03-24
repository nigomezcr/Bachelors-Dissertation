#ifndef __FUNCTIONS_HPP_
#define __FUNCTIONS_HPP_

#include <iostream>
#include "train.hpp"

const double g = -10;
const double R = 100;
const int NSTEPS = 500;
const double DT = 0.1;

// function declarations
void initial_conditions(Particle & body);
void compute_force(Particle & body);
void time_integration(Particle & body, const double & dt);
void start_integration(Particle & body, const double & dt);
void print(Particle & body, double time);

#endif // __FUNCTIONS_HPP_
