#ifndef __FUNCTIONS_HPP_
#define __FUNCTIONS_HPP_

#include <iostream>
#include "train.hpp"

const double g = -9.81;
const double R = 6.371e6;
const int NSTEPS = 2*2291;
const double DT = 1;

// function declarations
void initial_conditions(Particle & body);
void compute_force(Particle & body);
void time_integration(Particle & body, const double & dt);
void start_integration(Particle & body, const double & dt);
void print(Particle & body, double time);

#endif // __FUNCTIONS_HPP_
