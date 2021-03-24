#include "functions.hpp"
#include <cmath>


void initial_conditions(Particle & body)
{
  body.Ry = R;
  body.Rz = 0.;
  body.Vx = 0.;
  body.Vy = 0.;
  body.Vz = 0.;
  body.rad = 0.1;
  body.mass = 8000;
}

void compute_force(Particle & body)
{
  
  // reset force and potential energy
  body.Fx = body.Fy = body.Fz = 0.0;

  // gravitational force
    //Constant density case
//  body.Fy += body.mass*g*body.Ry/R;

    //Effective density case: analytic form
  double b = 2.37411572646472;
  double c = 0.986950462357128;
  double d = 0.5881369116467612;
  double num = 3*g*b*(2*R*R - pow(1 - c*body.Ry/R, d+1)*(c*c*(1+d)*(2+d)*body.Ry*body.Ry + 2*c*(1+d)*body.Ry*R + 2*R*R) );
  double den = pow(c,3)*body.Ry*body.Ry*(6 + 11*d + 6*d*d + d*d*d);
  body.Fy += body.mass*num/den;



}

void start_integration(Particle & body, const double & dt)
{
  body.Vx -= body.Fx*dt/(2*body.mass);
  body.Vy -= body.Fy*dt/(2*body.mass);
  body.Vz -= body.Fz*dt/(2*body.mass);
}

  void time_integration(Particle & body, const double & dt)
{
  // leap-frog
  body.Vx += body.Fx*dt/(body.mass);
  body.Vy += body.Fy*dt/(body.mass);
  body.Vz += body.Fz*dt/(body.mass);
  body.Rx += body.Vx*dt;
  body.Ry += body.Vy*dt;
  body.Rz += body.Vz*dt;
}


void print(Particle & body, double time)
{
  std::cout << time << "  "
            << body.Rx << "  "
            << body.Ry << "  "
            << body.Rz << "  "
            << body.Vx << "  "
            << body.Vy << "  "
            << body.Vz << "  "
            << "\n";
}
