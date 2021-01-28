"""
Description: Code to integrate TOV equations
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint as ode

#Relativistic equation of state
def epsilon(P):
    return 3*P    

#Parameters to TOV eqs
# in: cgs, mks, MeV, natural
G  = 6.6743e-11 
Ms =  1 #        1.989e+30 #in g
c  = 2.992e+8 #in cm/s
Pc =  1      #  3.3e+34   #in cgs   
alpha = 1.4829 # = G*Ms/c**2*10**(-5) in km
beta = 2.3289e-4

#def beta(Pc):
#    return 4*np.pi*Pc/(Ms*c*c)*

#Classical case, just to check
def Psi(y, r):
    P, M = y
    dydt = [-alpha*epsilon(P)*M, beta*epsilon(P)*r**2]
    return dydt


#Plotting functions
def Plot_mass(r, sol):
    plt.plot(r, sol[:, 1], 'g', label='M(r)')
    plt.legend(loc='best')
    plt.xlabel('r')
    plt.grid()

def Plot_pressure(r, sol):    
    plt.plot(r, sol[:, 0], 'r', label='P(r)')
    plt.legend(loc='best')
    plt.xlabel('r')
    plt.grid()
    
def Plot_density(r, sol):
    dener = np.zeros(101)
    for i in np.arange(101):
        dener[i] = epsilon(sol[i,0])

    plt.plot(r, dener, 'c', label='$\epsilon(r)$')
    plt.legend(loc='best')
    plt.xlabel('r')
    plt.grid()


def main():
    #Initial conditions
    y0 = [Pc, 0]

    #Plotting arrays
    radius = np.linspace(0, 15, 101)
    Solution = ode(Psi, y0, radius )

    #Plot
    Plot_density(radius, Solution)

#Output
main()







"""
//////////////////////////////////////////////////////////////////////////
Created on Fri Jan  1 13:32:11 2021                                    //
                                                                      //
@author: Nicolás Gómez                                               //
//////////////////////////////////////////////////////////////////////
"""

