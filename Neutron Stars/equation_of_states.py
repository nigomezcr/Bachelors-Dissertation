"""
Description: Code to plot the different equations of state of an Ideal Fermi Gas of n
"""
import numpy as np
import matplotlib.pyplot as plt

#constants
m = 1.67e-24
c= 2.998e10
hbar = 1.52-27
l_c = 2.1e-14
N_cgs = m*c**2/(np.pi**2*l_c**3) # ~1.635x10^49 erg/cm^3
N_mks = N_cgs*(1e-13)            # ~1.635x10^36 J/m^3
N = N_mks


#Energy density
def epsilon_NR(z):
    return 1/3*pow(z,3)

def epsilon_R(z):
    return 1/4*(pow(z,4) + pow(z,2))

def epsilon_NS(z):
    return 1/8*((2*z**3 + z)*np.sqrt(1+z*z) - np.arcsinh(z))


#Pressure
def Pressure_NR(z):
    return 1/15*pow(z,5)

def Pressure_R(z):
    return 1/12*pow(z,4)

def Pressure_NS(z):
    return 1/8*((2*z**3/3 - z)*np.sqrt(1+z*z) + np.arcsinh(z))


 #Equations of state 
def Equation_Of_State_NR(P):
    K = 5*pow(m/(15*np.pi*np.pi*l_c**3),0.4)
    g = 0.6
    return K*P**g    

def Equation_Of_State_R(P):
    return 3*P


#Plotting arrays
POW = 5
z_NS = np.zeros(9*POW)
for i in np.arange(POW):
    for j in np.arange(9):
        z_NS[i*9+j] = (j+1)*pow(10.,i-3.)  

POW = 2
z_NR = np.zeros(9*POW)
for i in np.arange(POW):
    for j in np.arange(9):
        z_NR[i*9+j] = (j+1)*pow(10.,i-3)  

POW = 2
z_R = np.zeros(9*POW)
for i in np.arange(POW):
    for j in np.arange(9):   
        z_R[i*9+j] = (j+1)*pow(10.,1)  


E_NR = epsilon_NR(z_NR)*N
E_R = epsilon_R(z_R)*N
E_NS = epsilon_NS(z_NS)*N

P_NR = Pressure_NR(z_NR)*N
P_R = Pressure_R(z_R)*N
P_NS = Pressure_NS(z_NS)*N

EOS_NR = Equation_Of_State_NR(P_NR)
EOS_R = Equation_Of_State_R(P_R)


##Plots
def Plot_pressure():
    plt.title('Pressure $P$')
    plt.xlabel('z')
    plt.ylabel('$P(z) \ [N/m^2]$')
    plt.loglog(z_NS, P_NS, 'k', label='Whole range')
    plt.loglog(z_R, P_R, color='m', linestyle='dashed', label='Relativistic')
    plt.loglog(z_NR, P_NR, color='r', linestyle='dashed', label='Non Relativistic')
    plt.legend()
    plt.savefig('Plots/pressure_n.pdf')

def Plot_energy():
    plt.title('Energy density $\epsilon$')
    plt.xlabel('z')
    plt.ylabel('$\epsilon(z) \ [J/m^3]$')
    plt.loglog(z_NS, E_NS, 'k', label='Whole range')
    plt.loglog(z_R, E_R, color='c', linestyle='dashed', label='Relativistic')
    plt.loglog(z_NR, E_NR, color='b', linestyle='dashed', label='Non Relativistic')
    plt.legend()
    plt.savefig('Plots/energy_density_n.pdf')

def Plot_eos():
    plt.title('Equation of state')
    plt.xlabel('$P \ [N/m^2]$')
    plt.ylabel('$\epsilon (P) \ [J/m^3]$')
    plt.loglog(P_NS, E_NS, color='black', label='Whole Range')
    plt.loglog(P_NR, E_NR, color='green', linestyle='dashed', label='Non Relativistic')
    plt.loglog(P_R, E_R, color='lime', linestyle='dashed', label='Relativistic')
    plt.legend()
    plt.savefig('Plots/equation_of_state_n.pdf')

#Output
Plot_eos()



"""
//////////////////////////////////////////////////////////////////////////
Created on Fri Jan  1 13:32:11 2021                                    //
                                                                      //
@author: Nicolás Gómez                                               //
//////////////////////////////////////////////////////////////////////
"""

