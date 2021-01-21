"""
Description: Code to plot the different equations of state of an Ideal Fermi Gas
"""
import numpy as np
import matplotlib.pyplot as plt

#constants
m = 0.939 # MeV
c= 1.
hbar = 1
l_c = hbar/(m*c)
N = m*c**2/(np.pi**2*l_c**3)


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

def Equation_Of_State_NS(P):
    return 0


#Plotting arrays
z_NR = np.arange(0, 0.5, 0.001)
z_R = np.arange(50, 100, 0.1)
z_NS = np.arange(0, 100, 0.01)

E_NR = epsilon_NR(z_NR)
E_R = epsilon_R(z_R)
E_NS = epsilon_NS(z_NS)

P_NR = Pressure_NR(z_NR)
P_R = Pressure_R(z_R)
P_NS = Pressure_NS(z_NS)

EOS_NR = Equation_Of_State_NR(P_NR)
EOS_R = Equation_Of_State_R(P_R)


##Plots
"""
plt.title('Energy density $\epsilon$')
plt.xlabel('z')
plt.ylabel('$\epsilon \ (MeV)$')
plt.plot(z_NS, E_NS, 'k', label='Whole range')
plt.plot(z_R, E_R, color='c', linestyle='dashed', label='Relativistic')
plt.plot(z_NR, E_NR, color='b', linestyle='dashed', label='Non Relativistic')
plt.legend()

plt.title('Pressure $P$')
plt.xlabel('z')
plt.ylabel('$P \ (MeV)$')
plt.plot(z_NS, P_NS, 'k', label='Whole range')
plt.plot(z_R, P_R, color='m', linestyle='dashed', label='Relativistic')
plt.plot(z_NR, P_NR, color='r', linestyle='dashed', label='Non Relativistic')
plt.legend()
"""

plt.title('Equation of state')
plt.xlabel('$P \ (MeV)$')
plt.ylabel('$\epsilon \ (MeV)$')
plt.plot(P_NS, E_NS, color='k', label='Whole Range')
plt.plot(P_NR, E_NR, color='lime', linestyle='dashed', label='Non Relativistic')
plt.plot(P_R, E_R, color='green', linestyle='dashed', label='Non Relativistic')

plt.legend()








"""
//////////////////////////////////////////////////////////////////////////
Created on Fri Jan  1 13:32:11 2021                                    //
                                                                      //
@author: Nicolás Gómez                                               //
//////////////////////////////////////////////////////////////////////
"""

