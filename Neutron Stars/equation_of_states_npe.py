"""
Description: Code to plot the different equations of state of an Ideal Fermi Gas of npe
"""
import numpy as np
import matplotlib.pyplot as plt

#constants
m_n, m_p, m_e = 1.67e-24, 1.6726219e-24, 9.10938e-28 # 939.5, 938.2 , 0.511 
l_n, l_p, l_e = 2.1e-14, 1.321e-13, 2.426e-10
c= 2.998e+10
hbar = 1.52-27

#Parameters
def N_cgs (m, l):
    return m*c**2/(np.pi**2*l**3) 

def N_mks (m, l):
    return N_cgs(m,l)*(1e-13)            

def N(m, l):
    return N_mks(m,l)

#Energy density
def epsilon_i(z, m, l):
    return N(m,l)/8*((2*z**3 + z)*np.sqrt(1+z*z) - np.arcsinh(z))

#Pressure
def Pressure_i(z, m, l):
    return N(m,l)/8*((2*z**3/3 - z)*np.sqrt(1+z*z) + np.arcsinh(z))

#Fermi momenta
def z_e (z):
    me = 0.511
    mp = 938.2
    return mp/me*z

def z_p (z):
    me = 0.511
    mp = 938.2
    mn = 939.5
    return np.sqrt((mn**2.*(1.+z**2.) - (me**2. + mp**2))**2 - 4*mp*me )/(2*mn*mp*np.sqrt(1.+z**2.))

def epsilon_NS(z):
    return epsilon_i(z, m_n, l_n) + epsilon_i(z_p(z), m_p, l_p) + epsilon_i(z_e(z_p(z)), m_e, l_e)

def Pressure_NS(z):
    return Pressure_i(z, m_n, l_n) + Pressure_i(z_p(z), m_p, l_p) + Pressure_i(z_e(z_p(z)), m_e, l_e)

#Plotting arrays
POW = 5
z_NS = np.zeros(9*POW)
for i in np.arange(POW):
    for j in np.arange(9):
        z_NS[i*9+j] = (j+1)*pow(10.,i-3.5)  


E_n = epsilon_i(z_NS, m_n, l_n)
E_npe = epsilon_NS(z_NS)
P_n = Pressure_i(z_NS, m_n, l_n)
P_npe = Pressure_NS(z_NS)


#Plotting functions
def Plot_energy(z, E):
    plt.title('Energy density $\epsilon$')
    plt.xlabel('z')
    plt.ylabel('$\epsilon \ [\epsilon_0]$')
    plt.loglog(z, E, color='c')
    plt.savefig('Plots/energy_density_npe.pdf')

def Plot_Presure(z, P):
    plt.title('Pressure $P$')
    plt.xlabel('z')
    plt.ylabel('$P \ [P_0]$')
    plt.loglog(z, P, 'C')
    plt.legend()
    plt.savefig('Plots/pressure_npe.pdf')

def Plot_EOS(P, E):
    plt.title('Equation of state (cgs)')
    plt.xlabel('$P \ [N/m^2]$')
    plt.ylabel('$\epsilon (P) \ [J/m^3]$')
    plt.loglog(P, E, color='c', label='npe ideal gas')
    plt.loglog(P_n, E_n, color='gray', linestyle='dashed', label='n ideal gas')
    plt.legend()
    plt.savefig('Plots/equation_of_state_npe.pdf')
    
    
#Output
Plot_EOS(P_npe, E_npe)







"""
//////////////////////////////////////////////////////////////////////////
Created on Fri Jan  1 13:32:11 2021                                    //
                                                                      //
@author: Nicolás Gómez                                               //
//////////////////////////////////////////////////////////////////////
"""

