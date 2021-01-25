"""
Description:
"""
import numpy as np
import matplotlib.pyplot as plt

#constants
m_n, m_p, m_e = 1.67e-24, 1.6726219e-24, 9.10938e-28 # 939.5, 938.2 , 0.511 
c= 2.998e10
hbar = 1.52-27
l_n = hbar/(m_n*c) 
N_cgs = m_n*c**2/(np.pi**2*l_n**3) # ~1.635x10^43 erg/cm^3
N_mks = N_cgs*(1e-13)            # ~1.635x10^30 J/m^3
N = N_cgs

print(N)
#Energy density
def epsilon_i(z):
    return N/8*((2*z**3 + z)*np.sqrt(1+z*z) - np.arcsinh(z))

#Pressure
def Pressure_i(z):
    return N/8*((2*z**3/3 - z)*np.sqrt(1+z*z) + np.arcsinh(z))

#Fermi momenta
def z_e (z):
    return 1












"""
//////////////////////////////////////////////////////////////////////////
Created on Fri Jan  1 13:32:11 2021                                    //
                                                                      //
@author: Nicolás Gómez                                               //
//////////////////////////////////////////////////////////////////////
"""

