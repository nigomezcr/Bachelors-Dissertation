"""
Description: Code to integrate TOV equations
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint as ode
from scipy.interpolate import interp1d as interpolate
from scipy.optimize import root_scalar as findroot


#Parameters to TOV eqs
# In: cgs, mks, MeV, comp
    #Physical constants
cgs, mks, MeV, comp = 0, 1, 2, 3
G  = (6.6743e-8, 6.6743e-11, 6.7e-45, 1)
c  = (2.9979e10, 2.99792e8 , 1      , 1) 
hb = (1.0546-27, 1.0546e-35, 1      , 1)

    #Constants of reference of the system
Ms = (1.9891e33, 1.9891e30, 1.1e60  , 1) 
Pc = (3.3e+34  , 3.3e+33  , 2.68e6  , 1) 

mn = (1.675e-24, 1.675e-27, 939.5   , 1)
mp = (1.673e-24, 1.673e-27, 938.2   , 0.9986)
me = (9.109e-28, 9.109e-31, 0.511   , 5.44e-4)

ln = (2.1e-14  , 2.1e-16  , 1.064e-3, 1)
lp = (1.321e-13, 1.321e-15, 1.066e-3, 1.001)
le = (2.4252-10, 2.425e-12, 1.956   , 1838.23)

     #Specify system being used
SYS = cgs

alpha = 1.47716 #in km
def beta(Pc):
    km = (10**(-5), 10**(-3))
    return 4*np.pi*Pc/(Ms[SYS]*c[SYS]**2)*(1/km[SYS])**3



#Equations of State

    #Non relativistic EOS for n
def Equation_Of_State_NR(P):
    K = 5*pow(mn[SYS]/(15*np.pi*np.pi*ln[SYS]**3),0.4)
    g = 0.6
    return K*P**g    

    #Relativistic equation of state for n
def Equation_Of_State_R(P):
    return 3*P    

def Equation_Of_State_NS(P):
    POW = 3
    z_i = np.zeros(9*POW)
    
    for i in np.arange(POW):
        for j in np.arange(9):
            z_i[i*9+j] = (j+1)*pow(10.,i-3)  
        
    z = np.insert(z_i, 0, 0)
    
    N = 1#*mn[SYS]*c[SYS]**2/(np.pi**2*ln[SYS]**3)
    P_NS = N/8*((2*z**3/3 - z)*np.sqrt(1+z*z) + np.arcsinh(z))
    E_NS = N/8*((2*z**3 + z)*np.sqrt(1+z*z) - np.arcsinh(z))
        
    
    f = interpolate(P_NS, E_NS)
    
    return f(P)
   
      
#Hydrostatic Equations

    #Classical case, just to check
def Newton(y, r, Equation_Of_State):
    P, M = y
    dydt = [-alpha*Equation_Of_State(P)*M, beta(Pc[SYS])*Equation_Of_State(P)*r**2]
    return dydt

    #TOV equations
def Tolman_Oppenheimer_Volkoff(y, r, Equation_Of_State):
    P, M = y
    dydt = [-alpha*(M + beta(Pc[SYS])*P*r**3)*(Equation_Of_State(P) + P)*(r**2 - 2*alpha*r*M)**(-1) , beta(Pc[SYS])*Equation_Of_State(P)*r**2]
    return dydt



#Plotting functions
def Plot_mass(r, sol):
    plt.plot(r, sol[:, 1], 'g') 
    #plt.legend(loc='best')
    plt.xlabel('r \ [km]')
    plt.ylabel('$M(r) \ [M_{odot}]$')
    plt.grid()

def Plot_pressure(r, sol):    
    plt.plot(r, sol[:, 0], 'r')
    plt.legend(loc='best')
    plt.xlabel('$r \ [km]$')
    plt.ylabel('$P(r) \ [P_c]$')
    plt.grid()
    
def Plot_density(r, sol, Equation_Of_State, N_data):
    dener = np.zeros(N_data+1)
    #ec = 7.531304e+35
    for i in np.arange(N_data):
        r[i] = "{:.2f}".format(r[i])
        dener[i] = "{:.6f}".format(Equation_Of_State(sol[i,0]))
        sol[i] = "{:.6f}".format(sol[i,0])
        print(r[i], sol[i,0], dener[i])

    """
    plt.title('Density in the NS.')
    plt.plot(r, dener, 'c', label='Newton case, Relativistic EOS')
    plt.legend(loc='best')
    plt.xlabel('$r [km]$')
    plt.ylabel('$\epsilon \ [\epsilon_0]$')
    plt.grid()
    plt.savefig('Plots/density_n_r.pdf')
    """
    
#Results
def main():
    #Initial conditions
    y0 = [Pc[comp], 0]
    
    #Plotting arrays
    rmin, rmax, dr = 0, 12.8, 0.1 
    rad_array = np.arange(rmin, rmax+dr, dr)
    Equation_Of_State = Equation_Of_State_R
    
    N = int((rmax - rmin)/dr)
    
    solution = ode(Newton, y0, rad_array, args=(Equation_Of_State,))
    
    
    
    
    
    
    #Plots
    Plot_density(rad_array, solution, Equation_Of_State, N)
    #print(Equation_Of_State(0))
    """
    
    
    
    
    dener = np.zeros(N)
    densup = 7.1727e-14
    R = 0
    
    
    for i in np.arange(N):
        dener[i] = Equation_Of_State_R(solution[i,0])
        R = rad_array[i] + dr
        if (dener[i] < densup):
            break
    
    r = np.arange(0, R, dr)
    dener = np.resize(dener, r.size)    
    
    print( 'Size or the arrays: ', r.size, dener.size )
    print( 'Radius: ', R)
    print( 'Density at the surface:' , dener[-1])

    """

#Output

main()







"""
//////////////////////////////////////////////////////////////////////////
Created on Fri Jan  1 13:32:11 2021                                    //
                                                                      //
@author: Nicolás Gómez                                               //
//////////////////////////////////////////////////////////////////////
"""

