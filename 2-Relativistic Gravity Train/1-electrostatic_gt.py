"""
Description: 2.1 Python code to compute traversal times for electrostatic tunnel.
"""
import numpy as np
import matplotlib.pyplot as plt


k = 8.988e9 #In Nm2/C2
m = [1,100,200,400,800] #in kg
q = np.array([1,2,3,5]) #in C
Q = 10e8 #in C
R = 10e3#6.371e6 #in m

qk = np.sqrt(k)*q
Qk = np.sqrt(k)*Q
Rho_Q = 3*Q/(4*np.pi*R**3) #9.28 #in C/m3
Rho_Qk = np.sqrt(k)*Rho_Q

    
#Physical functions
def V(Q, q, m, R):
    return np.sqrt(2*Q*q/(m*R))


def T(rho, q, m):
    return np.sqrt(3*np.pi*m/(rho*q))



#Plotting functions
def Plot_mass():
    
    x = np.arange(1,10,0.1)
    Rho = Rho_Qk*x

    Time1 = T(Rho, qk[0], m[1]) #in s
    Time2 = T(Rho, qk[0], m[2]) #in s
    Time3 = T(Rho, qk[0], m[3]) #in s
    Time4 = T(Rho, qk[0], m[4]) #in s

    plt.xlabel(' 'r'$\rho[C/m^3]$')
    plt.ylabel('$T [s]$')
    plt.plot(x, Time1, 'b', label='m=100kg')
    plt.plot(x, Time2, 'r', label='m=200kg')
    plt.plot(x, Time3, 'y', label='m=400kg')
    plt.plot(x, Time4, 'g', label='m=800kg')
    plt.legend(fontsize='x-large')
    plt.savefig('Plots/1-Electrostatic_times_masses_py.pdf')



def Plot_charge():
    
    x = np.arange(1,10,0.1)
    Rho = Rho_Qk*x

    Time1 = T(Rho, qk[0], m[1]) #in s
    Time2 = T(Rho, qk[1], m[1]) #in s
    Time3 = T(Rho, qk[2], m[1]) #in s
    Time4 = T(Rho, qk[3], m[1]) #in s

    plt.xlabel('$\\rho \ [\\rho_c]$')
    plt.ylabel('$T \ [s]$')
    plt.plot(x, Time1, 'b', label='q=1 C')
    plt.plot(x, Time2, 'r', label='q=2 C')
    plt.plot(x, Time3, 'y', label='q=3 C')
    plt.plot(x, Time4, 'g', label='q=5 C')
    plt.legend(fontsize='x-large')
    plt.savefig('Plots/1-Electrostatic_times_charges_py.pdf')

    


def main():
    print('Density = ', Rho_Qk, 'C/m^3')
    
    Velocity = V(Qk, qk[0], m[0], R)
    print("v = ", Velocity, "m/s")
    print("  = ", Velocity/3e8, 'c')

    #Plot_mass()
    Plot_charge()


main()




"""
//////////////////////////////////////////////////////////////////////////
Created on Fri Jan  1 13:32:11 2021                                    //
                                                                      //
@author: Nicolás Gómez                                               //
//////////////////////////////////////////////////////////////////////
"""
