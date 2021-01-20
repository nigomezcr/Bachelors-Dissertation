"""
Description:
"""
import numpy as np
import scipy.special as spy
import matplotlib.pyplot as plt

R1, R2, R3, R4 = 0.5, 1.0, 1.5, 2.0 

#effective potential
def V(x, R):
    E = 1 + 0.5*R*R
    return 0.5*E*x**2 - 0.125*x**4 
"""
x1 = np.arange(-0.5,0.51,0.01)
x2 = np.arange(-1  ,1.01,0.01)
x3 = np.arange(-1.5,1.51,0.01)
x4 = np.arange(-2  ,2.01,0.01)
y1 = V(x1, R1) - V(R1, R1)
y2 = V(x2, R2) - V(R2, R2)
y3 = V(x3, R3) - V(R3, R3)
y4 = V(x4, R4) - V(R4, R4)


plt.xlabel('x')
plt.ylabel('V_eff(x)')
plt.plot(x1, y1, 'b', label='R = 0.5')
plt.plot(x2, y2, 'r', label='R = 1')
plt.plot(x3, y3, 'g', label='R = 1.5')
plt.plot(x4, y4, 'y', label='R = 2')
plt.legend()

plt.savefig("Plots/effecive_potential.pdf")
"""

#Proper Trajectory of special relativity grav. train
def x(t, R):
    Rp = np.sqrt(4 + R*R)
    sn = spy.ellipj(np.pi*Rp*t , R/Rp )
    return sn[0]

"""
tau1 = np.arange(0,1.01,0.01)
tau2 = np.arange(0,1.01,0.01)
tau3 = np.arange(0,0.74,0.01)
tau4 = np.arange(0,0.306,0.01)
tau5 = np.arange(0,0.11,0.01)
Trajec1 = x(tau1, 0.000)
Trajec2 = x(tau2, 1)
Trajec3 = x(tau3, 4)
Trajec4 = x(tau4, 16)
Trajec5 = x(tau5, 64)


plt.xlabel('tau')
plt.ylabel('x (tau)')


plt.plot(tau1, Trajec1, 'b', label='R -> 0')
plt.plot(tau2, Trajec2, 'y', label='R = 1')
plt.plot(tau3, Trajec3, 'g', label='R = 4')
plt.plot(tau4, Trajec4, 'r', label='R = 16')
plt.plot(tau5, Trajec5, 'c', label='R = 64')
plt.legend()

plt.savefig("Plots/relativistic_trajectories.pdf")

"""

#Coordinate Trajectory 

def t(x, R):
    Rp = np.sqrt(4 + R*R)
    E = spy.ellipeinc(x/R, R/Rp)
    F = spy.ellipkinc(x/R, R/Rp)
    return Rp*E - 2/Rp*F

x = np.arange(0, 10, 0.01)
Traj = t(x, 4)/(2*np.pi)


plt.xlabel('t')
plt.ylabel('x (t)')

plt.plot(Traj, x)


"""
//////////////////////////////////////////////////////////////////////////
Created on Fri Jan  1 13:32:11 2021                                    //
                                                                      //
@author: Nicolás Gómez                                               //
//////////////////////////////////////////////////////////////////////
"""

