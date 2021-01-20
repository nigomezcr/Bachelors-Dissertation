"""
Description:
"""
import numpy as np
import matplotlib.pyplot as plt

def F(E, T):
    return 1/(np.exp((E-1)/T) + 1)

E1 = np.arange(0, 3, 0.05)
F1 = F(E1, 1)
F2 = F(E1, 0.1)
F3 = F(E1, 0.01)

plt.xlabel("$\epsilon / \mu$'")
plt.ylabel("$F(E,T)$")
plt.grid(True, 'major', 'both')

plt.plot(E1, F1, 'b', label='$k_B T / \mu = 1  $')
plt.plot(E1, F2, 'g', label='$k_B T / \mu = 0.1$')
plt.plot(E1, F3, 'r', label='$k_B T / \mu = 0.01$')

plt.legend()

plt.savefig('Plots/fermi_distribution.pdf')









"""
//////////////////////////////////////////////////////////////////////////
Created on Fri Jan  1 13:32:11 2021                                    //
                                                                      //
@author: Nicolás Gómez                                               //
//////////////////////////////////////////////////////////////////////
"""

