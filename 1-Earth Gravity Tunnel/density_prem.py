"""
Description: Prem density plot
"""
import numpy as np
import matplotlib.pyplot as plt


density_prem = open("Numerical Data/prem-density.csv")
temp_data = []

#Separate data in two
for line in density_prem:
    x, y = line.split('\t')
    temp_data += [ (float(x), float(y)) ]
density_prem.close()
data_set = np.array(temp_data)

x = np.zeros(55)
y = np.zeros(55)
for i in np.arange(55):
   x[i] = data_set[i][0]
   y[i] = data_set[i][1]
   
"""
plt.xlabel('r[km]')
plt.ylabel('Density [kg/m3]')
plt.plot(x,y,'-')
"""

# Now, let's plot the traversal times as a function of density
# for the case of constant density:
    
def times(x):
    T0 = 42.18
    return T0/np.sqrt(x)

n = np.arange(0, 5, 0.01)
Tn = times(n)
T_Earth = times(1)
T_Mars = times(0.71)
T_Sun = times(0.255)
T_Star = times(1e5)

T_EN =0
T_BH =0

plt.xlabel("Density [in units of Earth density]")
plt.ylabel("T [min]")
plt.plot(n, Tn, 'k')
#plt.plot(1, T_Earth, 'b*', label='Earth')
#plt.plot(0.255, T_Sun, 'y*', label='Sun')
#plt.plot(0.71, T_Mars, 'r*', label='Mars')
plt.legend()

plt.savefig("Plots/2-time_vs_density.pdf")













"""
//////////////////////////////////////////////////////////////////////////
Created on Fri Jan  1 13:32:11 2021                                    //
                                                                      //
@author: Nicolás Gómez                                               //
//////////////////////////////////////////////////////////////////////
"""

