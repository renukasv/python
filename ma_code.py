#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 17:02:57 2019

@author: renuka
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import math
from math import *
from mpl_toolkits.mplot3d import Axes3D

# initial conditions who won't change
ptheta = 0
rs = 2
t0 = 0
theta = np.pi/2


def func(f, T, rs):
    # f
    t, r, theta, phi, pt, pr, ptheta, pphi = f
    # derivatives
    dtdT = -((1-rs/r)**-1)*pt
    drdT = (1-rs/r)*pr
    dthetadT = ptheta/(r**2) 
    dphidT = pphi/((r**2)*(np.sin(theta)**2))
    dptdT = 0
    dprdT = ((ptheta**2)/(r**3)) + ((pphi**2)/(np.sin(theta)*(r**3))) - ((rs*(pr**2))/(2*(r**2))) - (((pt**2)*rs)/(2*(r**2)*((1-rs/r)**2)))
    dpthetadT = np.cos(theta)*(pphi**2)/((np.sin(theta)**3)*(r**2))
    dpphidT =  0
    return(dtdT, drdT, dthetadT, dphidT, dptdT, dprdT, dpthetadT, dpphidT)

# initial coordinates
r0 = 1.5*rs
r1 = 2*rs # 1.9 tourne puis va dedans
r2 = 10*rs
r3 = 30*rs
r4 = 10*rs
phi0 = 0
phi1 = np.pi/20
phi2 = np.pi/15
phi3 = 2*np.pi
phi4 = 3*np.pi
pr0 = 0 
pr1 = -0.3           
pr2 = -0.3
pr3 = -0.5
pr4 = -1
pphi0 = 4
pphi1 = 2 # 0.5 => photon goes right in the black hole, 2 => tourne autour 1 fois
pphi2 = 1
pphi3 = 2
pphi4 = 0.5

def time_momentum(r, pr, pphi):
    pt = -(sqrt((((1-(rs/r))**2)*(pr**2)) + (pphi**2)*(1-rs/r)/(r**2)))
    return(pt)
    
pt0 = time_momentum(r0, pr0, pphi0)
pt1 = time_momentum(r1, pr1, pphi1)
pt2 = time_momentum(r2, pr2, pphi2)
pt3 = time_momentum(r3, pr3, pphi3)
pt4 = time_momentum(r4, pr4, pphi4)

f0=[t0, r0, theta, phi0, pt0, pr0, ptheta, pphi0] # red
g0=[t0, r1, theta, phi1, pt1, pr1, ptheta, pphi1] # blue
h0=[t0, r2, theta, phi2, pt2, pr2, ptheta, pphi2] # green
i0=[t0, r3, theta, phi3, pt3, pr3, ptheta, pphi3] # yellow
j0=[t0, r4, theta, phi4, pt4, pr4, ptheta, pphi4] # turquoise

T = np.linspace(0, 1000, 9000)

xx=odeint(func, f0, T, args=(rs,)) # generate the solution, args=() is to add rs to the equations
yy=odeint(func, g0, T, args=(rs,))
zz=odeint(func, h0, T, args=(rs,))
ww=odeint(func, i0, T, args=(rs,))
vv=odeint(func, j0, T, args=(rs,))

# defining values from odeints
r00 = xx[:, 1]
theta00 = xx[:, 2]
phi00 = xx[:, 3]
pt00 = xx[:, 4]
pr00 = xx[:, 5]
pphi00 = xx[:, 7]

r11 = yy[:, 1]
theta11 = yy[:, 2]
phi11 = yy[:, 3]
pt11 = yy[:, 4]
pr11 = yy[:, 5]
pphi11 = yy[:, 7]

r22 = zz[:, 1]
theta22 = zz[:, 2]
phi22 = zz[:, 3]
pt22 = zz[:, 4]
pr22 = zz[:, 5]
pphi22 = zz[:, 7]


r33 = ww[:, 1]
theta33 = ww[:, 2]
phi33 = ww[:, 3]
pt33 = ww[:, 4]
pr33 = ww[:, 5]
pphi33 = ww[:, 7]


r44 = vv[:, 1]
theta44 = vv[:, 2]
phi44 = vv[:, 3]
pt44 = vv[:, 4]
pr44 = vv[:, 5]
pphi44 = vv[:, 7]

# checking H = 0
print("Hamiltonians")
    
def Hamiltonian(r, pt, pr, pphi):
    H = (-((1-rs/r)**-1)*(pt**2)/2 + (1-rs/r)*(pr**2)/2 + (pphi**2)/(2*(r**2)))
    if np.amax(H) < 10e-08:
        print("Your results are correct")
    else:
        print("Your results are wrong")
    return(H)    

def singularity(H, r):
    if (r).any < 1.5*rs:
        print(H)
    else:
        print("0")
        return(H, r)
        


    
#sing0=singularity(Hamiltonian(r00, pt00, pr00, pphi00), r00)        
#print(sing0)
        
print(Hamiltonian(r00, pt00, pr00, pphi00))   

print(Hamiltonian(r11, pt11, pr11, pphi11))   

print(Hamiltonian(r22, pt22, pr22, pphi22))   

print(Hamiltonian(r33, pt33, pr33, pphi33))

print(Hamiltonian(r44, pt44, pr44, pphi44))   
   
print("END")

# checking the thetas = pi/2
print("THE THETAS :")
print(theta00)
print(theta11)
print(theta22)
print(theta33)
print(theta44)
print("END")

# spherical to cartesian coordinates 
def sph2cart(r, phi, theta):
    X = r * np.cos(phi) * np.sin(theta)
    Y = r * np.sin(phi) * np.sin(theta)
    Z = r * np.cos(theta)
    return(X, Y, Z)

X0, Y0, Z0 = sph2cart(r00, phi00, theta00)
X1, Y1, Z1 = sph2cart(r11, phi11, theta11)
X2, Y2, Z2 = sph2cart(r22, phi22, theta22)
X3, Y3, Z3 = sph2cart(r33, phi33, theta33)
X4, Y4, Z4 = sph2cart(r44, phi44, theta44)


# Z values (checking if Z=0)
print("Z-Coordinates")
def zcoord(Z):
    if np.amax(Z) < 6e-8:
        print("Your results are correct")
    else:
       print("Your results are wrong") 
       return(Z)


print(zcoord(Z0))
print(zcoord(Z1))
print(zcoord(Z2))
print(zcoord(Z3))
print(zcoord(Z4))

print("END")

# r coordinates
print("R :")
print(r00)
print(r11)
print(r22)
print(r33)
print(r44)
print("END")

print("pt")
print(pt00)
print(pt11)
print(pt22)
print(pt33)
print(pt44)
print("END")
print("pphi")
print(pphi00)
print(pphi11)
print(pphi22)
print(pphi33)
print(pphi44)
print("END")

# black sphere
fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)


x = rs * np.outer(np.cos(u), np.sin(v))
y = rs * np.outer(np.sin(u), np.sin(v))
z = rs * np.outer(np.ones(np.size(u)), np.cos(v))

# plot black sphere
ax.plot_surface(x, y, z, rstride=4, cstride=4, color='k')

# plot trajectories
plt.plot(X0, Y0, Z0, 'r')
plt.plot(X1, Y1, Z1, 'b')
# plt.plot(X2, Y2, Z2, 'g')
plt.plot(X3, Y3, Z3, 'y')
plt.plot(X4, Y4, Z4, 'turquoise')

axlim = 5

# limits
ax.set_xlim(-axlim, axlim)
ax.set_ylim(-axlim, axlim)
ax.set_zlim(-axlim, axlim)
# labels
ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')

plt.show()