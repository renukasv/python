#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 09:46:23 2019

@author: renuka
"""

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as axes3d

fig = plt.figure()
ax = fig.gca(projection='3d')
phi = np.linspace(0, np.pi / 4, 100)
theta = np.linspace(4 * np.pi, 4 * np.pi, 100)
z = r * np.cos(theta)
r = 4
x = r * np.sin(theta) * np.cos(phi)
y = r * np.sin(theta) * sin(phi)
ax.plot(x, y, z, label='parametric curve')
ax.legend()
plt.show()