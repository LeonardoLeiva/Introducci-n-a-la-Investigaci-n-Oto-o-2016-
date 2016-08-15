#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import scipy.stats
from scipy.optimize import leastsq
from scipy import optimize as opt
from scipy.integrate import odeint
import time
from mpl_toolkits.mplot3d import Axes3D


A = np.load('densidades.npy')
B = np.load('chi_cuadrado.npy')
C = np.load('cantidadpasosaceptados.npy')
D = np.load('densidadesaceptadas.npy')
'''
A = np.load('d.npy')
B = np.load('e.npy')
C = np.load('f.npy')
D = np.load('g.npy')
'''
print A
print B
n = len(B)
N = np.linspace(1, n, num=n)
# grafico
fig = plt.figure()
fig.clf()
ax1 = fig.add_subplot(111)
ax1.plot(D[:, 0], D[:, 1], '.')
ax1.set_xlabel("Densidad de Materia")
ax1.set_ylabel("Densidad de Energia Oscura")
plt.legend(loc=4)
plt.savefig("mcmc.png")
plt.draw()
plt.show()
fig = plt.figure()
fig.clf()
ax2 = fig.add_subplot(111)
ax2.plot(N, B, '.')
ax2.set_xlabel("Densidad de Materia")
ax2.set_ylabel("Densidad de Energia Oscura")
plt.legend(loc=4)
plt.savefig("chicuad.png")
plt.draw()
plt.show()
print C

fig=plt.figure(1)
fig.clf()
ax=fig.add_subplot(111, projection='3d')
ax.set_aspect('equal')
ax.plot(D[:, 0], D[:, 1], D[:, 2], '.')
ax.set_xlabel('Densidad de Materia')
ax.set_ylabel('Densidad de Energia Oscura')
ax.set_zlabel('Chi Cuadrado')
plt.savefig("chicuadenelespacio.png")
plt.show()
