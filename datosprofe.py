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


def leer_archivo(nombre):
    '''
    lee el archivo
    nombre debe ser un str
    '''
    datos = np.loadtxt(nombre, usecols=(1, 2, 3, 4))
    d_m = datos[:, 0]
    d_de = datos[:, 1]
    no_se = datos[:, 2]
    chi_cuad = datos[:, 3]
    return z, mu, err_mu, chi_cuad

#A, B, C, D = leer_archivo('chi2_leonardo_table.txt')
A, B, C, D = np.load('chi2_leonardo_table.npy')
fig = plt.figure()
fig.clf()
ax1 = fig.add_subplot(111)
ax1.plot(A, B, '.')
ax1.set_xlabel("Densidad de Materia")
ax1.set_ylabel("Densidad de Energia Oscura")
plt.legend(loc=4)
plt.savefig("datosprofe.png")
plt.draw()
plt.show()
fig = plt.figure()
fig.clf()


D_min, D_max = -np.abs(D).max(), np.abs(D).max()
plt.subplot(1, 1, 1)
plt.imshow(D, cmap='RdBu', vmin=D_min, vmax=D_max,
           extent=[A.min(), A.max(), B.min(), B.max()],
           interpolation='nearest', origin='lower')
plt.title('image (interp. nearest)')
plt.colorbar()
plt.show()
