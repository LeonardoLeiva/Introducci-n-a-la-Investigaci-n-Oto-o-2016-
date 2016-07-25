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


A = np.load('d.npy')
B = np.load('e.npy')
C = np.load('f.npy')
A0 = np.load('d0.npy')
B0 = np.load('e0.npy')
C0 = np.load('f0.npy')
n = len(B)
print n
N = np.linspace(0, n, num=5000)
print len(N)
# grafico
fig = plt.figure()
fig.clf()
ax1 = fig.add_subplot(111)
ax1.plot(A[:, 0], A[:, 1])
ax1.plot(A0[:, 0], A0[:, 1])
ax1.set_xlabel("Densidad de Materia")
ax1.set_ylabel("Densidad de Energia Oscura")
plt.legend(loc=4)
plt.savefig("mcmc.png")
plt.draw()
plt.show()
print C, C0
