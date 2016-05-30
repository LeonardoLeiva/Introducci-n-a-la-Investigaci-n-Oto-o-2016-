#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
primer intento de chi cuadrado para los datos experimentales de supernovas
'''

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import scipy.stats
from scipy.optimize import leastsq
from scipy import optimize as opt
from scipy.integrate import odeint

# funciones estructurales


def w_param(w_0, w_a, z):
    W = w_0 + w_a * z / (z + 1)
    return W


def fwexp(p, z):
    d_m, d_de, w_0, w_a = p
    '''
    if z == 0:
        F = 1
    else:
        F = np.exp(3 * quad(f_i1, 0, z, args=(d_m, d_de)))
    F = []
    for i in range(len(z)):
        f = np.exp(3 * quad(f_i1, 0, z[i], args=(d_m, d_de)))
        if len(F) == 1:
            F = [f]
        else:
            F.append(f)
    '''
    F = np.exp(3 * quad(f_i1, 0, z, args=(d_m, d_de))[0])
    return F


def f_i1(w_0, w_a, z):
    G = (1 + w_param(w_0, w_a, z)) / (1 + z)
    return G


def hub(p, z):
    d_m, d_de, w_0, w_a = p
    H = np.sqrt(np.absolute(d_m * (1 + z) ** 3 + d_de * fwexp(p, z) + (1 - d_m - d_de) * (1 + z) ** 2))
    return H


def EDO(z, DL, p):
    d_m, d_de, w_0, w_a = p
    E = np.zeros(len(z))
    for i in range(len(z)):
        E[i - 1] = np.sqrt(1 + DL ** 2 * np.absolute((1 - d_m - d_de))) / hub(p, z[i -1])
    #E = np.sqrt(1 + DL ** 2 * (1 - d_m - d_de)) / hub(p, z)
    return E


def resolucion_EDO(p, DL_0, z_0, z_f):
    r = ode(EDO)
    r.set_integrator('dopri853')
    r.set_initial_value(DL_0, z_0)
    r.set_f_params(p)
    z = np.linspace(z_0, z_f, 100)
    DL = []
    for i in range(len(z)):
        r.integrate(z[i])
        if i == 1:
            DL.insert(0, r.y)
        else:
            DL.append(r.y)
    return DL


def res_EDO(p, DL_0, dDL_0, z_f, z_0=0., paso=100):
    init = DL_0, dDL_0
    z = np.linspace(z_0, z_f, paso)
    sol = odeint(EDO, init, z, args=(p,))
    return [z, sol[:, 0], sol[:, 1]]


# inicializacion
DL_0 = 0.
dDL_0 = 1.
z_0 = 0.
d_m0 = 0.32
d_de0 = 0.679
w_00 = - 1.
w_a0 = 0.
z_f = 1.
#z_f = np.linspace(z_0, 2, 20)
p = d_m0, d_de0, w_00, w_a0
print res_EDO(p, DL_0, dDL_0, z_f)
