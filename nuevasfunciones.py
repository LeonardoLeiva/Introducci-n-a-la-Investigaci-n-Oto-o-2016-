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


def w_param(w_0, w_a, z, modelo=0):
    if modelo == 0:
        W = - 1
    elif modelo == 1:
        W = w_0 + w_a * z / (z + 1)
    return W


def fwexp(p, z, modelo=0):
    F = np.exp(3 * quad(f_i1, 0, z, args=(p, modelo, ))[0])
    return F


def f_i1(p, z, modelo=0):
    if modelo == 0:
        G = 0
    elif modelo == 1:
        G = (1 + w_param(p[2], p[3], modelo, z)) / (1 + z)
    return G


def hub(p, z, modelo=0):
    d_m = p[0]
    d_de = p[1]
    H = np.sqrt(np.absolute(d_m * (1 + z) ** 3 + d_de * fwexp(p, z, modelo) + (1 - d_m - d_de) * (1 + z) ** 2))
    return H


def EDO(z, DL, p, modelo=0):
    d_m = p[0]
    d_de = p[1]
    E = np.zeros(len(z))
    for i in range(len(z)):
        E[i] = np.sqrt(1 - DL ** 2 * np.absolute(1 - d_m - d_de)) / hub(p, z[i], modelo)
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


def res_EDO(p, DL_0, z_f, z_0=0., paso=1000, modelo=0):
    init = DL_0
    z = np.linspace(z_0, z_f, paso)
    sol = odeint(EDO, init, z, args=(p, modelo, ))
    return sol[:, 0]


def Distancia(p, z, modelo=0, paso=100, DL_0=0):
    d = []
    for i in range(z.size):
        if z.size == 1:
            di = res_EDO(p, DL_0, z, 0, paso, modelo)[paso - 1]
        else:
            di = res_EDO(p, DL_0, z[i], 0, paso, modelo)[paso - 1]
        d.append(di)
    return np.asarray(d)


def D_L(p, z, modelo=0, paso=100, DL_0=0):
    '''
    d = []
    for i in range(z.size):
        if z.size == 1:
            di = res_EDO(p, DL_0, z, 0, paso, modelo)[paso - 1]
            dii = di * (z + 1)
        else:
            di = res_EDO(p, DL_0, z[i], 0, paso, modelo)[paso - 1]
            dii = di * (z[i] + 1)
        d.append(di)
    '''
    d_l = Distancia(p, z, modelo) * (1 + z)
    return np.asarray(d_l)


def mu_th(r, z, modelo=0):
    if modelo == 0:
        mu_0 = r[4]
        p = r[0], r[1], r[2], r[3]
    else:
        mu_0 = r[6]
        p = r[0], r[1], r[2], r[3], r[4], r[5]
    mu = 5 * np.log10(D_L(p, z, modelo))
    return mu
    '''
    mu = []
    for n in range(z.size):
        if z.size == 1:
            m = 5 * np.log10(D_L(r, z, modelo)) + mu_0
        else:
            m = 5 * np.log10(D_L(r, z[n], modelo)) + mu_0
        mu.append(m)
    return mu
    '''


# inicializacion
DL_0 = 0.
z_0 = 0.
d_m0 = 0.32
d_de0 = 0.68
w_00 = - 1.
w_a0 = 0.
z_f = 1.
mu_0 = 43.17
z_f0 = np.linspace(0.1, 2, 20)
p = d_m0, d_de0, w_00, w_a0
r = d_m0, d_de0, w_00, w_a0, mu_0
#print res_EDO(p, DL_0, z_f, z_0, modelo=0)
D = Distancia(p, z_f0)
mu = mu_th(r, z_f0)
print z_f0, mu
