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

# funciones estructurales


def mu_th(p, z):
    '''
    modulo de la distancia dependiente del redshift
    '''
    omega_m, omega_de, mu_00 = p
    N = z.size
    mu_0 = mu_00 * np.ones(N)
    '''
    distancia de luminosidad
    '''
    dist = []
    I = quad(f_integrar, 0, z, args=(omega_m, omega_de))
    for i in range(0, z.size):
        if z.size == 1:
            da = (1 + z) * I[0]
        else:
            da = (1 + z[i]) * integral(omega_m, omega_de, z[i])
        dist.append(da)
    p = 5 * np.log10(dist)
    mu = 5 * np.log10(dist) + mu_0
    return mu


def f_integrar(z, omega_m, omega_de):
    '''
    parametro de hubble normalizado por la constante de hubble
    '''
    h = omega_m * (z + 1) ** 3 + (1 - omega_m)
    #h = omega_m * (z + 1) ** 3 + omega_de + (z + 1) ** 2 * (1 - omega_m - omega_de)
    #print h
    return 1 / np.sqrt(h)


#
