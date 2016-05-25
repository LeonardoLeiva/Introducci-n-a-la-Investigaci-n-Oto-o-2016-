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


def leer_archivo(nombre):
    '''
    lee el archivo
    nombre debe ser un str
    '''
    datos = np.loadtxt(nombre, usecols=(1, 2, 3))
    z = datos[:, 0]
    mu = datos[:, 1]
    err_mu = datos[:, 2]
    return z, mu, err_mu


def mu_th(p, z):
    '''
    modulo de la distancia dependiente del redshift
    '''
    omega_m, omega_de, mu_00 = p
    N = z.size
    mu_0 = mu_00 * np.ones(N)
    p = 5 * np.log10(d_l(omega_m, omega_de, z))
    mu = 5 * np.log10(d_l(omega_m, omega_de, z)) + mu_0
    return mu


def d_l(omega_m, omega_de, z):
    '''
    distancia de luminosidad
    '''
    '''
    d = (1 + z) * integral(omega_m, omega_de, z)
    return d
    d = np.zeros(z.size)
    for i in range(z.size):
        d[i] = (1 + z[i]) * integral(omega_m, omega_de, z[i])
    return d
    '''
    dist = []
    for i in range(0, z.size):
        if z.size == 1:
            da = (1 + z) * integral(omega_m, omega_de, z)
        else:
            da = (1 + z[i]) * integral(omega_m, omega_de, z[i])
        dist.append(da)
    return dist


def integral(omega_m, omega_de, z):
    '''
    integral que incluye el parametro de hubble
    '''
    #print z
    I = quad(f_integrar, 0, z, args=(omega_m, omega_de)) # revisar args!!!!
    #print I
    return I[0]


def f_integrar(z, omega_m, omega_de):
    '''
    parametro de hubble normalizado por la constante de hubble
    '''
    h = omega_m * (z + 1) ** 3 + (1 - omega_m)
    #h = omega_m * (z + 1) ** 3 + omega_de + (z + 1) ** 2 * (1 - omega_m - omega_de)
    #print h
    return 1 / np.sqrt(h)


def residuo_modelo(p, z_exp, mu_exp):
    '''
    diferencia entre los valores del modelo y los experimentales
    '''
    #print mu_exp.shape
    #print mu_th(p, z_exp)
    err = mu_exp - mu_th(p, z_exp)
    return err


def optimizar(err, p0, z_exp, mu_exp):
    '''
    calcula los mejores parametros que se ajustan al modelo
    '''
    opt = leastsq(err, p0, args=(z_exp, mu_exp))
    return opt


def chi_cuadrado(p, x, y, f):
    S = np.sum((y - f(p, x)) ** 2)
    return S


#inicializacion
p0 = 0.3, 0.7, 43
z_exp, mu_exp, err_mu = leer_archivo('SnIa_data.txt')
b = residuo_modelo(p0, z_exp[0], mu_exp[0])
a = optimizar(residuo_modelo, p0, z_exp, mu_exp)
print a[0]
print chi_cuadrado(a[0], z_exp, mu_exp, mu_th), chi_cuadrado(p0, z_exp, mu_exp, mu_th)
np.save('resultados1.txt', a)
