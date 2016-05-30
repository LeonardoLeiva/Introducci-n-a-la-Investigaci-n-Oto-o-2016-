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
    '''
    if len(p) == 5:
        omega_m, omega_de, mu_00, w_0, w_a = p
        r = omega_m, omega_de, w_0, w_a
    else:
        omega_m, omega_de, mu_00 = p
        r = omega_m, omega_de
    '''
    omega_m, omega_de, mu_00 = p
    r = omega_m, omega_de
    N = z.size
    mu_0 = mu_00 * np.ones(N)
    #p = 5 * np.log10(d_l(omega_m, omega_de, z))
    mu = 5 * np.log10(d_l(r, z)) + mu_0
    return mu


def d_l(r, z):
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
    '''
    if len(r) == 4:
        omega_m, omega_de, w_0, w_a = r
    else:
        omega_m, omega_de = r
    '''
    omega_m, omega_de = r
    omega_k = 1 - omega_m - omega_de
    for i in range(0, z.size):
        if z.size == 1:
            da = (1 + z) * np.sin(np.sqrt(omega_k) * integral(r, z)) / np.sqrt(omega_k)
        else:
            da = (1 + z[i]) * np.sin(np.sqrt(omega_k) * integral(r, z[i])) / np.sqrt(omega_k)
        dist.append(da)
    return dist


def integral(r, z):
    '''
    integral que incluye el parametro de hubble
    '''
    #print z
    I = quad(f_integrar, 0, z, args=(r,)) # revisar args!!!!
    #print I
    return I[0]


def f_integrar(z, r):
    '''
    parametro de hubble normalizado por la constante de hubble
    '''
    omega_m, omega_de = r
    h = omega_m * (z + 1) ** 3 + omega_de + (1 - omega_m - omega_de) * (z + 1) ** 2
    return 1 / np.sqrt(h)


def fwexp(omega_m, omega_de, w_0, w_a, z):
    F = np.exp(3 * quad(f_i1, 0, z, args=(d_m, d_de))[0])
    return F

def f_i1(w_0, w_a, z):
    G = (1 + w_param(w_0, w_a, z)) / (1 + z)
    return G


def w_param(w_0, w_a, z):
    W = w_0 + w_a * z / (z + 1)
    return W


def residuo_modelo(p, z_exp, mu_exp):
    '''
    diferencia entre los valores del modelo y los experimentales
    '''
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
