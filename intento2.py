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
    mu = 5 * np.log10(D_l(r, z)) + mu_0
    return mu


def d_l(r, z, e=0.0000001):
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
    omega_k = np.absolute(1 - omega_m - omega_de)
    for i in range(0, z.size):
        if z.size == 1:
            if (omega_k > -e) and (omega_k < e):
                da = (1 + z) * integral(r, z)
            else:
                da = (1 + z) * np.sin(np.sqrt(omega_k) * integral(r, z)) / np.sqrt(omega_k)
        else:
            if (omega_k > -e) and (omega_k < e):
                da = (1 + z[i]) * integral(r, z[i])
            else:
                da = (1 + z[i]) * np.sin(np.sqrt(omega_k) * integral(r, z[i])) / np.sqrt(omega_k)
        dist.append(da)
    return dist


def D_l(r, z, e=0.00001):
    omega_m, omega_de= r
    #omega_m, omega_de, o= r
    omega_k = 1 - omega_m - omega_de
    dist = []
    for i in range(0, z.size):
        if omega_k > e:
            omega_k = np.absolute(omega_k)
            if z.size == 1:
                da = (1 + z) * np.sinh(np.sqrt(omega_k) * integral(r, z)) / np.sqrt(omega_k)
            else:
                da = (1 + z[i]) * np.sinh(np.sqrt(omega_k) * integral(r, z[i])) / np.sqrt(omega_k)
        elif omega_k < -e:
            omega_k = np.absolute(omega_k)
            if z.size == 1:
                da = (1 + z) * np.sin(np.sqrt(omega_k) * integral(r, z)) / np.sqrt(omega_k)
            else:
                da = (1 + z[i]) * np.sin(np.sqrt(omega_k) * integral(r, z[i])) / np.sqrt(omega_k)
        else:
            if z.size == 1:
                da = (1 + z) * integral(r, z)
            else:
                da = (1 + z[i]) * integral(r, z[i])
        dist.append(da)
    return dist


def integral(r, z):
    '''
    integral que incluye el parametro de hubble
    '''
    #print z
    I = quad(f_integrar, 0, z, args=(r,), limit=50) # revisar args!!!!
    #print I
    return I[0]


def f_integrar(z, r):
    '''
    parametro de hubble normalizado por la constante de hubble
    '''
    omega_m, omega_de = r
    '''
    omega_m, omega_de, o = r
    if o == 0
    '''
    h = omega_m * (z + 1) ** 3 + omega_de + (1 - omega_m - omega_de) * (z + 1) ** 2
    return 1 / np.sqrt(np.absolute(h))


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


def graficar_varios(z, mu, mu_err, p):
    '''
    solo se agregan barras de errores
    '''
    mu_nuevo = mu_th(p, z)
    fig = plt.figure()
    fig.clf()
    ax1 = fig.add_subplot(111)
    ax1.errorbar(z, mu, mu_err, xerr=None,
         fmt='.', ecolor='red', elinewidth=None, capsize=None,
         barsabove=False, lolims=False, uplims=False,
         xlolims=False, xuplims=False, errorevery=1,
         capthick=None)
    ax1.plot(z, mu_nuevo, '+')
    ax1.set_xlabel("Redshift")
    ax1.set_ylabel("Modulo de la Distancia")
    plt.legend(loc=4)
    plt.savefig("prueba.png")
    plt.draw()
    plt.show()


def leer_archivo(nombre):
    '''
    lee el archivo
    nombre debe ser un str
    '''
    #datos = np.loadtxt(nombre, dtype=[('f0',str),('f1',float),('f2',float),('f3',float)])
    datos = np.loadtxt(nombre, usecols=(1, 2, 3))
    z = datos[:, 0]
    mu = datos[:, 1]
    err_mu = datos[:, 2]
    return z, mu, err_mu


#inicializacion
p0 = 0.32, 0.68, 43
z_0 = 0.1
p1 = 0.32, 0.68
p2 = 0.267, 0.732, 43.175
zn = np.linspace(z_0, 2, 20)
z_exp, mu_exp, err_mu = leer_archivo('SnIa_data.txt')
b = residuo_modelo(p0, z_exp[0], mu_exp[0])
a = optimizar(residuo_modelo, p0, z_exp, mu_exp)
print a[0]
print chi_cuadrado(a[0], z_exp, mu_exp, mu_th), chi_cuadrado(p0, z_exp, mu_exp, mu_th)
print zn
print D_l(p1, zn)
graficar_varios(z_exp, mu_exp, err_mu, p0)
