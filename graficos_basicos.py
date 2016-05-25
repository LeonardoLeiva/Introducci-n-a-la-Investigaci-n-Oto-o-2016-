#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
graficos o:
'''

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import scipy.stats
from scipy.optimize import leastsq
from scipy import optimize as opt


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


def graficar(z, mu):
    '''
    grafica los datos experimentales
    '''
    fig = plt.figure()
    fig.clf()
    ax1 = fig.add_subplot(111)
    ax1.plot(z, mu, '+')
    ax1.set_xlabel("Redshift")
    ax1.set_ylabel("Modulo de la Distancia")
    plt.legend(loc=4)
    plt.savefig("grafico1.png")
    plt.draw()
    plt.show()


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
    plt.savefig("grafico2.png")
    plt.draw()
    plt.show()


def graficar_teorica(p, z):
    '''
    grafica el mu teorico segun la aproximacion de chi cuadrado de intento1.py
    '''
    mu = mu_th(p, z)
    fig1 = plt.figure()
    fig1.clf()
    ax1 = fig1.add_subplot(111)
    ax1.plot(z, mu, '+')
    ax1.set_xlabel("X")
    ax1.set_ylabel("W(x)")
    plt.savefig("a.png")
    plt.draw()
    plt.show()


def leer_resultado(nombre):
    '''
    lee los resultados obtenidos en intento1.py
    nombre debe ser str
    '''
    datos = np.loadtxt(nombre, usecols=(1, 2, 3))
    d_m = datos[:, 0]
    d_de = datos[:, 1]
    mu_0 = datos[:, 2]
    return d_m, d_de, mu_0


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


#inicializacion
dm = 0.267759
#p_0 = 0.3, 0.7, 43
#p_0 = 0.16581, 0.9403, 43.15
p_0 = dm, 1 - dm, 43.175
z_exp, mu_exp, mu_err = leer_archivo('SnIa_data.txt')
graficar_varios(z_exp, mu_exp, mu_err, p_0)
#graficar_teorica(p_0, z_exp)
