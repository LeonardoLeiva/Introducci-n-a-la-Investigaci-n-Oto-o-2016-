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


def graficar_errores(z, mu, mu_err):
    '''
    solo se agregan barras de errores
    '''
    fig = plt.figure()
    fig.clf()
    ax1 = fig.add_subplot(111)
    ax1.errorbar(z, mu, mu_err, xerr=None,
         fmt='.', ecolor='red', elinewidth=None, capsize=None,
         barsabove=False, lolims=False, uplims=False,
         xlolims=False, xuplims=False, errorevery=1,
         capthick=None)
    ax1.set_xlabel("Redshift")
    ax1.set_ylabel("Modulo de la Distancia")
    plt.legend(loc=4)
    plt.savefig("grafico2.png")
    plt.draw()
    plt.show()


#inicializacion
z_exp, mu_exp, mu_err = leer_archivo('SnIa_data.txt')
graficar_errores(z_exp, mu_exp, mu_err)
