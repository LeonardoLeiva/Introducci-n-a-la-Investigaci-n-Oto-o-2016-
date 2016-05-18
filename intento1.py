#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
primer intento de chi cuadrado para los datos experimentales de supernovas
'''

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
import scipy.stats
from scipy.optimize import leastsq
from scipy import optimize as opt

# funciones estructurales


def leer_archivo(nombre):
    '''
    lee el archivo
    nombre debe ser un str
    '''
    datos = np.loadtxt(nombre)
    z = datos[:, 1]
    mu = datos[:, 2]
    err_mu = datos[:, 3]
    return z, mu, err_mu


def mu_th(omega_m, omega_de, mu_0, z):
    '''
    modulo de la distancia dependiente del redshift
    '''
    mu = 5 * np.log10(D_l(omega_m, omega_de, z)) + mu_0
    return mu


def d_l(omega_m, omega_de, z):
    '''
    distancia de luminosidad
    '''
    D = (1 + z) * integral(omega_m, omega_de, z)
    return D


def integral(omega_m, omega_de, z):
    '''
    integral que incluye el parametro de hubble
    '''
    I = integrate.quad(f_integrar, 0, z, args=(omega_m, omega_de)) # revisar args!!!!
    return I


def f_integrar(omega_m, omega_de, z):
    '''
    parametro de hubble normalizado por la constante de hubble
    '''
    h = omega_m * (z + 1) ** 3 + omega_de + (z + 1) ** 2 * (1 - omega_m - omega_de)
    return 1 / h


def residuo_modelo(omega_m, omega_de, z_exp, mu_exp):
    '''
    diferencia entre los valores del modelo y los experimentales
    '''
    err = mu_exp - mu_th(omega_m, omega_de, z_exp)
    pass


def optimizar():
    '''
    calcula los mejores parametros que se ajustan al modelo
    '''
    pass


def chi_cuadrado(p, x, y, f):
    S = np.sum((y - f(p, x)) ** 2)
    return S
