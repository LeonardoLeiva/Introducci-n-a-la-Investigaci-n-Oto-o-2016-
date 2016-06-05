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
# la semilla!!!
np.random.seed(8)

# funciones estructurales previas (no mcmc)


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
        E[i] = np.sqrt(np.absolute(1 - DL ** 2 * np.absolute(1 - d_m - d_de))) / hub(p, z[i], modelo)
    return E


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
    d_l = Distancia(p, z, modelo) * (1 + z)
    return np.asarray(d_l)


def mu_th(r, z, modelo=0):
    if modelo == 0:
        #mu_0 = r[4]
        p = r[0], r[1]
    elif modelo == 1:
        #mu_0 = r[6]
        p = r[0], r[1], r[2], r[3]
    mu = 5 * np.log10(D_L(p, z, modelo))
    return mu


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


def marg_mu_0():
    '''
    marginaliza mu_0
    '''
    pass


# funciones mcmc


def prior(beta, p, model=0):
    '''
    probabilidad a priori que se le da a la adivinanza inicial. dos modelos
    '''
    if model == 1:
        d_m, d_de, w_0, w_a = beta
        a0, a1, a2, a3, b0, b1, b2, b3 = p
        s0 = ((d_m - a0) / b0) ** 2
        s1 = ((d_de - a1) / b1) ** 2
        s2 = ((w_0 - a2) / b2) ** 2
        s3 = ((w_a - a3) / b3) ** 2
        s = (s0 + s1 + s2 + s3) / 2.
        P = np.exp(-s) / (4 * np.pi ** 2 * b0 * b1 * b2 * b3)
        return P
    elif model == 0:
        d_m, d_de = beta
        a0, a1, b0, b1 = p
        s0 = ((d_m - a0) / b0) ** 2
        s1 = ((d_de - a1) / b1) ** 2
        s = (s0 + s1) / 2.
        P = np.exp(-s) / (2 * np.pi * b0 * b1)
        return P



def fill_prior(beta_grid, prior_p, model=0):
    '''
    rellena la grilla para la prob a priori. dos modelos al
    mismo tiempo
    '''
    if model == 1:
        d_m_grid, d_de_grid, w_0_grid, w_a_grid = beta_grid
        salida = np.zeros(d_m_grid.shape)
        ni, nj, nk, nl = d_m_grid.shape
        for i in range(ni):
            for j in range(nj):
                for k in range(nk):
                    for l in range(nl):
                        salida[i, j, k, l] = prior([d_m_grid[i, j, k, l],
                                                    d_de_grid[i, j, k, l],
                                                    w_0_grid[i, j, k, l],
                                                    w_a_grid[i, j, k, l]],
                                                   prior_p, 1)
    elif model == 0:
        d_m_grid, d_de_grid = beta_grid
        salida = np.zeros(d_m_grid.shape)
        ni, nj = d_m_grid.shape
        for i in range(ni):
            for j in range(nj):
                salida[i, j] = prior([d_m_grid[i, j], d_de_grid[i, j]],
                                     prior_p, 0)
    return salida


def likelihood(beta, datos, error, model=0):
    '''
    verosimilitud. se aborda en una misma funcion los dos modelos
    '''
    x, y = datos
    N = len(x)
    if model == 1:
        s = np.sum(y - mu_th(beta, x, modelo=1))
        L = (2 * np.pi * error ** 2) ** (-N / 2.) * np.exp(- s /
                                                           (2 * error ** 2))
    elif model == 0:
        s = np.sum(y - mu_th(beta, x, modelo=0))
        L = (2 * np.pi * error ** 2) ** (-N / 2) * np.exp(-s /
                                                          (2 * error ** 2))
    return L


def fill_likelihood(beta_grid, datos, error, model=0):
    '''
    rellena la grilla para verosimilitud de los parametros.
    '''
    if model == 1:
        d_m_grid, d_de_grid, w_0_grid, w_a_grid = beta_grid
        sal = np.zeros(d_m_grid.shape)
        ni, nj, nk, nl = d_m_grid.shape
        for i in range(ni):
            for j in range(nj):
                for k in range(nk):
                    for l in range(nl):
                        sal[i, j, k, l] = likelihood([d_m_grid[i, j, k, l],
                                                      d_de_grid[i, j, k, l],
                                                      w_0_grid[i, j, k, l],
                                                      w_a_grid[i, j, k, l]],
                                                     datos, error, 1)
    elif model == 0:
        d_m_grid, d_de_grid = beta_grid
        sal = np.zeros(d_m_grid.shape)
        ni, nj = d_m_grid.shape
        for i in range(ni):
            for j in range(nj):
                sal[i, j] = likelihood([d_m_grid[i, j], d_de_grid[i, j]],
                                       datos, error, 0)
    return sal


def paso_metropolis(p0, prior_params, datos, error=1, d=0.05, modelo=0):
    c = chi_cuadrado(p0, datos[0], datos[1], mu_th)
    print c
    if modelo == 0:
        x0, y0 = p0
        rx = np.random.uniform(low=-1, high=1)
        ry = np.random.uniform(low=-1, high=1)
        xp = x0 + d * rx
        yp = y0 + d * ry
        posterior_p0 = prior([x0, y0], prior_params) * likelihood([x0, y0], datos, error)
        posterior_pp = prior([xp, yp], prior_params) * likelihood([xp, yp], datos, error)
        if (posterior_pp / posterior_p0) > np.random.uniform(0, 1):
            p0 = [xp, yp]
    elif modelo == 1:
        x0, y0, u0, v0 = p0
        rx = np.random.uniform(low=-1, high=1)
        ry = np.random.uniform(low=-1, high=1)
        ux = np.random.uniform(low=-1, high=1)
        vy = np.random.uniform(low=-1, high=1)
        xp = x0 + d * rx
        yp = y0 + d * ry
        up = u0 + d * ru
        vp = v0 + d * rv
        posterior_p0 = prior([x0, y0, u0, v0], prior_params) * likelihood([x0, y0, u0, v0], datos, error)
        posterior_pp = prior([xp, yp, up, vp], prior_params) * likelihood([xp, yp, up, vp], datos, error)
        P = posterior_pp / posterior_p0
        print P
        if P > np.random.uniform(0, 1):
            p0 = [xp, yp, up, vp]
    print p0
    return p0


def monte_carlo(p0, prior_params, N, datos, error=1, d=0.05, modelo=0):
    if modelo == 0:
        muestra_met = np.zeros((N, 2))
        muestra_met[0] = [p0[0], p0[1]]
        rechazados = 0
        for i in range(1, N):
            muestra_met[i] = paso_metropolis(muestra_met[i-1], prior_params, datos, error, d, modelo) # intentar 0.1, 0.5, 1., 3., 10.
            if muestra_met[i][0] == muestra_met[i-1][0]:
                rechazados += 1
            print N
    elif modelo == 1:
        muestra_met = np.zeros((N, 4))
        muestra_met[0] = [p0[0], p0[1], p0[2], p0[3]]
        rechazados = 0
        for i in range(1, N):
            muestra_met[i] = paso_metropolis(muestra_met[i-1], prior_params, datos, d, error, d, modelo) # intentar 0.1, 0.5, 1., 3., 10.
            if muestra_met[i][0] == muestra_met[i-1][0]:
                rechazados += 1
            print N
    return muestra_met, rechazados


def chi_cuadrado(p, x, y, f):
    S = np.sum((y - f(p, x)) ** 2)
    return S


# inicializacion
mu_0 = 43.15
z, mu, mu_err = leer_archivo('SnIa_data.txt')
datos = z, mu - mu_0
beta_grid1 = np.mgrid[0.:1.:50j, 0.:1.:50j]
d_m_grid, d_de_grid = beta_grid1
adivinanza1 = [0.2, 0.3, 0.8, 0.5]
N = 1000
p0 = 0.5, 0.5
resultados = monte_carlo(p0, adivinanza1, N, datos)
print resultados[1]
