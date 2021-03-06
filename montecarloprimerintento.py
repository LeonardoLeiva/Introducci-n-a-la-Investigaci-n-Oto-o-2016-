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
import time
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
        s = np.sum((y - mu_th(beta, x, modelo=1)) ** 2)
        L = (2 * np.pi * error ** 2) ** (-N / 2.) * np.exp(- s /
                                                           (2 * error ** 2))
    elif model == 0:
        s = np.sum((y - mu_th(beta, x, modelo=0)) ** 2)
        L = (2 * np.pi * error ** 2) ** (-N / 2) * np.exp(-s /
                                                          (2 * error ** 2))
    return L, s


def paso_metropolis(p0, prior_params, datos, a, error=1, d=0.05, modelo=0):
    mu_0 = marginalizar_mu0(p0, datos)
    dats = datos[0], datos[1] - mu_0
    if modelo == 0:
        x0, y0 = p0
        rx = np.random.uniform(low=-1, high=1)
        ry = np.random.uniform(low=-1, high=1)
        xp = x0 + d * rx
        yp = y0 + d * ry
        while (xp < 0 or xp > 1) or (yp < 0 or yp > 1):
            rx = np.random.uniform(low=-1, high=1)
            ry = np.random.uniform(low=-1, high=1)
            xp = x0 + d * rx
            yp = y0 + d * ry
        L_0 = likelihood([x0, y0], dats, error)
        L_p = likelihood([xp, yp], dats, error)
        posterior_p0 = prior([x0, y0], prior_params) * L_0[0]
        posterior_pp = prior([xp, yp], prior_params) * L_p[0]
        P = posterior_pp / posterior_p0
        R = np.random.uniform(0, 1)
        if P > R:
            p0 = [xp, yp]
            a += 1
            print "se acepta"
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
        while xp < 0 or yp < 0 or up < 0 or vp < 0:
            rx = np.random.uniform(low=-1, high=1)
            ry = np.random.uniform(low=-1, high=1)
            ux = np.random.uniform(low=-1, high=1)
            vy = np.random.uniform(low=-1, high=1)
            xp = x0 + d * rx
            yp = y0 + d * ry
            up = u0 + d * ru
            vp = v0 + d * rv
        L_0 = likelihood([x0, y0, u0, v0], dats, error)
        L_p = likelihood([xp, yp, up, vp], dats, error)
        posterior_p0 = prior([x0, y0, u0, v0], prior_params) * L_0[0]
        posterior_pp = prior([xp, yp, up, vp], prior_params) * L_p[0]
        P = posterior_pp / posterior_p0
        R = np.random.uniform(0, 1)
        if P > R:
            p0 = [xp, yp, up, vp]
    return p0, L_0[1], a


def monte_carlo(p0, prior_params, N, datos, error=1, d=0.05, modelo=0):
    a = 0
    if modelo == 0:
        muestra_met = np.zeros((N, 2))
        chi_cuad = np.zeros(N)
        muestra_met[0] = [p0[0], p0[1]]
        for i in range(1, N):
            P_M = paso_metropolis(muestra_met[i-1], prior_params, datos, a, error, d, modelo)
            muestra_met[i] = P_M[0] # intentar 0.1, 0.5, 1., 3., 10.
            chi_cuad[i] = P_M[1]
            a = P_M[2]
            print("contador: "+str(i))
    elif modelo == 1:
        muestra_met = np.zeros((N, 4))
        chi_cuad = np.zeros((N, 2))
        muestra_met[0] = [p0[0], p0[1], p0[2], p0[3]]
        for i in range(1, N):
            P_M = paso_metropolis(muestra_met[i-1], prior_params, datos, d, error, d, modelo)
            muestra_met[i] = P_M[0] # intentar 0.1, 0.5, 1., 3., 10.
            chi_cuad[i-1] = P_M[1]
            a = P_M[2]
            print("contador: "+str(i))
    # guardar datos
    np.save('param_mcmc.npy', muestra_met)
    np.save('chi.npy', chi_cuad)
    np.save('rechazados.npy', a)
    '''
    plt.plot(muestra_met[:,0], muestra_met[:,1], marker='None', ls='-', lw=0.3, color='w')
    '''
    return muestra_met, chi_cuad, a


def chi_cuadrado(p, x, y, f):
    S = np.sum((y - f(p, x)) ** 2)
    return S


def marginalizar_mu0(r, datos, modelo=0):
    z, m_obs, error_obs = datos
    #A = np.sum((m_obs - mu_th(r, z, modelo)) ** 2 / error_obs ** 2)
    B = np.sum((m_obs - mu_th(r, z, modelo)) / error_obs ** 2)
    C = np.sum(1 / error_obs ** 2)
    mu_0 = B / C
    return mu_0


# inicializacion
mu_0 = 43.15
z, mu, mu_err = leer_archivo('SnIa_data.txt')
datos = z, mu, mu_err
'''
beta_grid1 = np.mgrid[0.:1.:50j, 0.:1.:50j]
d_m_grid, d_de_grid = beta_grid1
'''
adivinanza1 = [0.2, 0.3, 0.8, 0.5]
N = 5000
p0 = 0.999, 0.9999
t0 = time.time()
resultados = monte_carlo(p0, adivinanza1, N, datos)
tf = time.time()-t0
print("tiempo: "+str(tf))
A = np.load('param_mcmc.npy')
B = np.load('chi.npy')
C = np.load('rechazados.npy')
print("parametros: "+str(A))
print("chi cuadrados: "+str(B))
print("aceptados: "+str(C))
# grafico
fig = plt.figure()
fig.clf()
ax1 = fig.add_subplot(111)
#ax1.xlim(-0.2, 1.2)
#ax1.ylim(-0.2, 1.2)
ax1.plot(A[:,0], A[:,1], '.')
#ax1.plot(muestra_met[:,0], muestra_met[:,1], marker='None', ls='-', lw=0.3, color='w')
ax1.set_xlabel("Densidad de Materia")
ax1.set_ylabel("Densidad de Energia Oscura")
plt.legend(loc=4)
plt.savefig("mcmc.png")
plt.draw()
plt.show()
