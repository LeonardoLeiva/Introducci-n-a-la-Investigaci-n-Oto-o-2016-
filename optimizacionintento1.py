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


def hub(p, z, modelo=0):
    d_m = p[0]
    d_de = p[1]
    H = np.sqrt(np.absolute(d_m * (1 + z) ** 3 + d_de * + (1 - d_m - d_de) * (1 + z) ** 2))
    return H


def EDO(z, DL, p, modelo=0):
    d_m = p[0]
    d_de = p[1]
    E = np.zeros(len(z))
    for i in range(len(z)):
        E[i] = np.sqrt(np.absolute(1 - DL ** 2 * np.absolute(1 - d_m - d_de))) / hub(p, z[i], modelo)
    return E


def res_EDO(p, DL_0, z_f, z_0=0., paso=100, modelo=0):
    init = DL_0
    z = np.linspace(z_0, z_f, paso)
    sol = odeint(EDO, init, z, args=(p, modelo, ))
    return sol[:, 0]


def Distancia(p, z, modelo=0, paso=100, DL_0=0):
    #a = time.time()-t0
    d = []
    for i in range(z.size):
        if z.size == 1:
            di = res_EDO(p, DL_0, z, 0, paso, modelo)[paso - 1]
        else:
            di = res_EDO(p, DL_0, z[i], 0, paso, modelo)[paso - 1]
        d.append(di)
    #b = time.time()-t0
    #print("EDO: "+str(b-a))
    return np.asarray(d)


def D_L(p, z, modelo=0, paso=100, DL_0=0):
    d_l = Distancia(p, z, modelo) * (1 + z)
    return np.asarray(d_l)


def mu_th(r, z, modelo=0):
    p = r[0], r[1]
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


def prior(beta, p, model=0):
    '''
    probabilidad a priori que se le da a la adivinanza inicial. dos modelos
    '''
    d_m, d_de = beta
    a0, a1, b0, b1 = p
    s0 = ((d_m - a0) / b0) ** 2
    s1 = ((d_de - a1) / b1) ** 2
    s = (s0 + s1) / 2.
    P = np.exp(-s) / (2 * np.pi * b0 * b1)
    return P


def likelihood(beta, datos, error, model=0):
    '''
    verosimilitud. se aborda en una misma funcion los dos modelos
    '''
    x, y = datos
    N = len(x)
    s = np.sum((y - mu_th(beta, x, modelo=0)) ** 2)
    L = (2 * np.pi * error ** 2) ** (-N / 2) * np.exp(-s /
                                                      (2 * error ** 2))
    return L, s


def paso_metropolis(p0, prior_params, datos, L_0, a, error=1, d=0.05, modelo=0):
    mu_0 = marginalizar_mu0(p0, datos)
    dats = datos[0], datos[1] - mu_0
    x0, y0 = p0
    rx = np.random.uniform(low=-1, high=1)
    ry = np.random.uniform(low=-1, high=1)
    xp = x0 + d * rx
    yp = y0 + d * ry
    while xp < 0 or xp > 1:
        rx = np.random.uniform(low=-1, high=1)
        xp = x0 + d * rx
        print("densidad masa: "+str(xp))
    while yp < 0 or yp > 1:
        ry = np.random.uniform(low=-1, high=1)
        yp = y0 + d * ry
        #print("densidad energia oscura: "+str(yp))
    #L_0 = likelihood([x0, y0], dats, error)
    L_p = likelihood([xp, yp], dats, error)
    posterior_p0 = prior([x0, y0], prior_params) * L_0[0]
    posterior_pp = prior([xp, yp], prior_params) * L_p[0]
    if posterior_p0 > posterior_pp:
        P = posterior_pp / posterior_p0
        R = np.random.uniform(0, 1)
        if P > R:
            p0 = [xp, yp]
            a += 1
            L_0 = L_p
            #print "se acepta"
    else:
        p0 = [xp, yp]
        a += 1
        L_0 = L_p
        #print "se acepta"
    return p0, L_0, a


def monte_carlo(p0, prior_params, N, datos, error=1, d=0.05, modelo=0):
    a = 0
    muestra_met = np.zeros((N, 2))
    chi_cuad = np.zeros(N)
    muestra_met[0] = [p0[0], p0[1]]
    mu_0 = marginalizar_mu0(p0, datos)
    dats = datos[0], datos[1] - mu_0
    L_0 = likelihood(p0, dats, error)
    #posterior_p0 = prior(p0, prior_params) * L_0[0]
    chi_cuad[0] = L_0[1]
    for i in range(1, N):
        P_M = paso_metropolis(muestra_met[i-1], prior_params, datos, L_0, a, error, d, modelo)
        muestra_met[i] = P_M[0]
        L_0 = P_M[1]
        chi_cuad[i] = L_0[1]
        a = P_M[2]
        #posterior_p0 = P_M[4]
        print("contador: "+str(i))
    # guardar datos
    np.save('a.npy', muestra_met)
    np.save('b.npy', chi_cuad)
    np.save('c.npy', a)
    return muestra_met, chi_cuad, a


def marginalizar_mu0(r, datos, modelo=0):
    z, m_obs, error_obs = datos
    B = np.sum((m_obs - mu_th(r, z, modelo)) / error_obs ** 2)
    C = np.sum(1 / error_obs ** 2)
    mu_0 = B / C
    return mu_0


# inicializacion
z, mu, mu_err = leer_archivo('SnIa_data.txt')
datos = z, mu, mu_err
adivinanza1 = [0.2, 0.3, 0.8, 0.5]
N = 5000
p0 = 0.99999, 0.99999
t0 = time.time()
resultados = monte_carlo(p0, adivinanza1, N, datos)
tf=time.time()-t0
print("tiempo: "+str(tf))
A = np.load('a.npy')
B = np.load('b.npy')
C = np.load('c.npy')
# grafico
fig = plt.figure()
fig.clf()
ax1 = fig.add_subplot(111)
ax1.plot(A[:,0], A[:,1], '.')
ax1.set_xlabel("Densidad de Materia")
ax1.set_ylabel("Densidad de Energia Oscura")
plt.legend(loc=4)
plt.savefig("mcmc.png")
plt.draw()
plt.show()
