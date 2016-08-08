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
np.random.seed(888)


# funciones estructurales previas (no mcmc)


def hub(p, z, modelo=0):
    d_m = p[0]
    d_de = p[1]
    H = np.sqrt(np.absolute(d_m * (1 + z) ** 3 + d_de + (1 - d_m - d_de) * (1 + z) ** 2))
    return H


def EDO(z, DL, p, modelo=0):
    d_m = p[0]
    d_de = p[1]
    E = np.zeros(len(z))
    for i in range(len(z)):
        op = (1 + DL ** 2 * (1 - d_m - d_de))
        E[i] = np.sqrt(op) / hub(p, z[i], modelo)
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
    a0, b0, a1, b1 = p
    s0 = ((d_m - a0) / b0) ** 2
    s1 = ((d_de - a1) / b1) ** 2
    s = (s0 + s1) / 2.
    P = np.exp(-s) / (2 * np.pi * b0 * b1)
    return P


def likelihood(s, error, model=0):
    '''
    verosimilitud. se aborda en una misma funcion los dos modelos
    '''
    print s
    L = np.exp(-s / (2 * error ** 2))
    return L


def paso_metropolis(p0, pasos_aceptados, prior_params, datos, s_0, a, error=1, d=0.1, modelo=0):
    #mu_0 = marginalizar_mu0(p0, datos)
    #dats = datos[0], datos[1] - mu_0, datos[2]
    x0, y0 = p0
    rx = np.random.uniform(low=-1, high=1)
    ry = np.random.uniform(low=-1, high=1)
    xp = x0 + d * rx
    yp = y0 + d * ry
    while xp < 0 or xp > 1:
        rx = np.random.uniform(low=-1, high=1)
        xp = x0 + d * rx
        #print("densidad masa: "+str(xp))
    while yp < 0 or yp > 1:
        ry = np.random.uniform(low=-1, high=1)
        yp = y0 + d * ry
        #print("densidad energia oscura: "+str(yp))
    xi = xi_cuadrado([xp, yp], datos, mu_th)
    s_p = xi[0]
    if s_p > s_0:
        s = s_0 / s_p
        S = np.exp((s_0 - s_p) / 2)
        #S = np.exp(- s / 2)
        pp = prior([xp, yp], prior_params)
        p0 = prior([x0, y0], prior_params)
        L =  prior([xp, yp], prior_params) / prior([x0, y0], prior_params)
        #P = S * L
        P = S * L
        R = np.random.uniform(0, 1)
        if P > R:
            p_n = [xp, yp]
            a += 1
            s_n = s_p
            pasos_aceptados.append([xp, yp, s_p])
            #print "SI"
        else:
            s_n = s_0
            p_n = [x0, y0]
            print "NO"
    elif s_p <= s_0:
        p_n = [xp, yp]
        a += 1
        s_n = s_p
        pasos_aceptados.append([xp, yp, s_p])
        #print "SI"
    return p_n, s_n, a


def monte_carlo(p0, prior_params, N, datos, error=1, d=0.1, modelo=0):
    a = 0
    muestra_met = np.zeros((N, 2))
    pasos_aceptados = []
    chi_cuad = np.zeros(N)
    muestra_met[0] = [p0[0], p0[1]]
    #mu_0 = marginalizar_mu0(p0, datos)
    #dats = datos[0], datos[1] - mu_0, datos[2]
    xi = xi_cuadrado(p0, datos, mu_th)
    chi_cuad[0] = xi[0]
    paso_0 = [p0[0], p0[1], chi_cuad[0]]
    pasos_aceptados = [paso_0]
    #chi_cuad[0] = L_0[1]
    for i in range(1, N):
        P_M = paso_metropolis(muestra_met[i-1], pasos_aceptados, prior_params, datos, chi_cuad[i-1], a, error, d, modelo)
        muestra_met[i] = P_M[0]
        chi_cuad[i] = P_M[1]
        a = P_M[2]
        #posterior_p0 = P_M[4]
        print("contador: "+str(i))
    # guardar datos
    np.save('d.npy', muestra_met)
    np.save('e.npy', chi_cuad)
    np.save('f.npy', a)
    np.save('g.npy', pasos_aceptados)
    return muestra_met, chi_cuad, a


def marginalizar_mu0(r, datos, modelo=0):
    z, m_obs, error_obs = datos
    B = np.sum((m_obs - mu_th(r, z, modelo)) / error_obs ** 2)
    C = np.sum(1 / error_obs ** 2)
    mu_0 = B / C
    return mu_0


def chi_cuadrado(p, dat, f):
    x, y, err = dat
    S = np.sum(((y - f(p, x)) ** 2) / err ** 2)
    return S


def xi_cuadrado(p, dat, f, modelo=0):
    z, m_obs, error_obs = dat
    dif = m_obs - f(p, z, modelo)
    A = np.sum((dif ** 2) / error_obs ** 2)
    B = np.sum(dif / error_obs ** 2)
    C = np.sum(1 / error_obs ** 2)
    S = A - B ** 2 / C
    mu_0 = B / C
    return S, mu_0


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


def mejor_adivinanza(p0, datos):
    z_exp, mu_exp, err = datos
    mu_0 = marginalizar_mu0(p0, datos)
    a = leastsq(residuo_modelo, p0, args=(z_exp, mu_exp))
    return a[0]


# inicializacion
z, mu, mu_err = leer_archivo('SnIa_data.txt')
datos = z, mu, mu_err
adivinanza1 = [0.29, 0.2, 0.78, 0.3]
'''
p00 = 0.78, 0.29
mu_0 = marginalizar_mu0(p00, datos)
dats = datos[0], datos[1] - mu_0, datos[2]
print chi_cuadrado(p00, dats, mu_th)
'''
N = 5000
p0 = np.random.uniform(0, 1), np.random.uniform(0, 1)
t0 = time.time()
resultados = monte_carlo(p0, adivinanza1, N, datos, d=0.05)
tf=time.time()-t0
print("tiempo: "+str(tf))
A = np.load('d.npy')
B = np.load('e.npy')
C = np.load('f.npy')
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
