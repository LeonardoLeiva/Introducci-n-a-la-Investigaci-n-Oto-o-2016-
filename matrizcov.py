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
np.random.seed(888888)


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
        op = np.absolute((1 + DL ** 2 * (1 - d_m - d_de)))
        E[i] = np.sqrt(op) / hub(p, z[i], modelo)
    return E


def res_EDO(p, DL_0, z_f, z_0=0., paso=100, modelo=0):
    init = DL_0
    z = np.linspace(z_0, z_f, paso)
    sol = odeint(EDO, init, z, args=(p, modelo, ), mxstep=20)
    a = sol[paso - 1, 0]
    return a


def edo_varios(p, z, paso=50, modelo=0):
    '''
    requiere que z este ordenado
    '''
    m = len(z)
    mu = np.zeros(m)
    mu[0] = res_EDO(p, 0., z[0], 0., paso)
    for i in range(1, m):
        if z[i-1] == z[i]:
            mu[i] = mu[i-1]
        elif z[i] >= z[i-1]:
            mu[i] = res_EDO(p, mu[i-1], z[i], z[i-1], paso)
        else:
            print "error"
    return mu


def Distancia(p, z, modelo=0, paso=100, DL_0=0):
    if z.size == 1:
        d = res_EDO(p, DL_0, z, 0, 100, modelo)
    else:
        d = edo_varios(p, z, paso, modelo)
    return d


def D_L(p, z, modelo=0, paso=100, DL_0=0):
    d_l = Distancia(p, z, modelo) * (1 + z)
    return np.asarray(d_l)


def mu_th(r, z, modelo=0):
    p = r[0], r[1]
    mu = 5 * np.log10(D_L(p, z, modelo, 25))
    return mu


def leer_archivo(nombre):
    '''
    lee el archivo
    nombre debe ser un str
    '''
    datos = np.loadtxt(nombre, usecols=(1, 2, 3))
    #dtype = [('z', float), ('mu', float), ('error_mu', float)]
    z = datos[:, 0]
    mu = datos[:, 1]
    err_mu = datos[:, 2]
    n = np.argsort(z)
    k = len(z)
    zf = np.zeros(k)
    muf = np.zeros(k)
    err_muf = np.zeros(k)
    for i in range(k):
        m = n[i]
        zf[i] = z[m]
        muf[i] = mu[m]
        err_muf[i] = err_mu[m]
    return zf, muf, err_muf


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


def prior2(beta, p, model=0):
    '''
    probabilidad a priori que se le da a la adivinanza inicial. dos modelos
    '''
    '''
    mean, std, cov = p
    a0, a1 = mean[0], mean[1]
    b0, b1 = std[0], std[1]
    d_m, d_de = beta
    print d_m, a0, b0
    s0 = ((d_m - a0) / b0) ** 2
    s1 = ((d_de - a1) / b1) ** 2
    s = (s0 + s1) / 2.
    P = np.exp(-s) / (2 * np.pi * b0 * b1)
    '''
    P=1
    return P


def paso_metropolis(p0, pasos_aceptados, prior_params, datos, s_0, a, error=1, d=0.1, modelo=0):
    x0, y0 = p0
    rx = np.random.uniform(low=-1, high=1)
    ry = np.random.uniform(low=-1, high=1)
    xp = x0 + d * rx
    yp = y0 + d * ry
    print [xp, yp]
    while xp < 0 or xp > 1:
        rx = np.random.uniform(low=-1, high=1)
        xp = x0 + d * rx
    while yp < 0 or yp > 1:
        ry = np.random.uniform(low=-1, high=1)
        yp = y0 + d * ry
    xi = xi_cuadrado([xp, yp], datos, mu_th)
    s_p = xi[0]
    if s_p > s_0:
        s = s_0 / s_p
        print s_0 - s_p
        S = np.exp((s_0 - s_p) / 2)
        #S = np.exp(- s / 2)
        pp = prior([xp, yp], prior_params)
        p0 = prior([x0, y0], prior_params)
        L =  prior([xp, yp], prior_params) / prior([x0, y0], prior_params)
        P = S * L
        print S, L
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


def paso_metrop(p0, pasos_aceptados, prior_params, datos, s_0, a, error=1, d=0.1, modelo=0):
    x0, y0 = p0
    mean, std, cov = prior_params
    r = np.random.multivariate_normal(mean, cov)
    xp = x0 + d * r[0]
    yp = y0 + d * r[1]
    print xp, yp
    while xp < 0 or xp > 1:
        r = np.random.multivariate_normal(mean, cov)
        xp = x0 + d * r[0]
        yp = y0 + d * r[1]
    while yp < 0 or yp > 1:
        r = np.random.multivariate_normal(mean, cov)
        xp = x0 + d * r[0]
        yp = y0 + d * r[1]
    xi = xi_cuadrado([xp, yp], datos, mu_th)
    s_p = xi[0]
    if s_p > s_0:
        S = np.exp((s_0 - s_p) / 2)
        print s_0 - s_p
        pp = prior2([xp, yp], prior_params)
        p0 = prior2([x0, y0], prior_params)
        L =  pp / p0
        P = S * L
        print S, L
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


def monte_carlo(p0, prior_params, N, datos, precadena=0, error=1, d=0.1, modelo=0):
    a = 0
    muestra_met = np.zeros((N, 2))
    pasos_aceptados = []
    chi_cuad = np.zeros(N)
    muestra_met[0] = [p0[0], p0[1]]
    xi = xi_cuadrado(p0, datos, mu_th)
    chi_cuad[0] = xi[0]
    paso_0 = [p0[0], p0[1], chi_cuad[0]]
    pasos_aceptados = [paso_0]
    if precadena == 0:
        for i in range(1, N):
            P_M = paso_metropolis(muestra_met[i-1], pasos_aceptados, prior_params, datos, chi_cuad[i-1], a, error, d, modelo)
            muestra_met[i] = P_M[0]
            chi_cuad[i] = P_M[1]
            a = P_M[2]
            print("contador: "+str(i))
    else:
        beta = np.load('densidadesaceptadas.npy')
        a = beta[100:, 0]
        b = beta[100:, 1]
        A = np.asarray([a, b])
        mean = np.mean(A, axis=1)
        std = np.std(A, axis=1)
        cov = np.cov(A)
        prior_params = [mean, std, cov]
        for i in range(1, N):
            P_M = paso_metrop(muestra_met[i-1], pasos_aceptados, prior_params, datos, chi_cuad[i-1], a, error, d, modelo)
            muestra_met[i] = P_M[0]
            chi_cuad[i] = P_M[1]
            a = P_M[2]
            print("contador: "+str(i))
    # guardar datos
    np.save('a.npy', muestra_met)
    np.save('b', chi_cuad)
    np.save('c.npy', a)
    np.save('d.npy', pasos_aceptados)
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


# inicializacion
z, mu, mu_err = leer_archivo('SnIa_data.txt')
datos = z, mu, mu_err
adivinanza1 = [0.29, 0.2, 0.78, 0.3]
N = 100
p0 = np.random.uniform(0, 1), np.random.uniform(0, 1)
print p0
t0 = time.time()
resultados = monte_carlo(p0, adivinanza1, N, datos, precadena=0, d=0.05)
tf=time.time()-t0
print("tiempo: "+str(tf))
np.save('tiempo.npy', tf)
# grafico
A = np.load('d.npy')
B = np.load('b.npy')
C = np.load('c.npy')
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
