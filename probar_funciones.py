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


def D_L(p, z, modelo=0, paso=20, DL_0=0):
    d_l = Distancia(p, z, modelo) * (1 + z)
    return np.asarray(d_l)


def mu_th(r, z, modelo=0):
    p = r[0], r[1]
    mu = 5 * np.log10(D_L(p, z, modelo, 25))
    return mu


def mu_theo(r, z, modelo=0):
    p = r[0], r[1]
    print p
    graficar_varios(p)
    mu = 5 * np.log10(D_L(p, z, modelo))
    return mu


def xi_cuadrado(p, dat, f, modelo=0):
    z, m_obs, error_obs = dat
    dif = m_obs - f(p, z, modelo)
    A = np.sum((dif ** 2) / error_obs ** 2)
    B = np.sum(dif / error_obs ** 2)
    C = np.sum(1 / error_obs ** 2)
    S = A - B ** 2 / C
    mu_0 = B / C
    return S, mu_0


#graficos
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


def graficar():
    '''
    grafica los datos experimentales
    '''
    z, mu, err = leer_archivo('SnIa_data.txt')
    fig = plt.figure()
    fig.clf()
    ax1 = fig.add_subplot(111)
    ax1.plot(z, mu)
    ax1.set_xlabel("Redshift")
    ax1.set_ylabel("Modulo de la Distancia")
    plt.legend(loc=4)
    plt.savefig("gradatosexp.png")
    plt.draw()
    plt.show()


def graficar_varios(p):
    '''
    solo se agregan barras de errores
    '''
    dat = leer_archivo('SnIa_data.txt')
    z, mu, err = dat
    t0 = time.time()
    mu_nuevo = mu_th(p, z)
    tf=time.time()-t0
    print("tiempo: "+str(tf))
    s, mu_0 = xi_cuadrado(p, dat, mu_th)
    mu_fin = mu_nuevo + mu_0
    print ("chi2: "+str(s))
    fig = plt.figure()
    fig.clf()
    ax1 = fig.add_subplot(111)
    ax1.errorbar(z, mu, err, xerr=None,
         fmt='.', ecolor='red', elinewidth=None, capsize=None,
         barsabove=False, lolims=False, uplims=False,
         xlolims=False, xuplims=False, errorevery=1,
         capthick=None)
    ax1.plot(z, mu_fin, '+')
    ax1.set_xlabel("Redshift")
    ax1.set_ylabel("Modulo de la Distancia")
    plt.legend(loc=4)
    plt.savefig("graexpmasteo.png")
    plt.draw()
    plt.show()


def graficar_teorica(p):
    '''
    grafica el mu teorico segun la aproximacion de chi cuadrado de intento1.py
    '''
    dat = leer_archivo('SnIa_data.txt')
    z, mu, err = dat
    mu_n = mu_th(p, z)
    s, mu_0 = xi_cuadrado(p, dat, mu_th)
    mu_fin = mu_n + mu_0
    print ("chi2: "+str(s))
    fig1 = plt.figure()
    fig1.clf()
    ax1 = fig1.add_subplot(111)
    ax1.plot(z, mu_n, '+')
    ax1.set_xlabel("X")
    ax1.set_ylabel("W(x)")
    plt.savefig("grateor.png")
    plt.draw()
    plt.show()


#optimizar
def residuo_modelo(p, z_exp, mu_exp):
    '''
    diferencia entre los valores del modelo y los experimentales
    '''
    err = mu_exp - mu_th(p, z_exp)
    return err


def optimizar(p0):
    '''
    calcula los mejores parametros que se ajustan al modelo
    '''
    dat = leer_archivo('SnIa_data.txt')
    z, mu, err = dat
    opt = leastsq(residuo_modelo, p0, args=(z, mu))
    return opt


#inicializacion
#dm = input("")
#dde = input("")
p0 = 0.669222206370, 0.916308050753
#op = optimizar(p0)
#p1 = op[0]
#graficar()
graficar_varios(p0)
graficar_teorica(p0)
