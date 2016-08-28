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


def hub(p, z, modelo=0):
    if modelo == 0:
        d_m = p[0]
        d_de = p[1]
        H = np.sqrt(np.absolute(d_m * (1 + z) ** 3 + d_de + (1 - d_m - d_de) * (1 + z) ** 2))
    elif modelo == 1:
        d_m, w_0, w_1 = p
        H = np.sqrt(np.absolute(d_m * (1 + z) ** 3 + (1 - d_m) * (1 + z) ** (3 * (1 + w_0 - w_1)) * np.exp(3 * w_1 * z)))
    return H


def EDO(z, DL, p, modelo=0):
    d_m = p[0]
    d_de = p[1]
    E = np.zeros(len(z))
    for i in range(len(z)):
        if modelo == 0:
            op = np.absolute((1 + DL ** 2 * (1 - d_m - d_de)))
        elif modelo == 1:
            op = np.absolute((1 + DL ** 2 * (1 - d_m)))
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
    if modelo == 0:
        p = r[0], r[1]
    elif modelo == 1:
        p = r[0], r[1], r[2]
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
    if modelo == 0:
        d_m, d_de = beta
        a0, b0, a1, b1 = p
        s0 = ((d_m - a0) / b0) ** 2
        s1 = ((d_de - a1) / b1) ** 2
        s = (s0 + s1) / 2.
        P = np.exp(-s) / (2 * np.pi * b0 * b1)
        return P
    elif model == 1:
        '''
        d_m, d_de, w_0, w_a = beta
        a0, a1, a2, a3, b0, b1, b2, b3 = p
        s0 = ((d_m - a0) / b0) ** 2
        s1 = ((d_de - a1) / b1) ** 2
        s2 = ((w_0 - a2) / b2) ** 2
        s3 = ((w_a - a3) / b3) ** 2
        s = (s0 + s1 + s2 + s3) / 2.
        P = np.exp(-s) / (4 * np.pi ** 2 * b0 * b1 * b2 * b3)
        '''
        P = 1
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
    if modelo == 0:
        x0, y0 = p0
        rx = np.random.uniform(low=-1, high=1)
        ry = np.random.uniform(low=-1, high=1)
        xp = x0 + d * rx
        yp = y0 + d * ry
        while xp < 0 or xp > 1:
            rx = np.random.uniform(low=-1, high=1)
            xp = x0 + d * rx
        while yp < 0 or yp > 1:
            ry = np.random.uniform(low=-1, high=1)
            yp = y0 + d * ry
        Xp = xp, yp
        X0 = x0, y0
    elif modelo == 1:
        x0, y0, z0 = p0
        rx = np.random.uniform(low=-1, high=1)
        ry = np.random.uniform(low=-1, high=1)
        rz = np.random.uniform(low=-1, high=1)
        xp = x0 + d * rx
        yp = y0 + d * ry
        zp = z0 + d * rz
        while xp < 0 or xp > 1:
            rx = np.random.uniform(low=-1, high=1)
            xp = x0 + d * rx
        Xp = xp, yp, zp
        X0 = x0, y0, z0
    xi = xi_cuadrado(Xp, datos, mu_th)
    s_p = xi[0]
    if s_p > s_0:
        s = s_0 / s_p
        S = np.exp((s_0 - s_p) / 2)
        #S = np.exp(- s / 2)
        pp = prior(Xp, prior_params)
        p0 = prior(X0, prior_params)
        L =  1
        P = S * L
        R = np.random.uniform(0, 1)
        if P > R:
            p_n  = Xp
            a += 1
            s_n = s_p
            if modelo == 0:
                pasos_aceptados.append([xp, yp, s_p])
            elif modelo == 1:
                pasos_aceptados.append([xp, yp, zp, s_p])
            #print "SI"
        else:
            s_n = s_0
            p_n = X0
            print "NO"
    elif s_p <= s_0:
        p_n = Xp
        a += 1
        s_n = s_p
        if modelo == 0:
            pasos_aceptados.append([xp, yp, s_p])
        elif modelo == 1:
            pasos_aceptados.append([xp, yp, zp, s_p])
        #print "SI"
    return p_n, s_n, a


def paso_metrop(p0, pasos_aceptados, prior_params, datos, s_0, a, error=1, d=0.1, modelo=0):
    mean, std, cov = prior_params
    R = np.random.multivariate_normal(mean, cov)
    r = R - mean
    if modelo == 0:
        x0, y0 = p0
        xp = x0 + d * r[0]
        yp = y0 + d * r[1]
        while xp < 0 or xp > 1:
            R = np.random.multivariate_normal(mean, cov)
            r = R - mean
            xp = x0 + d * r[0]
            yp = y0 + d * r[1]
        while yp < 0 or yp > 1:
            R = np.random.multivariate_normal(mean, cov)
            r = R - mean
            xp = x0 + d * r[0]
            yp = y0 + d * r[1]
        Xp = xp, yp
    elif modelo == 1:
        x0, y0, z0 = p0
        xp = x0 + d * r[0]
        yp = y0 + d * r[1]
        zp = z0 + d * r[1]
        while xp < 0 or xp > 1:
            R = np.random.multivariate_normal(mean, cov)
            r = R - mean
            xp = x0 + d * r[0]
            yp = y0 + d * r[1]
            zp = z0 + d * r[1]
        Xp = xp, yp, zp
    xi = xi_cuadrado(Xp, datos, mu_th)
    s_p = xi[0]
    if s_p > s_0:
        S = np.exp((s_0 - s_p) / 2)
        pp = prior2(Xp, prior_params)
        p0 = prior2(X0, prior_params)
        L =  pp / p0
        P = S * L
        R = np.random.uniform(0, 1)
        if P > R:
            p_n = Xp
            a += 1
            s_n = s_p
            if modelo == 0:
                pasos_aceptados.append([xp, yp, s_p])
            elif modelo == 1:
                pasos_aceptados.append([xp, yp, zp, s_p])
            #print "SI"
        else:
            s_n = s_0
            p_n = X0
            print "NO"
    elif s_p <= s_0:
        p_n = Xp
        a += 1
        s_n = s_p
        if modelo == 0:
            pasos_aceptados.append([xp, yp, s_p])
        elif modelo == 1:
            pasos_aceptados.append([xp, yp, zp, s_p])
        #print "SI"
    return p_n, s_n, a


def monte_carlo(p0, prior_params, N, datos, error=1, d=0.1, modelo=0):
    mantener = input("desea conservar el paso"+str(d)+"? (Y=1/N=0) =")
    if mantener == 0:
        paso = input("nuevo paso? =")
        d = float(paso)
    primcadena = input("primera cadena? (Y=1/N=0) = ")
    numerodelacadena = input("numero de la cadena nueva? (ingresar un entero n >= 1) =")
    t0 = time.time()
    a = 0
    chi_cuad = np.zeros(N)
    xi = xi_cuadrado(p0, datos, mu_th, modelo)
    if modelo == 0:
        muestra_met = np.zeros((N, 2))
        muestra_met[0] = [p0[0], p0[1]]
        paso_0 = [p0[0], p0[1], chi_cuad[0]]
    elif modelo == 1:
        muestra_met = np.zeros((N, 3))
        muestra_met[0] = [p0[0], p0[1], p0[2]]
        paso_0 = [p0[0], p0[1], p0[2], chi_cuad[0]]
    pasos_aceptados = []
    chi_cuad[0] = xi[0]
    pasos_aceptados = [paso_0]
    if primcadena == 1:
        for i in range(1, N):
            P_M = paso_metropolis(muestra_met[i-1], pasos_aceptados, prior_params, datos, chi_cuad[i-1], a, error, d, modelo)
            muestra_met[i] = P_M[0]
            chi_cuad[i] = P_M[1]
            a = P_M[2]
            print("contador: "+str(i))
    elif primcadena == 0:
        beta = np.loadtxt('densidadesaceptadas'+str(numerodelacadena - 1)+'.dat')
        limp = limpiardatos(beta)
        if modelo == 0:
            am, bm, chim = limp
            A = np.asarray([am, bm])
        elif modelo == 1:
            am, bm, cm, chim = limp
            A = np.asarray([am, bm, cm])
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
    else:
        print "error del numero de la cadena"
    # guardar datos
    tf=time.time()-t0
    print("tiempo: "+str(tf))
    np.savetxt('tiempo'+str(numerodelacadena)+'.dat', np.asarray(tf).reshape(1,))
    np.savetxt('densidades'+str(numerodelacadena)+'.dat', muestra_met)
    np.savetxt('chi_cuadrado'+str(numerodelacadena)+'.dat', chi_cuad)
    np.savetxt('cantidadpasosaceptados'+str(numerodelacadena)+'.dat', np.asarray(a).reshape(1,))
    np.savetxt('densidadesaceptadas'+str(numerodelacadena)+'.dat', pasos_aceptados)
    return muestra_met, chi_cuad, a


def limpiardatos(datos, modelo=0, tol=3):
    if modelo == 0:
        a = datos[:, 0]
        b = datos[:, 1]
        chi = datos[:, 2]
        mean = np.mean(chi)
        std = np.std(chi)
        criterio = mean + tol * std
        n = len(chi)
        for i in range(n):
            if criterio >= chi[i]:
                break
        c = a[i:]
        d = b[i:]
        e = chi[i:]
        res = c, d, e
    elif modelo == 1:
        a = datos[:, 0]
        b = datos[:, 1]
        c = datos[:, 2]
        chi = datos[:, 3]
        mean = np.mean(chi)
        std = np.std(chi)
        criterio = mean + tol * std
        n = len(chi)
        for i in range(n):
            if criterio >= chi[i]:
                break
        d = a[i:]
        e = b[i:]
        f = c[i:]
        g = chi[i:]
        res = d, e, f, g
    return res


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
N = input("largo de la cadena? (numero entero) =")
N = int(N)
modelo = input("modelo? (numero entero) =")
if modelo == 0:
    p0 = np.random.uniform(0, 1), np.random.uniform(0, 1)
elif modelo == 1:
    p0 = np.random.uniform(0, 1), np.random.uniform(-1.5, 0.5), np.random.uniform(-2.0, 1.5)
resultados = monte_carlo(p0, adivinanza1, N, datos, d=0.05, modelo=modelo)
