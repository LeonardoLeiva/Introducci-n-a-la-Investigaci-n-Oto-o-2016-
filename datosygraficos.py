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
from mpl_toolkits.mplot3d import Axes3D


def graficar_densidades(D, n):
    fig = plt.figure()
    fig.clf()
    ax1 = fig.add_subplot(111)
    ax1.plot(D[:, 0], D[:, 1], '.')
    ax1.set_xlabel("Densidad de Materia")
    ax1.set_ylabel("Densidad de Energia Oscura")
    plt.legend(loc=4)
    plt.savefig("mcmctodasdensidades"+str(n)+".png")
    plt.draw()
    plt.show()


def graficar_densidadesburnin(D, n, tol=5):
    D = limpiardatos(D)
    fig = plt.figure()
    fig.clf()
    ax1 = fig.add_subplot(111)
    ax1.plot(D[:, 0], D[:, 1], '.')
    ax1.set_xlabel("Densidad de Materia")
    ax1.set_ylabel("Densidad de Energia Oscura")
    plt.legend(loc=4)
    plt.savefig("mcmcdensidades"+str(n)+".png")
    plt.draw()
    plt.show()


def graficar_chicuad(D, C, n):
    m = len(C)
    N = np.linspace(1, m, num=m)
    fig = plt.figure()
    fig.clf()
    ax2 = fig.add_subplot(111)
    ax2.plot(N, C, '.')
    ax2.set_xlabel("Paso de la Cadena")
    ax2.set_ylabel("Valor de Chi cuadrado")
    plt.legend(loc=4)
    plt.savefig("chicuad"+str(n)+".png")
    plt.draw()
    plt.show()
    fig=plt.figure(1)
    fig.clf()
    ax=fig.add_subplot(111, projection='3d')
    ax.set_aspect('equal')
    ax.plot(D[:, 0], D[:, 1], D[:, 2], '.')
    ax.set_xlabel('Densidad de Materia')
    ax.set_ylabel('Densidad de Energia Oscura')
    ax.set_zlabel('Chi Cuadrado')
    plt.savefig("chicuadenelespacio"+str(n)+".png")
    plt.show()


def datosimportantes(datos, tf):
    minimo = np.amin(datos, axis=0)
    print "Tiempo que tardo el codigo ="+str(tf)
    print "Minimos"
    print "Densidad de Materia = "+str(minimo[0])
    print "Densidad Energia Oscura = "+str(minimo[1])
    print "Chi Cuadrado = "+str(minimo[2])
    d_m = datos[:, 0]
    d_de = datos[:, 1]
    prom_dm = np.mean(d_m)
    std_dm = np.std(d_m)
    prom_dde = np.mean(d_de)
    std_dde = np.std(d_de)
    print "Promedios"
    print "Densidad de Materia = "+str(prom_dm)+"+-"+str(std_dm)
    print "Densidad Energia Oscura = "+str(prom_dde)+"+-"+str(std_dde)


def limpiardatos(datos, tol=3):
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
    return c, d, e


#inicializacion
numerodelacadena = input("numero de la cadena a graficar? (ingresar un entero n >= 1) =")
tf = np.loadtxt('tiempo'+str(numerodelacadena)+'.dat')
densidades = np.loadtxt('densidades'+str(numerodelacadena)+'.dat')
chicuad = np.loadtxt('chi_cuadrado'+str(numerodelacadena)+'.dat')
cantidadaceptados = np.loadtxt('cantidadpasosaceptados'+str(numerodelacadena)+'.dat')
datos = np.loadtxt('densidadesaceptadas'+str(numerodelacadena)+'.dat')
graficar_densidades(densidades, numerodelacadena)
graficar_chicuad(datos, chicuad, numerodelacadena)
datosimportantes(datos, tf)
