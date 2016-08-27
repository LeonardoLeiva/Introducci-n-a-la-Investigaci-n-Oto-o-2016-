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
    D = limpiardatos(D)[0]
    fig = plt.figure()
    fig.clf()
    ax1 = fig.add_subplot(111)
    ax1.plot(D[0], D[1], '.')
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
    dat = limpiardatos(datos)
    min_chi = np.amin(datos, axis=0)[2]
    chi = datos[:, 2]
    for i in range(len(chi)):
        if chi[i] == min_chi:
            minimo = datos[i, 0], datos[i, 1], datos[i, 2]
            break
    print "Tiempo que tardo el codigo ="+str(tf)
    print "Minimos"
    print "Densidad de Materia = "+str(minimo[0])
    print "Densidad Energia Oscura = "+str(minimo[1])
    print "Chi Cuadrado = "+str(minimo[2])
    d_m = dat[0]
    d_de = dat[1]
    prom_dm = np.mean(d_m)
    std_dm = np.std(d_m)
    prom_dde = np.mean(d_de)
    std_dde = np.std(d_de)
    print "Promedios"
    print "Densidad de Materia = "+str(prom_dm)+"+-"+str(std_dm)
    print "Densidad Energia Oscura = "+str(prom_dde)+"+-"+str(std_dde)


def regionesdeconfianza(datos, numerodelacadena):
    sigma_1 = 2.30
    sigma_2 = 6.18
    sigma_3 = 11.83
    mini = np.amin(datos, axis=0)
    d_m, d_de, chi = limpiardatos(datos, tol=1)[0]
    min_chi = mini[2]
    N_1 = [0]
    N_2 = [0]
    N_3 = [0]
    N_4 = [0]
    for i in range(len(chi)):
        if chi[i] == min_chi:
            minimo = d_m[i], d_de[i], chi[i]
    for j in range(len(chi)):
        if chi[j] <= min_chi + sigma_3:
            if chi[j] <= min_chi + sigma_2:
                if chi[j] <= min_chi + sigma_1:
                    N_1.append(j)
                else:
                    N_2.append(j)
            else:
                N_3.append(j)
        else:
            N_4.append(j)
    X1 = [minimo[0]]
    X2 = [minimo[1]]
    Y1 = [minimo[0]]
    Y2 = [minimo[1]]
    Z1 = [minimo[0]]
    Z2 = [minimo[1]]
    V1 = [minimo[0]]
    V2 = [minimo[1]]
    for k in N_1[1:]:
        X1.append(d_m[k])
        X2.append(d_de[k])
    for l in N_2[1:]:
        Y1.append(d_m[l])
        Y2.append(d_de[l])
    for m in N_3[1:]:
        Z1.append(d_m[m])
        Z2.append(d_de[m])
    for o in N_4[1:]:
        V1.append(d_m[o])
        V2.append(d_de[o])
    fig = plt.figure()
    fig.clf()
    ax1 = fig.add_subplot(111)
    ax1.plot(X1, X2, '.', label='68.3$\%$ de los datos')
    ax1.plot(Y1, Y2, '.', label='95.45$\%$ de los datos')
    ax1.plot(Z1, Z2, '.', label='99.73$\%$ de los datos')
    ax1.plot(V1, V2, '.')
    ax1.set_xlabel("Densidad de Materia")
    ax1.set_ylabel("Densidad de Energia Oscura")
    plt.legend(loc=4)
    plt.savefig("regionesdeconfianza"+str(numerodelacadena)+".png")
    plt.draw()
    plt.show()





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
    return [c, d, e], n


#inicializacion
numerodelacadena = input("numero de la cadena a graficar? (ingresar un entero n >= 1) =")
tf = np.loadtxt('tiempo'+str(numerodelacadena)+'.dat')
densidades = np.loadtxt('densidades'+str(numerodelacadena)+'.dat')
chicuad = np.loadtxt('chi_cuadrado'+str(numerodelacadena)+'.dat')
cantidadaceptados = np.loadtxt('cantidadpasosaceptados'+str(numerodelacadena)+'.dat')
datos = np.loadtxt('densidadesaceptadas'+str(numerodelacadena)+'.dat')
regionesdeconfianza(datos, numerodelacadena)
#graficar_densidades(densidades, numerodelacadena)
#graficar_densidadesburnin(datos, numerodelacadena)
graficar_chicuad(datos, chicuad, numerodelacadena)
datosimportantes(datos, tf)
