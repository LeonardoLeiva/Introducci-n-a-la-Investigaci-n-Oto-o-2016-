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

LP = input("letra? (Si = 1/No = 0) = ")
if int(LP) == 1:
     muestra_met = np.load('a.npy')
     chi_cuad = np.load('b.npy')
     a = np.load('c.npy')
     pasos_aceptados = np.load('d.npy')
     tf = np.load('e.npy')
elif int(LP) == 0:
    muestra_met = np.load('densidades.npy')
    chi_cuad = np.load('chi_cuadrado.npy')
    a = np.load('cantidadpasosaceptados.npy')
    pasos_aceptados = np.load('densidadesaceptadas.npy')
    tf = np.load('tiempo.npy')
else:
    print "error"
numerodelacadena = input("numero de la cadena? (ingresar un entero) =")
#guardar
np.savetxt('tiempo'+str(numerodelacadena)+'.dat', np.asarray(tf).reshape(1,))
np.savetxt('densidades'+str(numerodelacadena)+'.dat', muestra_met)
np.savetxt('chi_cuadrado'+str(numerodelacadena)+'.dat', chi_cuad)
np.savetxt('cantidadpasosaceptados'+str(numerodelacadena)+'.dat', np.asarray(a).reshape(1,))
np.savetxt('densidadesaceptadas'+str(numerodelacadena)+'.dat', pasos_aceptados)
