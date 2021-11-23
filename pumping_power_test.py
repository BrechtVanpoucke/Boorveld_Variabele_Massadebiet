# Imports
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import AutoLocator, LinearLocator, FormatStrFormatter
from numpy import vectorize

# Let's define the variables
Re = np.linspace(500, 4000, 50)  # Flow rate m^3/s
d = np.linspace(0.013, 0.04, 50)  # Pipe's diameter m
mu = 0.004  # Water viscosity Pa s @25 C
rho = 1030  # Water density kg/m^3
L = 300  # Pipe's length m
N = 120
Nu = 0.8


# This function calculates the velocity of the fluid
def vel(Re, d, mu=mu, rho=rho):
    v = (Re * mu) / (d * rho)
    return v


# This function calculates the friction coefficient Lambda
def lbd(Re):
    # If RE < 2300 flow is laminar, otherwise it is turbolent
    if Re <= 2300:
        lbda = 64 / Re
    else:
        lbda = 0.365 / (Re ** 0.25)
    return lbda


# This function calculates the pressure drop
def pressureDrop(Re, d, L=L):
    v = vel(Re, d)
    pDrop = 0.5 * rho * lbdv(Re) * v ** 2 * L / d
    return pDrop

def flowRate(Re, d, N=N):
    v = vel(Re, d)
    Q = v * N * d ** 2 * np.pi / 4
    return Q

def power(pDrop, Re, d, N=N, Nu=Nu):
    Q=flowRate(Re, d, N)
    P = pDrop * Q / Nu
    return P

# This function is used to plot the contour plot
def contour(X, Y, Z, N=20):
    plt.figure()
    CS = plt.contour(X, Y, Z, N)
    plt.clabel(CS, inline=1, fontsize=10)
    plt.title('contourplot P in functie van Re en d')
    plt.xlabel('Re')
    plt.ylabel('d (m)')
    plt.show()


# This function is used to make a 3d plot
def dplot(X, Y, Z):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(Re, d, P, rstride=1, cstride=1, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)
    ax.zaxis.set_major_locator(AutoLocator())
    ax.xaxis.set_major_locator(AutoLocator())
    ax.yaxis.set_major_locator(LinearLocator(6))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%2.f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%1.3f'))
    ax.set_xlabel('Re ',labelpad=10)
    ax.set_ylabel('d (m)',labelpad=10)
    ax.set_zlabel('P (W)',labelpad=10)
    ax.set_title('Vermogen pomp in functie van buisdiameter en Reynoldsgetal')
    #fig.colorbar(surf, shrink=0.5, aspect=10)
    plt.show()


# Vectorize function for calculating lambda. Sooo useful vectorize!!!
lbdv = vectorize(lbd)

# Grid
Re, d = np.meshgrid(Re, d)

# Calculate pressure drop
v = vel(Re, d)
Q = flowRate(v, d, N)
pDrop=pressureDrop(Re,d,L)
P = power(pDrop, Re, d, N, Nu)
#print(P)
#print(Re)

# Plot
dplot(Re, d, P)

# Contour plot
contour(Re, d, P)