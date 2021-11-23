# Imports
import pygfunction as gt
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from numpy import vectorize
from matplotlib.ticker import AutoMinorLocator

"""
# define variables
Re = 3500  # Reynoldsgetal
d = 0.013  # Pipe's diameter
mu = 0.004  # Monopropyleen  dynamische viscosity Pa s @25 C
rho = 1030  # Monopropyleen density kg/m^3
L = 300  # Pipe's lengte
e = 0.000007  # Pipe's ruwheid
n = 120 # aantal boorputten
nu = 0.8 # efficientie motor

# This function calculates the velocity of the fluid
def vel(Re, d):
    v = Re*mu / (rho*d)
    return v

#this function calculates the flow rate
def flr(v, d, n):
    Q = v * d**2 / 4 * np.pi*n
    return Q


# This function calculates the friction coefficient Lambda
def lbd(Re):

    # If RE < 2300 flow is laminar, otherwise it is turbulent
    if Re <= 2300:
        lbda = 64 / Re
    else:
        lbda = 0.316/Re**0.25
    return lbda

## Vectorize function for calculating lambda.
#lbdv = vectorize(lbd)

# This function calculates the pressure drop
def pressureDrop(Re, d, L=L, rho=rho):
    pDrop = 0.5 *  L / d * rho * lbda * v**2
    return pDrop

def PowerPump(pDrop, Q,  nu=nu):
    P = pDrop*Q/nu
    return P

# Calculate pressure drop
v=vel(Re, d)

# Calculate flow rate
Q=flr(v, d, n)

# Calculate f
lbda = lbd(Re)

# Calculate pressure drop
pdrop = pressureDrop(Re, d, L, rho)


#Calculate power
P=PowerPump(pdrop, Q,  nu)
"""

def main():
    # -------------------------------------------------------------------------
    # Simulation parameters
    # -------------------------------------------------------------------------

    # Number of boreholes
    nBoreholes = 5
    # Borehole dimensions
    D = 4.0             # Borehole buried depth (m)
    # Borehole length (m)
    H = 150.
    r_b = 0.075         # Borehole radius (m)
    B = 10             # Borehole spacing (m)

    # Pipe dimensions
    r_out = 0.018 / 2       # Pipe outer radius (m)
    r_in = 0.013 / 2        # Pipe inner radius (m)
    D_s = 0.055          # Shank spacing (m) from centre of borehole
    epsilon = 7.0e-6    # Pipe roughness (m)

    # Pipe positions
    # Single U-tube [(x_in, y_in), (x_out, y_out)]
    pos_pipes = [(-D_s, 0.), (D_s, 0.)]

    # Ground properties
    k_s = 3.5           # Ground thermal conductivity (W/m.K)

    # Grout properties
    k_g = 1.0           # Grout thermal conductivity (W/m.K)

    # Pipe properties
    k_p = 0.4           # Pipe thermal conductivity (W/m.K)

    # Fluid properties
    mu_f = 0.003971  # Fluid dynamic viscosity (kg/m.s)
    rho_f = 1030  # Fluid density (kg/m3)
    # Re=[500, 1000, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000, 3500, 4000, 4500, 5000, 5500, 6000]
    # Total fluid mass flow rate per borehole (kg/s), from 0.01 kg/s to 1 kg/s FOR Re it from 1500 to 3000
    m_flow_network =np.arange(1500, 3100, 100) * r_in * np.pi * mu_f / 2
    # The fluid is propylene-glycol (20 %) at 20 degC
    fluid = gt.media.Fluid('MPG', 20., 10)
    cp_f = fluid.cp     # Fluid specific isobaric heat capacity (J/kg.K)
    k_f = 0.455       # Fluid thermal conductivity (W/m.K)

    # -------------------------------------------------------------------------
    # Borehole field
    # -------------------------------------------------------------------------

    boreField = []
    bore_connectivity = []
    for i in range(nBoreholes):
        x = i*B
        borehole = gt.boreholes.Borehole(H, D, r_b, x, 0.)
        boreField.append(borehole)
        # Boreholes are connected in series: The index of the upstream
        # borehole is that of the previous borehole
        bore_connectivity.append(i - 1)

    # -------------------------------------------------------------------------
    # Evaluate the effective bore field thermal resistance
    # -------------------------------------------------------------------------

    # Initialize result array
    R = np.zeros((nBoreholes, len(m_flow_network)))
    for i in range(nBoreholes):
        for j in range(len(m_flow_network)):
            nBoreholes = i + 1
            # Boreholes are connected in series
            m_flow_borehole = m_flow_network[j]
            # Boreholes are single U-tube
            m_flow_pipe = m_flow_borehole

            # Pipe thermal resistance
            R_p = gt.pipes.conduction_thermal_resistance_circular_pipe(
                    r_in, r_out, k_p)
            # Fluid to inner pipe wall thermal resistance (Single U-tube)
            h_f = gt.pipes.convective_heat_transfer_coefficient_circular_pipe(
                    m_flow_pipe, r_in, mu_f, rho_f, k_f, cp_f, epsilon)
            R_f = 1.0/(h_f*2*np.pi*r_in)

            # Single U-tube, same for all boreholes in the bore field
            UTubes = []
            for borehole in boreField:
                SingleUTube = gt.pipes.SingleUTube(
                    pos_pipes, r_in, r_out, borehole, k_s, k_g, R_f + R_p)
                UTubes.append(SingleUTube)
            network = gt.networks.Network(
                boreField[:nBoreholes],
                UTubes[:nBoreholes],
                bore_connectivity=bore_connectivity[:nBoreholes])

            # Effective bore field thermal resistance
            R_field = gt.networks.network_thermal_resistance(
                network, m_flow_network[j], cp_f)
            # Add to result array
            R[i,j] = R_field

    # -------------------------------------------------------------------------
    # Plot bore field thermal resistances
    # -------------------------------------------------------------------------

    # Configure figure and axes
    fig = gt.utilities._initialize_figure()

    ax1 = fig.add_subplot(111)
    # Axis labels
    ax1.set_xlabel(r'$\dot{m}$ [kg/s]')
    ax1.set_ylabel(r'$R^*_{field}$ [m.K/W]')
    # Axis limits
    ax1.set_xlim([0., 0.2])
    ax1.set_ylim([0., 2.])

    gt.utilities._format_axes(ax1)

    # Bore field thermal resistances
    ax1.plot(m_flow_network, R[0,:], '-', label='1 borehole')
    ax1.plot(m_flow_network, R[2,:], '--', label='3 boreholes')
    ax1.plot(m_flow_network, R[4,:], '-.', label='5 boreholes')
    ax1.legend()
    # Adjust to plot window
    plt.tight_layout()
    plt.show()
    print(R[0,:])

#pos = [(-0.06, 0.),(0.,0.), (0.06, 0.)]
#R, Rd = gt.pipes.thermal_resistances(pos, 0.01, 0.075, 2., 1., 0.1, J=2)
main()



# print('v=',v,'Q=',Q,'lambda=',lbda,'pressure drop =',pdrop,'P=',P,'R=',R,'Rd=',Rd)

