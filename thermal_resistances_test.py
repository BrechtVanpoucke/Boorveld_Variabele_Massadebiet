# -*- coding: utf-8 -*-
""" Example of calculation of grout and ground temperatures using the multipole
     method.

     The thermal resistances of a borehole with two pipes are evaluated using
     the multipole method of Claesson and Hellstrom (2011). Based on the
     calculated thermal resistances, the heat flows from the pipes required to
     obtain pipe temperatures of 1 degC are evaluated. The temperatures in and
     around the borehole with 2 pipes are then calculated. Results are verified
     against the results of Claesson and Hellstrom (2011).

"""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import AutoMinorLocator
import pygfunction as gt


def main():
 """
    # -------------------------------------------------------------------------
     # Simulation parameters
    # -------------------------------------------------------------------------

     # Borehole dimensions
     r_b = 0.070         # Borehole radius (m)
     # Pipe dimensions
     n_p = 2             # Number of pipes
     # Pipe outer radius (m)
     rp_out = 0.02*np.ones(n_p)

     # Pipe positions
     # Single U-tube [(x_1, y_1), (x_2, y_2)]
     pos_pipes = [(0.03, 0.00), (-0.03, 0.02)]

     # Ground properties
     k_s = 2.5           # Ground thermal conductivity (W/m.K)

     # Grout properties
     k_g = 1.5           # Grout thermal conductivity (W/m.K)

     # Fluid properties
     # Fluid to outer pipe wall thermal resistance (m.K/W)
     R_fp = 1.2/(2*np.pi*k_g)*np.ones(n_p)

     # Borehole wall temperature (degC)
     T_b = 0.0

     # Fluid temperatures (degC)
     T_f = np.array([1., 1.])

     # Path to validation data
     filePath = './ClaHel11_multipole_temperature.txt'

     # Thermal resistances for J=3
     R_Claesson = 0.01*np.array([25.592, 1.561, 25.311])

     # Number of multipoles per pipe
     J = 3

     # -------------------------------------------------------------------------
     # Evaluate the internal thermal resistances
     # -------------------------------------------------------------------------

     # Thermal resistances
     (R, Rd) = gt.pipes.thermal_resistances(pos_pipes, rp_out, r_b, k_s, k_g,
                                            R_fp, J=3)
     print(50*'-')
     print('Thermal resistance:\t100*R11\t100*R12\t100*R22')
     print('Claesson and Hellstrom:\t{:.3f}\t{:.3f}\t{:.3f}'.format(
             100*R_Claesson[0], 100*R_Claesson[1], 100*R_Claesson[2]))
     print('Present:\t\t{:.3f}\t{:.3f}\t{:.3f}'.format(
             100*R[0,0], 100*R[0,1], 100*R[1,1]))
     print(50*'-')

     # Heat flows
     Q = np.linalg.solve(R, T_f - T_b)

     # -------------------------------------------------------------------------
     # Temperatures along y=0.
     # -------------------------------------------------------------------------

     # Grid points to evaluate temperatures
     x = np.linspace(-0.1, 0.1, num=200)
     y = np.zeros_like(x)

     # Evaluate temperatures using multipole method
     (T_f, T, it, eps_max) = gt.pipes.multipole(pos_pipes, rp_out, r_b, k_s,k_g, R_fp, T_b, Q, J, x_T=x, y_T=y)

     #load validation data
     data = np.loadtxt(filePath, skiprows=1)

     # Configure figure and axes
     fig = gt.utilities._initialize_figure()

     ax1 = fig.add_subplot(111)
     # Axis labels
     ax1.set_xlabel(r'x (m)')
     ax1.set_ylabel(r'$T(x,0)$')
     # Axis limits
     ax1.set_xlim([-0.1, 0.1])
     ax1.set_ylim([-0.2, 1.2])
     # Show grid
     ax1.grid()
     gt.utilities._format_axes(ax1)

     ax1.plot(x, T, label='pygfunction')
     ax1.plot(data[:,0], data[:,1], 'ko',
              label='Claesson and Hellstrom (2011)')
     ax1.legend(loc='upper left')

     # Adjust to plot window
     plt.tight_layout()

     # -------------------------------------------------------------------------
     # Temperatures in -0.1 < x < 0.1, -0.1 < y < 0.1
     # -------------------------------------------------------------------------

     # Grid points to evaluate temperatures
     N_xy = 200
     x = np.linspace(-0.1, 0.1, num=N_xy)
     y = np.linspace(-0.1, 0.1, num=N_xy)
     X, Y = np.meshgrid(x, y)

     # Evaluate temperatures using multipole method
     (T_f, T, it, eps_max) = gt.pipes.multipole(pos_pipes, rp_out, r_b, k_s,
                                               k_g, R_fp, T_b, Q, J,
                                               x_T=X.flatten(),
                                               y_T=Y.flatten())

     # Configure figure and axes
     fig = gt.utilities._initialize_figure()

     ax1 = fig.add_subplot(111)
     # Axis labels
     ax1.set_xlabel('x (m)')
     ax1.set_ylabel('y (m)')
     # Axis limits
     plt.axis([-0.1, 0.1, -0.1, 0.1])
     plt.gca().set_aspect('equal', adjustable='box')
     gt.utilities._format_axes(ax1)

     # Borehole wall outline
     borewall = plt.Circle((0., 0.), radius=r_b,
                           fill=False, linestyle='--', linewidth=2.)
     ax1.add_patch(borewall)
     # Pipe outlines
     for n in range(n_p):
         pipe = plt.Circle(pos_pipes[n], radius=rp_out[n],
                           fill=False, linestyle='-', linewidth=4.)
         ax1.add_patch(pipe)
     # Temperature contours
     CS = ax1.contour(X, Y, T.reshape((N_xy, N_xy)),
                      np.linspace(-0.2, 1.0, num=7))
     plt.clabel(CS, inline=1, fontsize=10)

     # Adjust to plot window
     plt.tight_layout()
     return
"""
#variabelen toewijzen
r_out = 0.018 / 2       # Pipe outer radius (m)
r_in = 0.013 / 2        # Pipe inner radius (m)
mu_f = 0.003971  # Fluid dynamic viscosity (kg/m.s)
rho_f = 1030  # Fluid density (kg/m3)
k_f = 0.455
fluid = gt.media.Fluid('MPG', 20., 10)
cp_f = fluid.cp  # Fluid specific isobaric heat capacity (J/kg.K)
epsilon = 7.0e-6    # Pipe roughness (m)
Re = np.arange(1500, 4600, 10)
m_flow_pipe = Re * r_in * np.pi * mu_f / 2 #massadebiet voor Re van 1500 tot 4600
Pr = cp_f * mu_f / k_f

# Initialize result array fdarcy
fdarcy = np.zeros((1, len(m_flow_pipe)))
for j in range(len(m_flow_pipe)):
    m_flow_borehole = m_flow_pipe[j]
    # Fluid to inner pipe wall thermal resistance (Single U-tube)
    fdarcy1 = gt.pipes.fluid_friction_factor_circular_pipe(m_flow_borehole, r_in, mu_f, rho_f, epsilon, tol=1.0e-6)
    fdarcy[0, j] = fdarcy1

# Initialize result array h_f
h_f = np.zeros((1, len(m_flow_pipe)))
for j in range(len(m_flow_pipe)):
    m_flow_borehole = m_flow_pipe[j]
    # Fluid to inner pipe wall thermal resistance (Single U-tube)
    h_f1 = gt.pipes.convective_heat_transfer_coefficient_circular_pipe(m_flow_borehole, r_in, mu_f, rho_f, k_f, cp_f, epsilon)
    h_f[0, j] = h_f1
"""
# Initialize result array Nu
Pr = cp_f * mu_f / k_f
Nu = gt_pipes__Nusselt_number_turbulent_flow(Re, Pr, fdarcy[0,:])
print(Nu)
"""
def _Nusselt_number_turbulent_flow(Re, Pr, fDarcy):

    """
    An empirical equation developed by Volker Gnielinski (1975)
    [#Gnielinski1975]_ based on experimental data for turbulent flow in pipes.
    Cengel and Ghajar (2015, pg. 497) [#Nusselt-CengelGhajar2015]_ say that the
    Gnielinski equation should be preferred for determining the Nusselt number
    in the transition and turbulent region.
    .. math::
        	\\text{Nu} = \\dfrac{(f/8)(\\text{Re}-1000)\\text{Pr}}
        	{1 + 12.7(f/8)^{0.5} (\\text{Pr}^{2/3}-1)} \\;\\;\\;
        	\\bigg(
            \\begin{array}{c}
                0.5 \leq \\text{Pr} \leq 2000 \\\\
                3 \\times 10^5 <  \\text{Re} < 5 \\times 10^6
            \\end{array}
            \\bigg)
    .. note::
        This equation does not apply to Re < 3000.
    Parameters
    ----------
    Re : float
        Reynolds number.
    Pr : float
        Prandlt Number.
    fDarcy : float
        Darcy friction factor.
    Returns
    -------
    Nu : float
        The Nusselt number
    References
    ------------
    .. [#Gnielinski1975] Gnielinski, V. (1975). Neue Gleichungen für
        den Wärme- und den Stoffübergang in turbulent durchströmten Rohren und
        Kanälen. Forschung im Ingenieurwesen, 41(1), 8–16.
        https://doi.org/10.1007/BF02559682
    .. [#Nusselt-CengelGhajar2015] Çengel, Y.A., & Ghajar, A.J. (2015). Heat
        and mass transfer: fundamentals & applications (Fifth edition.).
        McGraw-Hill.
    """

    # Warn the user if the Reynolds number is out of bounds, but don't break


    Nu = 0.125 * fDarcy * (Re - 1.0e3) * Pr / \
        (1.0 + 12.7 * np.sqrt(0.125*fDarcy) * (Pr**(2.0/3.0) - 1.0))
    return Nu

# plotting the points
plot1 = plt. figure(1)
plt.plot(Re, h_f[0, :])
plot2 = plt. figure(2)
NuArray = [_Nusselt_number_turbulent_flow(4000,Pr,fdarcy[0,i]) for i in range(len(Re))]
plt.plot(Re, NuArray)
# naming the x axis
plt.xlabel('Re')
# naming the y axis
plt.ylabel('Nu')
# giving a title to my graph
plt.title('transition zone')

# function to show the plot
plt.show()
print("h=",h_f[0,:], "f=", fdarcy[0,:])
# Main function
if __name__ == '__main__':
    main()