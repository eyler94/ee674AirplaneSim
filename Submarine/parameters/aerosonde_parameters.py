import sys
sys.path.append('..')
import numpy as np
from tools.tools import Euler2Quaternion

######################################################################################
                #   Initial Conditions
######################################################################################
#   Initial conditions for MAV
pn0 = 0.  # initial north position
pe0 = 0.  # initial east position
pd0 = 100.0  # initial down position
u0 = 1.  # initial velocity along body x-axis
v0 = 0.  # initial velocity along body y-axis
w0 = 0.  # initial velocity along body z-axis
phi0 = 0.  # initial roll angle
theta0 =  0.  # initial pitch angle
psi0 = 0.0  # initial yaw angle
p0 = 0  # initial roll rate
q0 = 0  # initial pitch rate
r0 = 0  # initial yaw rate
Va0 = np.sqrt(u0**2+v0**2+w0**2)
#   Quaternion State
e = Euler2Quaternion(phi0, theta0, psi0)
e0 = e.item(0)
e1 = e.item(1)
e2 = e.item(2)
e3 = e.item(3)


######################################################################################
                #   Physical Parameters
######################################################################################
mass = 18.826 #kg
Ix = 0.0727 #kg m^2
Iy = 1.77
Iz = 1.77
Ixz = 0.
gravity = 9.8
rho = 997
length = 1.391
radius = 0.076
xfin = 0.537
cg = np.array([-0.012,0.,0.0048])
xg = 0.#-0.012
yg = 0.
zg = 0.#0048

######################################################################################
                #   Sub Parameters
######################################################################################
Xudot = 4.21e-1
Xpdot =
Yvdot =
Zwdot =
Yvdot = Zwdot
Mwdot =
Nvdot = -Mwdot
Yrdot = Nvdot
Zqdot = Mwdot
Mqdot =
Nrdot = Mqdot
Xwq = Zwdot
Xvr = -Yvdot
Yura = Xudot
Ypq = -Zqdot
Zuqa = -Xudot
Zrp = Yrdot
Muwa = -(Zwdot-Xudot)
Mrp = (Xpdot-Nrdot)
Nuva = -(Xudot-Yvdot)




# # S_wing = 0.55
# # b = 2.8956
# # c = 0.18994
# # S_prop = 0.2027
# # e = 0.9
# # AR = (b**2) / S_wing
#
# ######################################################################################
#                 #   Longitudinal Coefficients
# ######################################################################################
# C_L_0 = 0.23
# C_D_0 = 0.043
# C_m_0 = 0.0135
# C_L_alpha = 5.61
# C_D_alpha = 0.03
# C_m_alpha = -2.74
# C_L_q = 7.95
# C_D_q = 0.0
# C_m_q = -38.21
# C_L_delta_e = 0.13
# C_D_delta_e = 0.0135
# C_m_delta_e = -0.99
# M = 50.0
# alpha0 = 0.47
# epsilon = 0.16
# C_D_p = 0.0
#
#
# ######################################################################################
#                 #   Lateral Coefficients
# ######################################################################################
# C_Y_0 = 0.0
# C_ell_0 = 0.0
# C_n_0 = 0.0
# C_Y_beta = -0.98
# C_ell_beta = -0.13
# C_n_beta = 0.073
# C_Y_p = 0.0
# C_ell_p = -0.51
# C_n_p = 0.069
# C_Y_r = 0.0
# C_ell_r = 0.25
# C_n_r = -0.095
# C_Y_delta_a = 0.075
# C_ell_delta_a = 0.17
# C_n_delta_a = -0.011
# C_Y_delta_r = 0.19
# C_ell_delta_r = 0.0024
# C_n_delta_r = -0.069
#
# ######################################################################################
#                 #   Propeller thrust / torque parameters (see addendum by McLain)
# ######################################################################################
# C_prop = 1.0
# S_prop = 0.2027
# k_motor = 80
# kTp = 0.
# kOmega = 0.
#
# # Prop parameters
# D_prop = 20*(0.0254)     # prop diameter in m
#
# # Motor parameters
# K_V = 145.                   # from datasheet RPM/V
# KQ = (1. / K_V) * 60. / (2. * np.pi)  # KQ in N-m/A, V-s/rad
# R_motor = 0.042              # ohms
# i0 = 1.5                     # no-load (zero-torque) current (A)
#
#
# # Inputs
# ncells = 12.
# V_max = 3.7 * ncells  # max voltage for specified number of battery cells
#
# # Coeffiecients from prop_data fit
# C_Q2 = -0.01664
# C_Q1 = 0.004970
# C_Q0 = 0.005230
# C_T2 = -0.1079
# C_T1 = -0.06044
# C_T0 = 0.09357
#
# ######################################################################################
#                 #   Calculation Variables
# ######################################################################################
# #   gamma parameters pulled from page 36 (dynamics)
# gamma = Ix * Iz - (Ixz**2)
# gamma1 = (Ixz * (Ix - Iy + Iz)) / gamma
# gamma2 = (Iz * (Iz - Iy) + (Ixz**2)) / gamma
# gamma3 = Iz / gamma
# gamma4 = Ixz / gamma
# gamma5 = (Iz - Ix) / Iy
# gamma6 = Ixz / Iy
# gamma7 = ((Ix - Iy) * Ix + (Ixz**2)) / gamma
# gamma8 = Ix / gamma
#
# #   C values defines on pag 62
# C_p_0         = gamma3 * C_ell_0      + gamma4 * C_n_0
# C_p_beta      = gamma3 * C_ell_beta   + gamma4 * C_n_beta
# C_p_p         = gamma3 * C_ell_p      + gamma4 * C_n_p
# C_p_r         = gamma3 * C_ell_r      + gamma4 * C_n_r
# C_p_delta_a    = gamma3 * C_ell_delta_a + gamma4 * C_n_delta_a
# C_p_delta_r    = gamma3 * C_ell_delta_r + gamma4 * C_n_delta_r
# C_r_0         = gamma4 * C_ell_0      + gamma8 * C_n_0
# C_r_beta      = gamma4 * C_ell_beta   + gamma8 * C_n_beta
# C_r_p         = gamma4 * C_ell_p      + gamma8 * C_n_p
# C_r_r         = gamma4 * C_ell_r      + gamma8 * C_n_r
# C_r_delta_a    = gamma4 * C_ell_delta_a + gamma8 * C_n_delta_a
# C_r_delta_r    = gamma4 * C_ell_delta_r + gamma8 * C_n_delta_r
