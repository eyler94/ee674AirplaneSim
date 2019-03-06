"""
compute_ss_model
    - Chapter 5 assignment for Beard & McLain, PUP, 2012
    - Update history:  
        2/4/2019 - RWB
"""
import sys
sys.path.append('..')
import numpy as np
from scipy.optimize import minimize
from tools.tools import Euler2Quaternion, Quaternion2Euler
# from tools.transfer_function import transfer_function
import parameters.aerosonde_parameters as MAV
from parameters.simulation_parameters import ts_simulation as Ts
from control import matlab
import parameters.aerosonde_parameters as MAV


def compute_tf_model(mav, trim_state, trim_input):
    # trim values
    rho = MAV.rho
    Va = trim_state.item(3)
    Vg = np.sqrt(trim_state.item(3)**2 + trim_state.item(4)**2 + trim_state.item(5)**2)
    S = MAV.S_wing
    q_bar = -0.5*rho*Va**2.0*S

    b = MAV.b

    a_phi_1 = q_bar*b*MAV.C_p_p*b/2./Va
    a_phi_2 = -q_bar*b*MAV.C_p_delta_a

    T_phi_delta_a = matlab.tf([a_phi_2],[1., a_phi_1, 0.])

    T_chi_phi = matlab.tf([MAV.gravity/Vg],[1., 0.])

    a_beta_1 = -rho*Va*S/2/MAV.mass*MAV.C_Y_beta
    a_beta_2 = rho*Va*S/2/MAV.mass*MAV.C_Y_delta_r

    T_beta_delta_r = matlab.tf([a_beta_2],[1., a_beta_1])

    c = MAV.c
    a_theta_1 = q_bar*c/MAV.Jy*MAV.C_m_q*c/2./Va
    a_theta_2 = q_bar*c/MAV.Jy*MAV.C_m_alpha
    a_theta_3 = -q_bar * c / MAV.Jy * MAV.C_m_delta_e

    T_theta_delta_e = matlab.tf([a_theta_3],[1., a_theta_1, a_theta_2])

    T_h_theta = matlab.tf([Va],[1., 0])

    e = np.array([trim_state.item(6), trim_state.item(7), trim_state.item(8), trim_state.item(9)])
    [phi,theta,psi] = Quaternion2Euler(e)

    T_h_Va = matlab.tf([theta],[1., 0])

    a_V_1 = rho*Va*S/MAV.mass*(MAV.C_D_0+(MAV.C_D_alpha*mav._alpha)+(MAV.C_D_delta_e*trim_input.item(1)))+rho*MAV.S_prop/MAV.mass*MAV.C_prop*Va
    a_V_2 = rho*MAV.S_prop/MAV.mass*MAV.C_prop*MAV.k_motor**2.*trim_input.item(1)
    a_V_3 = MAV.gravity*np.cos(theta) # Didn't include chi_star because in trim, chi_star should be zero.

    T_Va_delta_t = matlab.tf([a_V_2],[1., a_V_1]) # Didn't include dt_bar in the num of the tf because in trim it should be zero.

    T_Va_theta = matlab.tf([-a_V_3],[1., a_V_1]) # Didn't include d_theta_bar in the num of the tf because in trim it should be zero.

    tf = open('./tf.txt','w')

    tf.write("T_phi_delta_a: " + str(T_phi_delta_a) + "\n")
    tf.write("T_chi_phi: " + str(T_chi_phi) + "\n")
    tf.write("T_theta_delta_e: " + str(T_theta_delta_e) + "\n")
    tf.write("T_h_theta: " + str(T_h_theta) + "\n")
    tf.write("T_h_Va: " + str(T_h_Va) + "\n")
    tf.write("T_Va_delta_t: " + str(T_Va_delta_t) + "\n")
    tf.write("T_Va_theta: " + str(T_Va_theta) + "\n")
    tf.write("T_beta_delta_r: " + str(T_beta_delta_r) + "\n")

    tf.close()

    return T_phi_delta_a, T_chi_phi, T_theta_delta_e, T_h_theta, T_h_Va, T_Va_delta_t, T_Va_theta, T_beta_delta_r

def compute_ss_model(mav, trim_state, trim_input):
    x_euler = euler_state(trim_state)

    f_x = df_dx(mav,trim_state,trim_input)
    f_u = df_du(mav,trim_state,trim_input)

    dt_dq = np.zeros([12,13])
    dt_dq[0:6,0:6] = np.eye(6)
    dt_dq[9:,10:]=np.eye(3)
    # dt_dq[6:9,6:10]=

    return A_lon, B_lon, A_lat, B_lat

def euler_state(x_quat):
    # convert state x with attitude represented by quaternion
    # to x_euler with attitude represented by Euler angles
    e = np.array([x_quat.item(6), x_quat.item(7), x_quat.item(8), x_quat.item(9)])
    [phi, theta, psi] = Quaternion2Euler(e)

    x_euler = np.array([[x_quat.item(0)],
                        [x_quat.item(1)],
                        [x_quat.item(2)],
                        [x_quat.item(3)],
                        [x_quat.item(4)],
                        [x_quat.item(5)],
                        [phi],
                        [theta],
                        [psi],
                        [x_quat.item(10)],
                        [x_quat.item(11)],
                        [x_quat.item(12)]
                        ])

    return x_euler

def quaternion_state(x_euler):
    # convert state x_euler with attitude represented by Euler angles
    # to x_quat with attitude represented by quaternions

    e = Euler2Quaternion(x_euler.item(6),x_euler.item(7),x_euler.item(8))
    x_quat = np.array([[x_euler.item(0)],
                        [x_euler.item(1)],
                        [x_euler.item(2)],
                        [x_euler.item(3)],
                        [x_euler.item(4)],
                        [x_euler.item(5)],
                        [e.item(0)],
                        [e.item(1)],
                        [e.item(2)],
                        [e.item(3)],
                        [x_euler.item(10)],
                        [x_euler.item(11)],
                        [x_euler.item(12)]
                        ])

    return x_quat

def f_euler(mav, x_euler, input):
    # return 12x1 dynamics (as if state were Euler state)
    # compute f at euler_state
    return f_euler

def df_dx(mav, x_euler, input):
    # take partial of f_euler with respect to x_euler
    eps = 0.01
    A = np.zeros((12,12))
    f_m = mav._forces_moments(input)
    f_at_x_quat = mav._derivatives(x_euler,f_m)
    f_at_x = euler_state(f_at_x_quat)
    for i in range(0,12):
        x_eps = np.copy(x_euler)
        x_eps[i][0] += eps
        f_at_x_eps_quat = mav._derivatives(x_eps,f_m)
        f_at_x_eps = euler_state(f_at_x_eps_quat)
        df_dxi = (f_at_x_eps - f_at_x)/eps
        A[:,i] = df_dxi[:,0]
    return A

def df_du(mav, x_euler, delta):
    # take partial of f_euler with respect to delta
    eps = 0.01
    B = np.zeros((12,4))

    return B

def dT_dVa(mav, Va, delta_t):
    # returns the derivative of motor thrust with respect to Va
    return dThrust

def dT_ddelta_t(mav, Va, delta_t):
    # returns the derivative of motor thrust with respect to delta_t
    return dThrust