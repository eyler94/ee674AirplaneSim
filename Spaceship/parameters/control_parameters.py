import sys
sys.path.append('..')
import numpy as np
# import chap5.transfer_function_coef as TF
import parameters.aerosonde_parameters as MAV

gravity = MAV.gravity
sigma = 0.05
Va0 = MAV.Va0

#----------roll loop-------------
wn_phi = 5.
zeta_phi = 0.8
a_phi_2 = 130.6
a_phi_1 = 22.6

roll_kp = wn_phi**2./a_phi_2
roll_kd = (2.*zeta_phi*wn_phi-a_phi_1)/a_phi_2

#----------course loop-------------
zeta_chi = 0.707
wn_chi = 1.

course_kp = 2.*zeta_chi*wn_chi*Va0/gravity
course_ki = wn_chi**2.*Va0/gravity

#----------sideslip loop-------------
sideslip_ki = 0 # Fix this later
sideslip_kp = 1 #Fix this later, too

# #----------yaw damper-------------
# yaw_damper_tau_r =
# yaw_damper_kp =

#----------pitch loop-------------
wn_theta = 1
zeta_theta = 0.707

a_theta_1 = 5.288
a_theta_2 = 99.7#
# #----------pitch loop-------------
# pitch_kp =
# pitch_kd =
# K_theta_DC =
#
# #----------altitude loop-------------
# altitude_kp =
# altitude_ki =
# altitude_zone =
#
# #---------airspeed hold using throttle---------------
# airspeed_throttle_kp =
# airspeed_throttle_ki =
a_theta_3 = -36.02

pitch_kp = (wn_theta**2.-a_theta_2)/a_theta_3
pitch_kd = (2.*zeta_theta*wn_theta - a_theta_1)/a_theta_3
K_theta_DC = pitch_kp*a_theta_3/wn_theta**2.

#----------altitude loop-------------
zeta_h = 0.707
wn_h = 1

altitude_kp = (2.*zeta_h*wn_h)/(K_theta_DC*Va0)
altitude_ki = (wn_h**2.)/(K_theta_DC*Va0)
altitude_zone = 10

#---------airspeed hold using throttle---------------
wn_v = 1
zeta_v = 0.707
a_v_1 = 0.1888
a_v_2 = 1.896

airspeed_throttle_kp = (2.*zeta_v*wn_v-a_v_1)/a_v_2
airspeed_throttle_ki = wn_v**2./a_v_2
