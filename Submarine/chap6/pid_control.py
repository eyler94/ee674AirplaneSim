"""
pid_control
    - Beard & McLain, PUP, 2012
    - Last Update:
        2/6/2019 - RWB
"""
import sys
import numpy as np
sys.path.append('..')

class pid_control:
    def __init__(self, kp=0.0, ki=0.0, kd=0.0, Ts=0.01, sigma=0.05, limit=1.0):
        self.kp = kp
        self.ki = ki
        self.kd = kd
        self.Ts = Ts
        self.limit = limit
        self.beta = sigma
        self.integrator = 0.0

        self.y_dot = 0.0
        self.y_d1 = 0.0

        self.error_dot = 0.0
        self.error_d1 = 0.0

        self.error_delay_1 = 0.0
        self.error_dot_delay_1 = 0.0
        # gains for differentiator
        self.a1 = (2.0 * sigma - Ts) / (2.0 * sigma + Ts)
        self.a2 = 2.0 / (2.0 * sigma + Ts)

    def update(self, y_ref, y, reset_flag=False):

        error = y_ref - y

        if reset_flag is True:
            if error > np.pi:
                error = error - 2.*np.pi
            elif error < -np.pi:
                error = error + 2.*np.pi

        diffError = True # Set true if you want derivative to act on error instead of y


        self.integrateError(error)
        # self.differentiateError(error)
        self.differentiateY(y)

        if diffError is True:
            u_unsat = self.kp*error + self.ki*self.integrator + self.kd*self.error_dot
        else:
            u_unsat = self.kp*error + self.ki*self.integrator - self.kd*self.y_dot

        u_sat = self._saturate(u_unsat)

        self.integratorAntiWindup(u_sat, u_unsat)

        return u_sat

    def update_with_rate(self, y_ref, y, ydot, reset_flag=False):
        error = y_ref - y

        self.integrateError(error)

        u_unsat = self.kp * error + self.ki * self.integrator - self.kd * ydot

        u_sat = self._saturate(u_unsat)

        self.integratorAntiWindup(u_sat, u_unsat)

        return u_sat

    def integrateError(self, error):
        self.integrator = self.integrator + (self.Ts/2.)*(error + self.error_d1)
        self.error_d1 = error

    # def differentiateError(self, error):
    #     self.error_dot = self.beta*self.error_dot + (1.-self.beta)*((error-self.error_d1)/self.Ts)
    #     self.error_d1 = error

    def differentiateY(self, y):
        self.y_dot = self.beta*self.y_dot + (1.-self.beta)*((y-self.y_d1)/self.Ts)
        self.y_d1 = y

    def integratorAntiWindup(self, u_sat, u_unsat):
        if self.ki != 0.0:
            self.integrator = self.integrator + self.Ts/self.ki*(u_sat-u_unsat)

    def _saturate(self, u):
        # saturate u at +- self.limit
        if u >= self.limit:
            u_sat = self.limit
        elif u <= -self.limit:
            u_sat = -self.limit
        else:
            u_sat = u
        return u_sat
#
# class pi_control:
#     def __init__(self, kp=0.0, ki=0.0, Ts=0.01, limit=1.0):
#         self.kp = kp
#         self.ki = ki
#         self.Ts = Ts
#         self.limit = limit
#         self.integrator = 0.0
#         self.error_delay_1 = 0.0
#
#     def update(self, y_ref, y):
#         return u_sat
#
#     def _saturate(self, u):
#         # saturate u at +- self.limit
#         if u >= self.limit:
#             u_sat = self.limit
#         elif u <= -self.limit:
#             u_sat = -self.limit
#         else:
#             u_sat = u
#         return u_sat
#
# class pd_control_with_rate:
#     # PD control with rate information
#     # u = kp*(yref-y) - kd*ydot
#     def __init__(self, kp=0.0, kd=0.0, limit=1.0):
#         self.kp = kp
#         self.kd = kd
#         self.limit = limit
#
#     def update(self, y_ref, y, ydot):
#         return u_sat
#
#     def _saturate(self, u):
#         # saturate u at +- self.limit
#         if u >= self.limit:
#             u_sat = self.limit
#         elif u <= -self.limit:
#             u_sat = -self.limit
#         else:
#             u_sat = u
#         return u_sat