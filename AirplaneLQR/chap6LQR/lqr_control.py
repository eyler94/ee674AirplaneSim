"""
pid_control
    - Beard & McLain, PUP, 2012
    - Last Update:
        2/6/2019 - RWB
"""
import sys
import numpy as np
sys.path.append('..')

class lqr_control:
    def __init__(self, A, B, C, K, x_eq, y_eq, u_eq, limit, Ki, Ts, throttle_flag = False):
        ### Basic stuff needed for LQR
        self.A = A
        self.B = B
        self.C = C
        self.K = K
        self.x_eq = x_eq
        self.y_eq = y_eq
        self.u_eq = u_eq
        self.Kr = -np.linalg.inv(C@np.linalg.inv(A-B@K)@B) # Why negative?

        ### Integrator stuff
        self.Ki = Ki
        self.int_error = 0.
        self.error_d1 = 0.

        ### Saturation protection
        self.Ts = Ts
        self.limit = limit
        self.throttle_flag = throttle_flag

    def update(self, r, y, x, wrap_flag=False):

        error = r-y
        if wrap_flag is True:
            while error.item(1) >= np.pi:
                # print("Correct positive")
                error[1,0] = error.item(1) - 2.*np.pi
            while error.item(1) <= -np.pi:
                # print("Correct negative")
                error[1,0] = error.item(1) + 2.*np.pi
        self.integrateError(error)

        u_unsat = self.u_eq + self.Kr@(r-self.y_eq) - self.K@(x-self.x_eq) #- self.Ki@(self.int_error)

        u_sat = self._saturate(u_unsat)

        self.integratorAntiWindup(u_sat, u_unsat)

        return u_sat

    def integrateError(self, error):
        self.int_error = self.int_error + (self.Ts/2.)*(error + self.error_d1)
        self.error_d1 = error

    def integratorAntiWindup(self, u_sat, u_unsat):
        # if self.Ki != 0.0:
        self.int_error = self.int_error + self.Ts/self.Ki*(u_sat-u_unsat)

    def _saturate(self, u):
        # saturate u at +- self.limit
        u_sat = np.copy(u)
        for spot in range(0,2):
            if u.item(spot) >= self.limit.item(spot):
                u_sat[spot][0] = self.limit.item(spot)
            elif spot == 1 and self.throttle_flag == True:
                if u.item(spot) <= 0.0:
                    u_sat[spot][0] = 0.0
                else:
                    u_sat[spot][0] = u.item(spot)
            elif u.item(spot) <= -self.limit.item(spot):
                u_sat[spot][0] = -self.limit.item(spot)
            else:
                u_sat[spot][0] = u.item(spot)
        return u_sat