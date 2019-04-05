"""
mav_dynamics
    - this file implements the dynamic equations of motion for MAV
    - use unit quaternion for the attitude state

"""
import sys
sys.path.append('..')
import numpy as np

# load message types
from message_types.msg_state import msg_state

import parameters.aerosonde_parameters as MAV
from tools.tools import Quaternion2Rotation, Quaternion2Euler

class mav_dynamics:
    def __init__(self, Ts):
        self._ts_simulation = Ts
        # set initial states based on parameter file
        # _state is the 13x1 internal state of the aircraft that is being propagated:
        # _state = [pn, pe, pd, u, v, w, e0, e1, e2, e3, p, q, r]
        # We will also need a variety of other elements that are functions of the _state and the wind.
        # self.true_state is a 19x1 vector that is estimated and used by the autopilot to control the aircraft:
        # true_state = [pn, pe, h, Va, alpha, beta, phi, theta, chi, p, q, r, Vg, wn, we, psi, gyro_bx, gyro_by, gyro_bz]
        self._state = np.array([[MAV.pn0],  # (0)
                               [MAV.pe0],   # (1)
                               [MAV.pd0],   # (2)
                               [MAV.u0],    # (3)
                               [MAV.v0],    # (4)
                               [MAV.w0],    # (5)
                               [MAV.e0],    # (6)
                               [MAV.e1],    # (7)
                               [MAV.e2],    # (8)
                               [MAV.e3],    # (9)
                               [MAV.p0],    # (10)
                               [MAV.q0],    # (11)
                               [MAV.r0]])   # (12)
        self.u_dot = 0.
        self.v_dot = 0.
        self.w_dot = 0.
        self.e0_dot = 0.
        self.e1_dot = 0.
        self.e2_dot = 0.
        self.e3_dot = 0.
        self.p_dot = 0.
        self.q_dot = 0.
        self.r_dot = 0.
        # store wind data for fast recall since it is used at various points in simulation
        self._wind = np.array([[0.], [0.], [0.]])  # wind in NED frame in meters/sec
        self._update_velocity_data()
        # store forces to avoid recalculation in the sensors function
        self._forces = np.array([[0.], [0.], [0.]])
        self._Va = MAV.u0
        self._alpha = 0
        self._beta = 0
        # initialize true_state message
        self.msg_true_state = msg_state()

    ###################################
    # public functions
    def update_state(self, delta, wind):
        '''
            Integrate the differential equations defining dynamics, update sensors
            delta = (delta_a, delta_e, delta_r, delta_t) are the control inputs
            wind is the wind vector in inertial coordinates
            Ts is the time step between function calls.
        '''
        # get forces and moments acting on rigid body
        forces_moments = self._forces_moments(delta)

        # Integrate ODE using Runge-Kutta RK4 algorithm
        time_step = self._ts_simulation
        k1 = self._derivatives(self._state, forces_moments)
        k2 = self._derivatives(self._state + time_step/2.*k1, forces_moments)
        k3 = self._derivatives(self._state + time_step/2.*k2, forces_moments)
        k4 = self._derivatives(self._state + time_step*k3, forces_moments)
        self._state += time_step/6.0 * (k1 + 2.0*k2 + 2.0*k3 + k4)

        # normalize the quaternion
        e0 = self._state.item(6)
        e1 = self._state.item(7)
        e2 = self._state.item(8)
        e3 = self._state.item(9)
        normE = np.sqrt(e0**2+e1**2+e2**2+e3**2)
        self._state[6][0] = self._state.item(6)/normE
        self._state[7][0] = self._state.item(7)/normE
        self._state[8][0] = self._state.item(8)/normE
        self._state[9][0] = self._state.item(9)/normE

        # update the airspeed, angle of attack, and side slip angles using new state
        self._update_velocity_data(wind)

        # update the message class for the true state
        self._update_msg_true_state()

    ###################################
    # private functions
    def _derivatives(self, state, forces_moments):
        """
        for the dynamics xdot = f(x, u), returns f(x, u)
        """
        # extract the states
        pn = state.item(0)
        pe = state.item(1)
        pd = state.item(2)
        u = state.item(3)
        v = state.item(4)
        w = state.item(5)
        e0 = state.item(6)
        e1 = state.item(7)
        e2 = state.item(8)
        e3 = state.item(9)
        p = state.item(10)
        q = state.item(11)
        r = state.item(12)
        #   extract forces/moments
        fx = forces_moments.item(0)
        fy = forces_moments.item(1)
        fz = forces_moments.item(2)
        k = forces_moments.item(3)
        m = forces_moments.item(4)
        n = forces_moments.item(5)
        # print("forces and moments:", forces_moments)

        # position kinematics
        pn_dot = (e1**2+e0**2-e2**2-e3**2)*u + 2*(e1*e2-e3*e0)*v + 2*(e1*e3+e2*e0)*w
        pe_dot = 2*(e1*e2+e3*e0)*u + (e2**2+e0**2-e1**2-e3**2)*v + 2*(e2*e3-e1*e0)*w
        pd_dot = 2*(e1*e3-e2*e0)*u + 2*(e2*e3+e1*e0)*v + (e3**2+e0**2-e1**2-e2**2)*w

        # position dynamics
        mass = MAV.mass
        xg = MAV.xg
        yg = MAV.yg
        zg = MAV.zg
        self.u_dot = fx/mass + v*r - w*q + xg*(q**2 + r**2) - yg*(p*q-self.r_dot) - zg*(p*r+self.q_dot)
        self.v_dot = fy/mass + w*p - u*r + yg*(r**2 + p**2) - zg*(q*r-self.p_dot) - xg*(q*p+self.r_dot)
        self.w_dot = fz/mass + u*q - v*p + zg*(p**2 + q**2) - xg*(r*p-self.q_dot) - yg*(r*q+self.p_dot)
        # u_dot = (r*v-q*w)+fx/mass
        # v_dot = (p*w-r*u)+fy/mass
        # w_dot = (q*u-p*v)+fz/mass

        # rotational kinematics
        e0_dot = 0.5*(-p*e1-q*e2-r*e3)
        e1_dot = 0.5*(p*e0+r*e2-q*e3)
        e2_dot = 0.5*(q*e0-r*e1+p*e3)
        e3_dot = 0.5*(r*e0+q*e1-p*e2)

        # rotatonal dynamics
        Ix = MAV.Ix
        Iy = MAV.Iy
        Iz = MAV.Iz
        self.p_dot = -((Iz-Iy)*q*r-mass*(yg*(self.w_dot-u*q+v*p)-zg*(self.v_dot-w*p+u*r)))/Ix + k/Ix
        self.q_dot = -((Ix-Iz)*r*p-mass*(zg*(self.u_dot-v*r+w*q)-xg*(self.w_dot-u*q+v*p)))/Iy + m/Iy
        self.r_dot = -((Iy-Ix)*p*q-mass*(xg*(self.v_dot-w*p+u*r)-yg*(self.u_dot-v*r+w*q)))/Iz + n/Iz
        # p_dot = MAV.gamma1*p*q-MAV.gamma2*q*r + MAV.gamma3*k+MAV.gamma4*n
        # r_dot = MAV.gamma7*p*q - MAV.gamma1*q*r + MAV.gamma4*k+MAV.gamma8*n
        # q_dot = MAV.gamma5*p*r - MAV.gamma6*(p**2-r**2) + m/Iy

        # collect the derivative of the states
        x_dot = np.array([[pn_dot, pe_dot, pd_dot, self.u_dot, self.v_dot, self.w_dot,
                           e0_dot, e1_dot, e2_dot, e3_dot, self.p_dot, self.q_dot, self.r_dot]]).T
        return x_dot

    def _update_velocity_data(self, wind=np.zeros((6,1))):
        # Split wind into components
        self._ur = self._state.item(3)-wind.item(0)  # u - uw
        self._vr = self._state.item(4)-wind.item(1)  # v - vw
        self._wr = self._state.item(5)-wind.item(2)  # w - ww
        # compute airspeed
        self._Va = np.sqrt(self._ur**2 + self._vr**2 + self._wr**2)
        # compute angle of attack
        self._alpha = np.arctan(self._wr/self._ur)
        # compute sideslip angle
        self._beta = np.arcsin(self._vr/self._Va)

    def _forces_moments(self, delta):
        """
        return the forces on the UAV based on the state, wind, and control surfaces
        :param delta: np.matrix(delta_a, delta_e, delta_r, delta_t)
        :return: Forces and Moments on the UAV np.matrix(Fx, Fy, Fz, Ml, Mn, Mm)
        """
        assert delta.shape == (5,1)
        dspl = delta[0,0] # Left Stern Plane when looking stern to bow
        dspr = delta[1,0] # Right Stern Plane
        drt = delta[2,0] # Top Rudder
        drb = delta[3,0] # Bottom Rudder
        dm = delta[4,0] # Motor

        pn = self._state.item(0)
        pe = self._state.item(1)
        pd = self._state.item(2)
        u = self._state.item(3)
        v = self._state.item(4)
        w = self._state.item(5)
        e0 = self._state.item(6)
        e1 = self._state.item(7)
        e2 = self._state.item(8)
        e3 = self._state.item(9)
        p = self._state.item(10)
        q = self._state.item(11)
        r = self._state.item(12)


        ##Forces and Moments
        #Hydrostatic Forces
        phi, theta, psi = Quaternion2Euler(self._state[6:10])
        xg = MAV.xg
        yg = MAV.yg
        zg = MAV.zg
        Weight = MAV.mass*MAV.gravity
        Bouyancy = MAV.rho*MAV.gravity*MAV.length*MAV.radius**2*np.pi*self._submerged()
        Xhs = -(Weight-Bouyancy)*np.sin(theta)
        Yhs = (Weight-Bouyancy)*np.cos(theta)*np.sin(phi)
        Zhs = (Weight-Bouyancy)*np.cos(theta)*np.cos(phi)
        Khs = -MAV.yg * Weight * np.cos(theta) * np.cos(phi) - zg * np.cos(theta) * np.sin(phi)
        Mhs = -zg*Weight*np.sin(theta) - xg*Weight*np.cos(theta)*np.cos(phi)
        Nhs = -yg*Weight*np.cos(theta)*np.sin(phi) - zg*Weight*np.sin(theta)

        #Added Mass
        Xa = Xudot*self.u_dot + Xwq*w*q + Wqq*q**2 + Xvr*v*r + Xrr*r**2
        Ya = Yvdot*self.v_dot + Yrdot*self.r_dot + Yura*u*r + Ywp*w*p + Ypq*p*q
        Za = Zwdot*self.w_dot + Zqdot*self.q_dot + Zuqa*u*q + Zvp*v*p + Zrp*r*p
        Ka = Kpdot*self.p_dot
        Ma = Mwdot*self.w_dot + Mqdot*self.q_dot + Muwa*u*w + Mvp*v*p + Mrp*r*p + Muqa*u*q
        Na = Nvdot*self.v_dot + Nrdot*self.r_dot + Nuva*u*v + Nwp*w*p + Mpq*p*q + Nura*u*r


        #Hydrodynamic forces

        #Control surfaces

        #Forces and Moments summation
        fx = Xhs
        fy = Yhs
        fz = Zhs
        Mx = Khs
        My = Mhs
        Mz = Nhs

        self._forces[0] = fx
        self._forces[1] = fy
        self._forces[2] = fz
        return np.array([[fx, fy, fz, Mx, My, Mz]]).T

    def _update_msg_true_state(self):
        # update the class structure for the true state:
        #   [pn, pe, h, Va, alpha, beta, phi, theta, chi, p, q, r, Vg, wn, we, psi, gyro_bx, gyro_by, gyro_bz]
        phi, theta, psi = Quaternion2Euler(self._state[6:10])
        self.msg_true_state.pn = self._state.item(0)
        self.msg_true_state.pe = self._state.item(1)
        self.msg_true_state.h = -self._state.item(2)
        self.msg_true_state.Va = self._Va
        self.msg_true_state.alpha = self._alpha
        self.msg_true_state.beta = self._beta
        self.msg_true_state.phi = phi
        self.msg_true_state.theta = theta
        self.msg_true_state.psi = psi
        self.msg_true_state.Vg = np.sqrt((self._state.item(3))**2+(self._state.item(4))**2+(self._state.item(5))**2)
        Vg = np.array([self._state.item(3),self._state.item(4),self._state.item(5)])
        Vg_M = np.linalg.norm(Vg)
        Vg_h = np.array([self._state.item(3),self._state.item(4),0.])
        Vg_h_M = np.linalg.norm(Vg_h)
        Va_h = np.array([self._ur,self._vr,0.])
        Va_h_M = np.linalg.norm(Va_h)
        self.msg_true_state.gamma = np.arccos(Vg.dot(Vg_h)/(Vg_M*Vg_h_M))
        num = Vg_h.dot(Va_h)
        den = (Vg_h_M*Va_h_M)
        frac = np.round(num/den,8)
        self.msg_true_state.chi = psi + self._beta + np.arccos(frac)
        self.msg_true_state.p = self._state.item(10)
        self.msg_true_state.q = self._state.item(11)
        self.msg_true_state.r = self._state.item(12)
        self.msg_true_state.wn = self._wind.item(0)
        self.msg_true_state.we = self._wind.item(1)

    def _submerged(self):
        depth = self._state.item(2)
        if depth >= MAV.radius*2.:
            print("deep")
            return 1.
        elif depth >= 0:
            print("shallow")
            return depth/MAV.radius*2.
        else:
            print("outta water")
            return 0.