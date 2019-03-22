import numpy as np
from math import sin, cos, atan, atan2
import sys

sys.path.append('..')
from message_types.msg_autopilot import msg_autopilot

class path_follower:
    def __init__(self):
        self.chi_inf = np.radians(70)  # approach angle for large distance from straight-line path
        self.k_path = 0.6  # proportional gain for straight-line path following
        self.k_orbit = 0.6  # proportional gain for orbit following
        self.gravity = 9.8
        self.autopilot_commands = msg_autopilot()  # message sent to autopilot

    def update(self, path, state):
        if path.flag=='line':
            self._follow_straight_line(path, state)
        elif path.flag=='orbit':
            self._follow_orbit(path, state)
        return self.autopilot_commands

    def _follow_straight_line(self, path, state):
        chi_q = np.arctan2(path.line_direction.item(1),path.line_direction.item(0))
        cc = np.cos(chi_q)
        sc = np.sin(chi_q)
        Rpi = np.array([[cc, sc, 0.],\
                        [-sc, cc, 0.],\
                        [0., 0., 1.]])
        p = np.array([[state.pn, state.pe, -state.h]]).T
        ep = Rpi@(p-path.line_origin)
        epy = ep.item(1)
        self.autopilot_commands.airspeed_command = 25.0
        self.autopilot_commands.course_command = chi_q - self.chi_inf*2./np.pi*np.arctan(self.k_path*epy)
        self.autopilot_commands.altitude_command = 100.
        self.autopilot_commands.phi_feedforward = 0.

    def _follow_orbit(self, path, state):
        self.autopilot_commands.airspeed_command = 25.0
        self.autopilot_commands.course_command = 0.
        self.autopilot_commands.altitude_command = 100.
        self.autopilot_commands.phi_feedforward = 0.

    def _wrap(self, chi_c, chi):
        while chi_c-chi > np.pi:
            chi_c = chi_c - 2.0 * np.pi
        while chi_c-chi < -np.pi:
            chi_c = chi_c + 2.0 * np.pi
        return chi_c

