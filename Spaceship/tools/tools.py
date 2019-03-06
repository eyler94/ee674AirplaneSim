import numpy as np
from math import atan2, asin

def Euler2Quaternion(phi, theta, psi):
    cth = np.cos(theta/2)
    cph = np.cos(phi/2)
    cps = np.cos(psi/2)

    sth = np.sin(theta/2)
    sph = np.sin(phi/2)
    sps = np.sin(psi/2)

    e0 = cps*cth*cph+sps*sth*sph
    e1 = cps*cth*sph-sps*sth*cph
    e2 = cps*sth*cph+sps*cth*sph
    e3 = sps*cth*cph+cps*sth*sph
    e = np.array([e0,e1,e2,e3])
    return e


def Quaternion2Euler(e):
    e0 = e.item(0)
    e1 = e.item(1)
    e2 = e.item(2)
    e3 = e.item(3)

    phi = atan2(2*(e0*e1 + e2*e3),(e0**2+e3**2-e1**2-e2**2))
    theta = asin(2*(e0*e2-e1*e3))
    psi = atan2(2*(e0*e3 + e1*e2),(e0**2+e1**2-e2**2-e3**2))
    return [phi,theta,psi]

def Quaternion2Rotation(e):
    e0 = e.item(0)
    e1 = e.item(1)
    e2 = e.item(2)
    e3 = e.item(3)

    R = np.array([[e0**2 + e1**2 - e2**2 - e3**2, 2*(e1*e2 - e0*e3), 2*(e1*e3 + e0*e2)],
                  [2*(e1*e2 + e0*e3), e0**2 - e1**2 + e2**2 - e3**2, 2*(e2*e3 - e0*e1)],
                  [2*(e1*e3 - e0*e2), 2*(e2*e3 + e0*e1), e0**2 - e1**2 - e2**2 + e3**2]
                  ])
    return R