def deriv(X, t):
    x, y, z, xdot, ydot, zdot = X
    r1 = np.sqrt((x      + mu)**2 + y**2 + z**2)
    r2 = np.sqrt((x - 1. + mu)**2 + y**2 + z**2)

    term_1 = x + 2. * ydot
    term_2 = -(1.-mu) * (x + mu) / r1**3
    term_3 =     -mu  * (x - 1. + mu) / r2**3
    xddot  = term_1 + term_2 + term_3

    term_1 = -2. * xdot
    term_2 = 1. - (1.-mu)/r1**3 - mu/r2**3 
    yddot  = term_1 + y * term_2

    term_1 = (1. - mu)/r1**3 + mu/r2**3  # should be plus???
    zddot  = -z * term_1

    return np.array([xdot, ydot, zdot, xddot, yddot, zddot])


class Sat(object):
    def __init__(self, X0, T0, nu12):
        self.X0 = X0
        self.pos0 = X0[:3]
        self.v0   = X0[3:]
        self.T0 = T0
        self.nu1, self.nu2 = nu12       


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint as ODEint
from mpl_toolkits.mplot3d import Axes3D

# From "Three-Dimensional, Periodic 'Halo' Orbits,
# Kathleen Connor Howell, Celestial Mechanics 32 (1984) 53-71 

pi, twopi = np.pi, 2*np.pi
mu = 0.04

# starting points:
x0     =   [0.723268, 0.729988, 0.753700, 0.777413, 0.801125, 0.817724]
y0     = 6*[0.0]
z0     =   [0.040000, 0.215589, 0.267595, 0.284268, 0.299382, 0.313788]
xdot0  = 6*[0.0]
ydot0  =   [0.198019, 0.397259, 0.399909, 0.361870, 0.312474, 0.271306]
zdot0  = 6*[0.0]

X0s    = np.array(zip(x0, y0, z0, xdot0, ydot0, zdot0))

Thalf0s = [1.300177, 1.348532, 1.211253, 1.101099, 1.017241, 0.978653]
T0s     = [2.0*x for x in Thalf0s]

nu1s    = [1181.69,    51.07839,  4.95816,  1.101843,  0.94834,  1.10361]
nu2s    = [   0.98095, -0.90203, -0.40587, -0.420200, -1.58429, -2.09182]
nu12s   = zip(nu1s, nu2s)

n_half  = 200
fractional_times  = np.linspace(0.0, 1.0, 2*n_half+1)

rtol, atol = 1E-12, 1E-12

sats   = []
for X0, T0, nu12 in zip(X0s, T0s, nu12s):
    sat = Sat(X0, T0, nu12)
    sat.n_half  = n_half
    sat.t = sat.T0 * fractional_times
    sat.rtol, sat.atol = rtol, atol    
    sats.append(sat)

for sat in sats:
    answer, info = ODEint(deriv, sat.X0, sat.t,
                          rtol=sat.rtol, atol=sat.atol,
                          full_output = True )
    sat.answer   = answer
    sat.mid    = answer[sat.n_half]
    sat.mid    = answer[sat.n_half]
    sat.info     = info

if 1 == 1:
    xL2, xL1 = 0.74091, 1.21643  # lazy!
    fig = plt.figure(figsize=[10, 8])
    ax = fig.add_subplot(1, 1, 1, projection='3d')

    for sat in sats:
        x,  y,  z  = sat.answer.T[:3]
        ax.plot(x, y, z)

    ax.plot([0.0-mu], [0], [0], 'ob', markersize=20)
    ax.plot([1.0-mu], [0], [0], 'og', markersize=12)
    ax.plot([xL2], [0], [0], 'ok', markersize=8)
    ax.plot([xL1], [0], [0], 'ok', markersize=8)

    ax.set_xlim(0.7, 1.25)
    ax.set_ylim(-0.225, 0.225)
    ax.set_zlim(-0.15, 0.40)
    ax.text(xL1, 0, -0.05, "L1", fontsize=14, horizontalalignment='center')
    ax.text(xL2, 0, -0.05, "L2", fontsize=14, horizontalalignment='center')

    nplot    = 80
    thetas   = np.linspace(0, twopi, nplot+1)[:-1]
    azimuths = -90 + 10.0 * np.cos(thetas)

    fnames = []
    for i, azim in enumerate(azimuths):
        fname = "haloz_3D_" + str(10000+i)[1:]
        ax.elev, ax.azim = 0, azim
        plt.savefig(fname)
        fnames.append(fname)

    # tight cropping
#    for i in range(len(fnames)):
#        fname_in  = "haloz_3D_" + str(10000+i)[1:]
#        fname_out = "haloz_3D_crop_" + str(10000+i)[1:] + ".png"
#        img = plt.imread(fname_in + ".png")
#        plt.imsave(fname_out, img[200:-175, 240:-190])