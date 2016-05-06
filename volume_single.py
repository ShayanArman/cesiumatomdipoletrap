# volume_single.py

import math                    # http://docs.python.org/library/math.html
import numpy                   # numpy.scipy.org/
import datetime                # to measure the time consumed for simulation
# import matplotlib
# from mpl_toolkits.mplot3d import Axes3D
# from matplotlib import mlab
# from matplotlib import pyplot


def ax(x, beta):
    # Formula for the acceleration in the x direction here
    return (beta*x)


def ay(y, beta, rpow2, yratio):
    # Formula for the acceleration in the x direction here
    if (y > 0):
        return (beta*((w0pow2-rpow2/yratio)/(yR**2))*y-g)
    else:
        return (-g)
    # return (beta*((w0pow2-rpow2/yratio)/(yR**2))*y)
    # return (-g)


def az(z, beta):
    # Formula for the acceleration in the x direction here
    return (beta*z)


# Speed of light [m/s]
c = 2.9979e8
# Atomic mass of cesium [kg]
m = 2.2069e-25
global g
# Gravitational acceleration at latitude 43.47 degree [m/s^2]
g = 9.805
# Wavelength of dipole trap laser [m]
wavelength = 9.356e-7
# Light frequency of dipole trap 936.5 nm [/s]
omega = 2*numpy.pi*c/wavelength
# Light frequency of D1 [/s]
omegaD1 = 2*numpy.pi*c/894.59295986e-9
# Light frequency of D2 [/s]
omegaD2 = 2*numpy.pi*c/852.34727582e-9
# Natural linewidth of D1 [/s]
GammaD1 = 2*numpy.pi*4.5612e6
# Natural linewidth of D2 [/s]
GammaD2 = 2*numpy.pi*5.2227e6

# U_D1 = alphaD1 * I
alphaD1 = -1.5*numpy.pi*(c**2)*(omegaD1**(-3))*(GammaD1/(omegaD1-omega)+GammaD1/(omegaD1+omega))
# U_D2 = alphaD2 * I
alphaD2 = -1.5*numpy.pi*(c**2)*(omegaD2**(-3))*(GammaD2/(omegaD2-omega)+GammaD2/(omegaD2+omega))
# a = alpha * (dI/dr)
alpha = -((1/3)*alphaD1 + (2/3)*alphaD2)/m

# Laser power [W]
P = 0.1
# Beam waist at the point of focus = mode field radius [m]
w0 = 2.75e-6
global w0pow2
# Sqare of beam waist [m^2]
w0pow2 = w0*w0
# Peak intensity at the point of focus [W/^2]
I0 = 2*P/(numpy.pi*w0pow2)
global yR
# Rayleigh distance [m]
yR = numpy.pi*w0pow2/wavelength

# Initial values
# x0 = 0.0 # initial position in x [m]
# initial position in x [m]
pos_vec = [0, 0.007, 0]
x0 = pos_vec[0]
# initial position in y [m]
y0 = pos_vec[1]
# initial position in z [m]
z0 = pos_vec[2]
v_vec = [-0.06786763250028102, 0.04507199789615998, 0.0013829422806097343]
# initial velocity in x [m]
vx0 = v_vec[0]
# initial velocity in y [m]
vy0 = v_vec[1]
# initial velocity in z [m]
vz0 = v_vec[2]
# time [s]
t = 0
# time step [s]
dt = 0.00001
# dt = 0.000001 # time step [s]
# vstep = [0.001, 0.0001] # velocity step [m]
# Square of beam radius at position y
wypow2 = 0
# Square of atom's distance from axis
rpow2 = 0
beta = 0
vxlist1 = []
vylist1 = []
vzlist1 = []

# fh = open("data.txt","w")
vxflag = 1
vzflag = 1

t1 = datetime.datetime.now()
t = 0
dt = 0.00001
x = x0
y = y0
z = z0
vx = vx0
vy = vy0
vz = vz0
result = "TRUE"
while (y > -0.0003):
    x = x + vx*dt
    y = y + vy*dt
    z = z + vz*dt

    # Square of atom's distance from axis
    rpow2 = x**2 + z**2
    yratio = 1+(y/yR)**2
    # Square of beam radius at position y
    wypow2 = w0pow2 * yratio
    # Peak intensity at position y
    I0y = I0/yratio
    if (y > 0):
        beta = -alpha*(2*I0y/wypow2)*numpy.exp(-rpow2/wypow2)
    else:
        # 5.625e-11 = (7.5um)^2
        if (rpow2 > 5.625e-11):
            result = "FALSE"
            # Break the while loop
            y = -0.002
            vxlist1.append(vx0)
            vylist1.append(vy0)
            vzlist1.append(vz0)
            vxflag = 0
        else:
            beta = -alpha*(2*I0/w0pow2)*numpy.exp(-rpow2/w0pow2)

    ax_val = ax(x, beta)
    vx = vx + dt*ax_val
    # vx = vx + dt*ax(x,beta)

    ay_val = ay(y, beta, rpow2, yratio)
    vy = vy + dt*ay_val
    # vy = vy + dt*ay(y,beta,rpow2,yR)

    az_val = az(z, beta)
    vz = vz + dt*az_val
    # vz = vz + dt*az(z,beta)

    if (y < 0.0035):
        dt = 0.000001

    t = t + dt

t2 = datetime.datetime.now()
print result
print(t2 - t1)

for k in range(0, len(vzlist1)):
    print('vx={0:.6f}\tvz={1:.6f}\tvy={2:.6f}'.format(vxlist1[k], vzlist1[k], vylist1[k]))

# fig = pyplot.figure()
# ax = fig.gca(projection='3d')
# ax.plot_surface(vxlist,vzlist,vylist)
# pyplot.show()
# fh.close()
