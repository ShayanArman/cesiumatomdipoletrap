# %matplotlib inline
import math                    # http://docs.python.org/library/math.html
import numpy as np                   # numpy.scipy.org/
import matplotlib.pyplot as pyplot # matplotlib.sourceforge.net
from mpl_toolkits import mplot3d


# Initial values
t = 0 # time [s]
dt = 0.0001 # time step [s]
mass = 10e-9
waist_size_init = 2.75e-6
lambda_dip = 935.6e-9
lambda_d1 = 894.59e-9
lambda_d2 = 852.35e-9
gamma_d1 = 2*math.pi*4.56e+6
gamma_d2 = 2*math.pi*5.22e+6
alpha = -((3/2)*math.pi*2.998e+8)*(1/w_d1**3)*((gamma_d1/(w_d1 + w)) + (gamma_d1/(w_d1 + w)))
Po = 10
yr = math.pi * (waist_size_init**2)/lambda_

class Cesium(object):
    def __init__(self, init_pos_vec, init_vel_vec):
        self.position = init_pos_vec
        self.velocity = init_vel_vec

    def acc(self, x, y, z):
        """
        Returns an acceleration vector for the particle in position x, y, z

        """
        pos_intensity_vector = self.intensity_gradient(x, y, z)
        alpha_mass = -alpha/mass

        acceleration = [(alpha_mass * i) for i in pos_intensity_vector]
        # Add gravity into the acceleration in the Y direction
        acceleration[1] += -9.805
        return acceleration

    def intensity_gradient(self, x, y, z):
        """
        Formula for intensity gradient at point x, y, z.

        """
        num = (x**2 + z**2)
        b_w_temp = (1 + (y/yr)**2)
        beam_width = (waist_size_init**2) * b_w_temp
        exp_parameter = -num / (beam_width)
        gaussian_dist = math.exp(exp_parameter)
        temp = -(4*Po*gaussian_dist) / (math.pi * beam_width)

        gradient_x = temp * (x / beam_width)
        gradient_y = temp * (y / ((yr**2) * (b_w_temp))) * (1 + exp_parameter)
        gradient_z = temp * (z / beam_width)

        return [gradient_x, gradient_y, gradient_z]

    def update_position(self, time_delta):
        """
        Update the position of the atom based on its velocity.

        """
        # position_x = position_x + velocity_x*dt
        # position_x = position_x + position_delta
        position_delta = map(lamda v_direction: v_direction*time_delta, self.velocity)
        self.position = map(lambda x, y: x + y, position_delta, self.position)


init_pos_vector = [0.0, 0.0, 0.0]
init_vel_vector = [0.0, 0.0, 0.0]
c = Cesium(init_pos_vector, init_vel_vector)

position_trace = [[init_pos_vector[0]], [init_pos_vector[1]], [init_pos_vector[2]]]
def update_position_trace(new_positions):
    position_trace[0].append(new_position[0])
    position_trace[1].append(new_position[1])
    position_trace[2].append(new_position[2])

for i in range(1,1200):
    self.update_position(dt)
    # Update position based on current velocity
    update_position_trace(self.position)
    # Update velocity based on new position in space

    acceleration = c.acc(c.x, c.y, c.z)
    c.vx = c.vx + acceleration[x]*dt;
    c.vy = c.vy + acceleration[y]*dt;
    c.vz = c.vz + acceleration[z]*dt;
    if c.y <= 0 and i > 1:
        break
    t = t + dt;

fig = pyplot.figure()
ax = pyplot.axes(projection='3d')
ax.plot3D(x_values, z_values, y_values, 'red')
pyplot.show()

print('t = {0:.5f} s'.format(t))
print('x = {0:.5f} m'.format(x))
print('y = {0:.5f} m'.format(y))
print('z = {0:.5f} m'.format(z))
print('vx = {0:.5f} m/s'.format(vx))
print('vy = {0:.5f} m/s'.format(vy))
print('vz = {0:.5f} m/s'.format(vz))

