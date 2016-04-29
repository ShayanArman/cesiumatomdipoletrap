# %matplotlib inline
# import matplotlib.pyplot as pyplot
import numpy as np
import math
import random


# Initial values
# Speed of light [m/s]
light_speed = 2.9979e8
# Gravitational acceleration at latitude 43.47 degree [m/s^2]
gravity_acc = 9.805
# Atomic mass of cesium [kg]
mass = 2.2069e-25
# Boltzmann Constant KB
KB = 1.38065e-23
# Temperature
T = 2e-5

MAX_BOLT_MU = 0
MAX_BOLT_SIGMA = math.sqrt(KB * T / mass)

# Minimum and maximum velocity of Cesium atoms
v_min = 0
v_max = 0.25

# Wavelength of dipole trap laser [m]
wavelength = 9.356e-7
# Light frequency of dipole trap; 936.5 nm [/s]
omega = 2*np.pi*light_speed/wavelength
# Light frequency of D1 [/s]
omega_d1 = 2*np.pi*light_speed/894.59295986e-9
# Light frequency of D2 [/s]
omega_d2 = 2*np.pi*light_speed/852.34727582e-9
# Natural linewidth of D1 [/s]
gamma_d1 = 2*np.pi*4.5612e6
# Natural linewidth of D2 [/s]
gamma_d2 = 2*np.pi*5.2227e6

# Laser power [W]
Po = 0.1
# Beam waist at the point of focus = mode field radius [m]
waist_size_init = 2.75e-6
# time step [s]
dt = 0.000001

alpha_d1 = (1/(omega_d1**3))*(
    gamma_d1*((1/(omega_d1 - omega)) + (1/(omega_d1 + omega))))
alpha_d2 = (1/(omega_d2**3))*(
    gamma_d2*((1/(omega_d2 - omega)) + (1/(omega_d2 + omega))))
alpha = ((3/2.0)*math.pi*(light_speed**2.0))*(
    ((1/3.0)*alpha_d1 + (2/3.0)*alpha_d2) / mass)
yr = math.pi * (waist_size_init**2) / wavelength


class Cesium(object):
    def __init__(self, init_pos_vec, init_vel_vec):
        self.position = init_pos_vec
        self.velocity = init_vel_vec
        # If result is true, then atom at init pos and init vel will make it
        # into the fiber.
        self.result = False
        self.position_trace = [
            [self.position[0]], [self.position[1]], [self.position[2]]]

    def acc(self):
        """
        Returns an acceleration vector for the particle in position x, y, z

        """
        pos_intensity_vector = self.intensity_gradient()

        acceleration = [(alpha * i) for i in pos_intensity_vector]
        # Add gravity into the acceleration in the Y direction
        acceleration[1] += -gravity_acc
        return acceleration

    def intensity_gradient(self):
        """
        Formula for intensity gradient at point x, y, z.

        """
        num = (self.position[0]**2 + self.position[2]**2)
        beam_width_temp = (1 + (self.position[1]/yr)**2)
        beam_width = (waist_size_init**2) * beam_width_temp
        exp_parameter = -num / (beam_width)
        gaussian_dist = math.exp(exp_parameter)
        temp = -(4*Po*gaussian_dist) / (math.pi * beam_width)

        gradient_x = temp * (self.position[0] / beam_width)
        gradient_y = (
            temp *
            (self.position[1] / ((yr**2) * (beam_width_temp))) *
            (1 + exp_parameter))
        gradient_z = temp * (self.position[2] / beam_width)
        return [gradient_x, gradient_y, gradient_z]

    def update_position(self):
        """
        Update the position of the atom based on its velocity.

        """
        # position_x = position_x + velocity_x*dt
        # position_x = position_x + position_delta[0]
        position_delta = map(lambda v_direction: v_direction*dt, self.velocity)
        self.position = map(
            lambda pos_init, pos_delta: pos_init + pos_delta, self.position,
            position_delta)

    def update_velocity(self):
        """
        Update velocity of the atom based on its acceleration.

        """
        # velocity_x = velocity_x + acceleration_x*dt
        # velocity_x = velocity_x + vel_delta[0]
        acceleration_vector = self.acc()
        vel_delta = map(
            lambda acc_direction: acc_direction*dt,
            acceleration_vector)
        self.velocity = map(
            lambda v_init, v_delta: v_init + v_delta, self.velocity, vel_delta)
        return acceleration_vector

    def update_position_trace(self, time, acceleration):
        self.position_trace[0].append(self.position[0])
        self.position_trace[1].append(self.position[1])
        self.position_trace[2].append(self.position[2])

        print 't={time}\tx={x}\tvx={vx}\tax={ax}\ty={y}\tvy={vy}\tay={ay}'.format(
            time=time, x=round(self.position[0], 6), y=round(self.position[1], 6),
            vx=round(self.velocity[0], 4), vy=round(self.velocity[1], 4),
            ax=round(acceleration[0], 4), ay=round(acceleration[1], 4))

    def run(self, trace_position=False):
        time = 0

        while self.position[1] > 0:
            # Update position based on current velocity
            self.update_position()
            # Update velocity based on new position in space
            acc_v = self.update_velocity()

            if trace_position:
                # Record the positions in order to graph the spiral
                self.update_position_trace(time, acc_v)

            time += dt
        distance_from_zero = math.sqrt(
            self.position[0]**2 + self.position[2]**2)
        if distance_from_zero < .00014:
            self.result = True
        return time


class Tests:
    def maxwell_boltzmann_test(self):
        """
        Given a single position, it goes through a number of random velocity
        vectors, and returns the percentage of atoms that reaches the fiber.

        """
        init_pos_vector = [0.00005, 0.05, 0.00034]
        num_successes = 0
        num_trials = 150
        for x in range(num_trials):
            velocity = np.random.normal(MAX_BOLT_MU, MAX_BOLT_SIGMA, 3).tolist()
            c = Cesium(init_pos_vector, velocity)
            c.run()

            if c.result is True:
                num_successes += 1

        # Return percentage of atoms that made it into the fiber at this position.
        return (num_successes / num_trials) * 100


Tests().maxwell_boltzmann_test()
# init_pos_vector = [0.00005, 0.05, 0.00034]
# init_vel_vector = [0.0, 0.0, 0.0]
# c = Cesium(init_pos_vector, init_vel_vector)
# time = c.run(trace_position=True)

# PLOTTING
# fig = pyplot.figure()
# ax = pyplot.axes(projection='3d')
# ax.plot3D(position_trace[0], position_trace[2], position_trace[1], 'red')
# pyplot.show()

# INITIAL CESIUM
# print('t = {0:.5f} s'.format(time))
# print('x = {0:.5f} m'.format(c.position[0]))
# print('y = {0:.5f} m'.format(c.position[1]))
# print('z = {0:.5f} m'.format(c.position[2]))
# print('vx = {0:.5f} m/s'.format(c.velocity[0]))
# print('vy = {0:.5f} m/s'.format(c.velocity[1]))
# print('vz = {0:.5f} m/s'.format(c.velocity[2]))
