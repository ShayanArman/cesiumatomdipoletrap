# %matplotlib inline
# import matplotlib.pyplot as pyplot
import numpy as np
import math
import pdb
import volume_single

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
FIBER_RADIUS_SQUARED = (7.5e-6)**2
# FIBER_BUFFER_DISTANCE = 0.00001
FIBER_BUFFER_DISTANCE = 0.00003

# Laser power [W]
Po = 0.1
# Beam waist at the point of focus = mode field radius [m]
waist_size_init = 2.75e-6
# time step [s]
DELTA_TIME = 0.000001
# DELTA_TIME = 0.00001

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
        self.result = True
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
        position_delta = map(lambda v_direction: v_direction*DELTA_TIME, self.velocity)
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
            lambda acc_direction: acc_direction*DELTA_TIME,
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
        """
        Given the initial position of the atom. Iterate over delta time,
        updating its position, and decide if the atom ends up inside the fiber.

        """
        time = 0
        dt = DELTA_TIME
        while self.position[1] > -FIBER_BUFFER_DISTANCE:
            # Update position based on current velocity
            self.update_position()
            # Update velocity based on new position in space
            acc_v = self.update_velocity()

            if trace_position:
                # Record the positions in order to graph the spiral
                self.update_position_trace(time, acc_v)
            if self.position[1] < 0 and (
                    self.position[0]**2 + self.position[2]**2 > FIBER_RADIUS_SQUARED):
                self.result = False
                break

            if (self.position[1] < 0.0035):
                dt = 0.000001

            time += dt
        return time


class Tests:
    def single_atom_test(self, init_pos, init_vel, trace_position=True):
        """
        Test whether a Cesium atom with the init position and init velocity
        will make it into the fiber.

        """
        c = Cesium(init_pos, init_vel)
        c.run(trace_position)
        return c.result

    def maxwell_boltzmann_test(self, init_position, num_trials=10):
        """
        Given a single position, it goes through a number of random velocity
        vectors, and returns the percentage of atoms that reach the fiber.

        """
        num_successes = 0.0
        for x in range(int(num_trials)):
            velocity = np.random.normal(MAX_BOLT_MU, MAX_BOLT_SIGMA, 3).tolist()
            cesium = Cesium(init_position, velocity)
            cesium.run()
            # result = volume_single.CesiumAtom().enters_fiber(init_position, velocity)

            if cesium.result is True:
                num_successes += 1

        # Return percentage of atoms that made it into the fiber at this position.
        return (num_successes / num_trials) * 100

    def y_range_test(self):
        """
        Test to see if atoms loaded in the MOT with y ranges 6.5 mm to 7.5 mm
        enters the fiber.
        z is constant at 0.

        """
        delta_y = 0.0001
        delta_x = 0.00003
        y_posit = .0068
        x_posit = 0.00033
        results = {}
        # while y_posit < .0076:
        #     x_posit = 0.0
        while x_posit < .0006:
            init_pos = [x_posit, y_posit, 0.0]
            print x_posit, y_posit, self.maxwell_boltzmann_test(init_pos, 2000)
            x_posit += delta_x
        # y_posit += delta_y

    def test_mine_taehyun(self, p_vector, v_vector):
        this = self.single_atom_test(p_vector, v_vector)
        taehyun = volume_single.CesiumAtom().enters_fiber(p_vector, v_vector)
        return this, taehyun


class Colormap:
    def create_percentage_matrix(self, filename='maxwell_boltzmann_test'):
        """
        Read the csv of the position, velocity, percentage
        lines and create a matrix of only the percentages or values.
        The rows are the y values, and the columns are the x values.
        A row by column matrix is returned.

        """
        color_map_matrix = []
        with open(filename, 'r') as file:
            nums = file.readline().strip().split(', ')
            current_y = float(nums[1])
            index = 0
            color_map_matrix.append([float(nums[-1])])
            for line in file:
                nums = line.strip().split(', ')
                percent = float(nums[-1])
                y = float(nums[1])
                if current_y != float(nums[1]):
                    index += 1
                    current_y = y
                    color_map_matrix.append([percent])
                else:
                    color_map_matrix[index].append(percent)
        return color_map_matrix

# Tests().y_range_test()
print Tests().test_mine_taehyun([0.00057, 0.0074, 0.0], [-0.015747562062037146, -0.022824595387545336, -0.00018671093214203338])
# print Tests().test_mine_taehyun([0.0005700000000000001, 0.007400000000000001, 0.0], [-0.015747562062037146, -0.022824595387545336, -0.00018671093214203338])
# print Tests().test_mine_taehyun([0.00048, 0.0072, 0], [-0.012938235379657806, -0.003702158386303568, 0.0003317082568845725])


# PLOTTING
# fig = pyplot.figure()
# ax = pyplot.axes(projection='3d')
# ax.plot3D(position_trace[0], position_trace[2], position_trace[1], 'red')
# pyplot.show()

# print('t = {0:.5f} s'.format(time))
# print('x = {0:.5f} m'.format(c.position[0]))
# print('y = {0:.5f} m'.format(c.position[1]))
# print('z = {0:.5f} m'.format(c.position[2]))
# print('vx = {0:.5f} m/s'.format(c.velocity[0]))
# print('vy = {0:.5f} m/s'.format(c.velocity[1]))
# print('vz = {0:.5f} m/s'.format(c.velocity[2]))
