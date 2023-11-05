#!/usr/bin/env python3
import math as math_lib
from mpl_toolkits import mplot3d
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys


"""
This file is mainly responsible for simulating particles that 
interact with eachother. Actual representation of this happens
elsewhere, we just take in current coordinates, previous vectors, and
adjust them for current time.

The more often you call the interact function for here, the greater the 
accuracy you get for the model. However, the slower it is to compute.
"""

"""
Start by defining properites of a particle. Do note, that interaction
between multiple particles is NOT handled by the move() function, rather
it only applies vectors to coordinates.

You shall first calculate the gravity vector between two particles from
their masses & time of 1 frame. Then break the gravity vector into
x, y, and z vectors, and add those to particle-bound x-, y-, and z vector.
"""
class particle:
    def __init__(self, mass, initial_x, initial_y, initial_z,
                 initial_vector_x, initial_vector_y, initial_vector_z):
        self.mass = mass
        self.coordinates = [initial_x, initial_y, initial_z]
        self.vectors = [initial_vector_x, 
                        initial_vector_y, 
                        initial_vector_z]

    def move(self) -> None:
        for i in range(0, 3, 1):
            tmp = self.coordinates[i]
            self.coordinates[i] += self.vectors[i]
            if (tmp != 0):
                if (self.coordinates[i] == self.vectors[i]):
                    print("BUG: %f + %f = %f" %
                          (tmp, self.coordinates[i], self.vectors[i]))
                    exit(0)
    
    def x(self) -> float:
        return self.coordinates[0]

    def y(self) -> float:
        return self.coordinates[1]

    def z(self) -> float:
        return self.coordinates[2]

    def xvect(self) -> float:
        return self.vectors[0]

    def yvect(self) -> float:
        return self.vectors[1]

    def zvect(self) -> float:
        return self.vectors[2]

"""
Physics

NOTE: We are dealing with 3-dimensional space, so instead of the
      traditional trigonometrical funcs, we're using a bit altered variants

    A random ascii diagram for how to calculate distance to X/Y/Z axis, which
    is base of all calculus performed here
     |      /
     |     / \
     |    /   \
     |   /     \
     |  /       \
     |-/---------*
     |/__________|_
"""

class physics:
    # Constants
    #
    gravitational_constant = 6.67 * (10 ** -11)

    # Calculate the distance from point A(x, y, z) to point B(x, y, z)
    # this is just square root of (a^2 + b^2 + c^2) since pythagorean theorem
    # goes as:
    #   c^2 = a^2 + b^2 + ... n^2
    #
    # We could go substracting x1 - x2, y1 - y2, z1 - z2 but that's lame and 
    # doesn't give us neat single number distance, rather just 3 different
    # distances on separate axis, whilst we need absolute distance for
    # physics stuff for example.
    #
    def distance(size_x, size_y, size_z) -> float:
        return math_lib.sqrt((size_x ** 2) + (size_y ** 2) + (size_z ** 2))

    # Caclculate how quickly particles 1 and 2 start to accelerate towards
    # eachother at their current distance from eachother.
    #
    # Return [particle1 acceleration, particle2 acceleration]
    def gravity_accelerations(p1, p2, distance):
        # Calculate mass-distance relation
        mdr = p1.mass / (distance ** 2)
        field1_strength = physics.gravitational_constant * mdr
        mdr2 = p2.mass / (distance ** 2)
        field2_strength = physics.gravitational_constant * mdr2
        return [p1.mass * field2_strength, p2.mass * field1_strength]
    
    # Convert acceleration to average velocity of 1 frame
    # NOTE:
    #   this function doesn't care about which direction the acceleration is
    #   applied to. I ***strongly*** recommend FIRST splitting gravity to
    #   X, Y, Z vectors, and then operating with each of those independently
    #
    def acceleration_to_average_velocity(acceleration, framesize=1):
        final_velocity = acceleration + (acceleration * framesize)
        return final_velocity / 2.0

    # Break a some-directional gravity vector into x,y,z vectors
    # NOTE: Below is horribly utterly oversimplified and wrong
    def gvec_to_directionals(xsize, ysize, zsize, relation):
        gvec_x = xsize * relation
        gvec_y = ysize * relation
        gvec_z = zsize * relation
        return [gvec_x, gvec_y, gvec_z]

    # Calculate a convenient X, Y, Z vector group for gravitational pull
    # between particles 1 and 2
    #
    def gravity_vectors(p1, p2, framesize=1):
        xsize = p1.x() - p2.x()
        ysize = p1.y() - p2.y()
        zsize = p1.z() - p2.z()
        sizes = [xsize, ysize, zsize]
        nc = 0
        for i in range(len(sizes)):
            if (sizes[i] < 0):
                sizes[i] *= -1
            elif (sizes[i] == 0):
                nc += 1
        if (nc < 2):
            distance = physics.distance(xsize, ysize, zsize)
        else:
            distance = xsize + ysize + zsize
        print("Distance: %f" % distance)
        if (distance == 0):
            print("Collision at %f / %f / %f" % (p1.x(), p1.y(), p1.z()))
            exit(0)

        p1_gvec, p2_gvec = physics.gravity_accelerations(p1, p2, distance)

        # Since angles between p1, p2 haven't changed, the relation between
        # gravity triangle sides, and distance triangle sizes is linear. We
        # know all the angles are same/identical.
        #
        # We know the relation between hypotenuse with gravity vector
        # versus the distance, so we can apply that to all the remaning sides
        # of the triangles.
        #
        # p1_gvec_d = physics.gvec_to_directionals(p1, p2, distance, p1_gvec)
        relation = p1_gvec / distance
        p1_gvec_d = physics.gvec_to_directionals(xsize, ysize, zsize, relation)

        for i in range(0, 3, 1):
            p1_gvec_d[i] = physics.acceleration_to_average_velocity(
                                                            p1_gvec_d[i], 
                                                            framesize)

        # p2_relation = p2_gvec / distance
        # p2_gvec_d = gvec_to_directionals(xsize, ysize, zsize, relations)
        # for i in range(0, 3, 1):
        #   p2_gvec_d[i] = ...

        # We now have vectors, finally adjust their directions
        for i in range(0, 3, 1):
           if (p1.coordinates[i] > p2.coordinates[i]):
               if (p1_gvec_d[i] > 0):
                   p1_gvec_d[i] *= -1
           else:
               if (p1_gvec_d[i] < 0):
                   p1_gvec_d[i] *= -1
        return p1_gvec_d

    def infodump(prev, current, vector):
        if (prev != 0):
            if (current == vector):
                print("BUG: %f + %f = %f" % (
                    prev, vector, current))
                exit(0)
        print("%f -> %f - vector: %f" % (
            prev,
            current,
            vector))

    # Run a single frame of simulation
    def run_frame(*particles, framesize=1):
        if (len(particles) < 1):
            return
        elif (len(particles) == 1):
            particles[0].move()
            return
        for c in range(len(particles)):
            for o in range(0, len(particles), 1):
                if (c == o):
                    continue
                gvecs = physics.gravity_vectors(particles[c], particles[o],
                                                framesize)
                print("gravity vectors: ")
                print(gvecs)
                for i in range(0, 3, 1):
                    particles[c].vectors[i] += gvecs[i]
        print("-" * 80)
        for i in range(len(particles)):
            prevx = particles[i].x()
            prevy = particles[i].y()
            prevz = particles[i].z()
            particles[i].move()
            print("particle %d moved" % i)
            physics.infodump(prevx, particles[i].x(), particles[i].xvect())
            physics.infodump(prevy, particles[i].y(), particles[i].yvect())
            physics.infodump(prevz, particles[i].z(), particles[i].zvect())


class simulation:
    def __init__(self, *particles, ticks=5000, tick_size=1):
        self.particles = particles
        self.xlines = []
        self.ylines = []
        self.zlines = []
        self.ticks = ticks
        self.framesize = tick_size

    def run(self):
        for p in self.particles:
            self.xlines.append([])
            self.ylines.append([])
            self.zlines.append([])
        for i in range(0, self.ticks, 1):
            print("*" * 80)
            print("Frame: %d" % i)
            print("*" * 80)
            physics.run_frame(*self.particles, framesize=self.framesize)
            for p in range(len(self.particles)):
                self.xlines[p].append(self.particles[p].x() / 1000)
                self.ylines[p].append(self.particles[p].y() / 1000)
                self.zlines[p].append(self.particles[p].z() / 1000)

    def draw_graph(self):
        try:
            plt.ion()
            fig = plt.figure()
            for i in range(len(self.xlines[0])):
                ax = plt.axes(projection='3d')
                for p in range(len(self.particles)):
                    ax.plot3D(self.xlines[p][:i], self.ylines[p][:i], self.zlines[p][:i], colour(p))
                plt.draw()
                plt.pause(0.005)
                plt.clf()
        except Exception as E:
            print(E)

def colour(x):
    colours = ["red", "blue", "gray", "yellow", "green"]
    return colours[x % len(colours)]

tick_sz=1
tick_cnt=9000
if (len(sys.argv) == 2):
    tick_cnt = int(sys.argv[1])

# particle(mass, x, y, z, xvector, yvector, zvector)
p1 = particle(5.3 * (10 ** 13), 100124124, -934124996, 999966435, 900000, 90124, 525667)
p2 = particle(5.3 * (10 ** 13), 942141242, 112414720, 535363430, -4387860, 0, 0)
p3 = particle(5.3 * (10 ** 18), 500000, 100000, 100000, -50000, 5000, 5000)

sim = simulation(p1, p2, p3, ticks=tick_cnt, tick_size=tick_sz)
sim.run()
sim.draw_graph()

