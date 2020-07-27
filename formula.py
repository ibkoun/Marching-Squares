import math
import random
import numpy as np
    

# Line equation slope form : y = mx + b
# Line equation standard form : Ax + By = C where A = m, B = -1 and C = -b (y = mx + b => mx - y = - b)
class Segment:
    def __init__(self, start, end):
        self.start = start
        self.end = end
        self.vector = self.end - self.start
        dx = self.end[0] - self.start[0]
        dy = self.end[1] - self.start[1]
        if math.isclose(dx, 0):
            self.A = 1
            self.B = 0
            self.C = start[0]
        else:
            self.A = dy / dx
            self.B = -1
            self.C = self.A * self.start[0] + self.B * self.start[1]
        assert math.isclose(self.A * self.start[0] + self.B * self.start[1] - self.C, 0)

    # Find x knowing y.
    def x(self, y):
        dy = self.end[1] - self.start[1]
        if math.isclose(dy, 0):
            return random.random() * (self.end[0] - self.start[0]) + self.start[0]
        return (self.C - self.B * y) / self.A

    # Find y knowing x.
    def y(self, x):
        dx = self.end[0] - self.start[0]
        if math.isclose(dx, 0):
            return random.random() * (self.end[1] - self.start[1]) + self.start[1]
        return (self.C - self.A * x) / self.B

    # point = start + (end - start)t
    def t(self, point):
        x, y = point
        dx = self.end[0] - self.start[0]
        dy = self.end[1] - self.start[1]
        tx = (x - self.start[0]) / dx if not math.isclose(dx, 0) else 0
        ty = (y - self.start[1]) / dy if not math.isclose(dy, 0) else 0
        return tx, ty

    def squared_distance_from_point(self, point):
        x, y = point
        return math.pow(self.A * x + self.B * y - self.C, 2) / (math.pow(self.A, 2) + math.pow(self.B, 2))

    def distance_from_point(self, point):
        x, y = point
        return abs(self.A * x + self.B * y - self.C) / math.sqrt(math.pow(self.A, 2) + math.pow(self.B, 2))

    # Check if the segment intersects a point.
    def intersects_point(self, point):
        tx, ty = self.t(point)
        dx = self.end[0] - self.start[0]
        dy = self.end[1] - self.start[1]
        if math.isclose(dy, 0):  # Horizontal segment.
            return math.isclose(self.start[1], point[1]) and 0 <= tx <= 1
        elif math.isclose(dx, 0):  # Vertical segment.
            return math.isclose(self.start[0], point[0]) and 0 <= ty <= 1
        elif math.isclose(tx, ty):
            return 1 >= tx >= 0 <= ty <= 1
        return False

    # Check if the segment intersects with another.
    def intersects_segment(self, segment):
        # Calculate the angle between the two vectors.
        angle = angle_between(self.vector, segment.vector)

        # Check if the two segments are collinear.
        if math.isclose(angle, 0) or math.isclose(angle, 180):
            # Check whether one of the endpoints in a segment belongs to the other.
            return self.intersects_point(segment.start) or self.intersects_point(segment.end)\
                or segment.intersects_point(self.start) or segment.intersects_point(self.end)\

        # Check if the two segments have an intersection point.
        x = (self.B * segment.C - segment.B * self.C) / (self.B * segment.A - segment.B * self.A)
        y = (self.A * segment.C - segment.A * self.C) / (self.A * segment.B - segment.A * self.B)
        return self.intersects_point([x, y]) and segment.intersects_point([x, y])

    def intersection_point(self, segment):
        x = (self.B * segment.C - segment.B * self.C) / (self.B * segment.A - segment.B * self.A)
        y = (self.A * segment.C - segment.A * self.C) / (self.A * segment.B - segment.A * self.B)
        return np.array([x, y])


def resize_vector(vector, magnitude):
    return (magnitude / math.sqrt(np.dot(vector, vector))) * vector


def rotate_vector(vector, angle):
    radians = math.radians(angle)
    rotation = np.array([[math.cos(radians), -math.sin(radians)], [math.sin(radians), math.cos(radians)]])
    return np.dot(rotation, vector)


# Project vector u on vector v.
def project_vector(u, v):
    return (np.dot(u, v) / np.dot(v, v)) * v


# Calculate the angle between vector u and vector v.
def angle_between(u, v, signed=False):
    mu = math.sqrt(np.dot(u, u))
    mv = math.sqrt(np.dot(v, v))

    # Calculate the signed angle necessary to rotate u to v.
    if signed:
        ua = math.degrees(math.acos(u[0] / mu))
        va = math.degrees(math.acos(v[0] / mv))
        return va - ua
    return math.degrees(math.acos(round(np.dot(u, v) / (mu * mv), 2)))


# TODO: Case when math.sqrt(math.pow(b, 2) - 4 * a * c) < 0
def quadratic(a, b, c):
    discriminant = math.sqrt(math.pow(b, 2) - 4 * a * c)
    v1 = (-b - discriminant) / (2 * a)
    v2 = (-b + discriminant) / (2 * a)
    return v1, v2


# Return the coordinates of a point on the outline of a circle as a function of two points along a segment.
def point_on_circumference(radius, center, start, end):
    ax = math.pow(end[0] - start[0], 2)
    ay = math.pow(end[1] - start[1], 2)
    a = ax + ay

    bx = -2 * (end[0] - start[0]) * (center[0] - start[0])
    by = -2 * (end[1] - start[1]) * (center[1] - start[1])
    b = bx + by

    cx = math.pow((center[0] - start[0]), 2)
    cy = math.pow((center[1] - start[1]), 2)
    c = cx + cy - math.pow(radius, 2)
    t1, t2 = quadratic(a, b, c)

    delta = end - start
    p1 = start + t1 * delta
    p2 = start + t2 * delta
    delta1 = end - p1
    delta2 = end - p2
    distance1 = np.dot(delta1, delta1)
    distance2 = np.dot(delta2, delta2)
    return p1 if distance1 < distance2 else p2
