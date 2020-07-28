import numpy as np
import math
import formula
from random import random
from abc import ABC, ABCMeta, abstractmethod


class Graphic2D(metaclass=ABCMeta):
    @classmethod
    def __subclasscheck__(cls, subclass):
        return (hasattr(subclass, 'draw') and
                callable(subclass.draw) and
                hasattr(subclass, 'redraw') and
                callable(subclass.draw) or
                NotImplemented)

    def draw(self, canvas, **kwargs):
        raise NotImplementedError

    def redraw(self, canvas, **kwargs):
        raise NotImplementedError


class Collider2D(metaclass=ABCMeta):
    @classmethod
    def __subclasscheck__(cls, subclass):
        return (hasattr(subclass, 'contains_point') and
                callable(subclass.contains_point) and
                hasattr(subclass, 'collides_circle') and
                callable(subclass.collides_circle) or
                hasattr(subclass, 'overlaps_circle') and
                callable(subclass.overlaps_circle) and
                hasattr(subclass, 'confines_circle') and
                callable(subclass.confines_circle) and
                NotImplemented)

    # Check if a coordinate is within the boundaries.
    def contains_point(self, point):
        raise NotImplementedError

    # Check for a collision with a circle.
    def collides_circle(self, circle):
        raise NotImplementedError

    # Check if a circle crosses the boundaries.
    def overlaps_circle(self, circle):
        raise NotImplementedError

    # Check if a circle is within the boundaries.
    def confines_circle(self, circle):
        raise NotImplementedError


class Shape2D(ABC, Graphic2D, Collider2D):
    def __init__(self, x, y, rotation=float(0)):
        self._center = np.array([x, y])
        self._rotation = rotation
        self._item = None

    def get_center(self):
        return np.copy(self._center)

    def set_center(self, x, y):
        self._center = np.array([x, y])

    def get_rotation(self):
        return self._rotation

    def set_rotation(self, rotation):
        self._rotation = rotation

    def get_item(self):
        return self._item

    def set_item(self, item):
        self._item = item

    def squared_distance_from_point(self, point):
        vector = point - self._center
        return np.dot(vector, vector)

    def distance_from_point(self, point):
        vector = point - self._center
        return math.sqrt(np.dot(vector, vector))

    @abstractmethod
    def randomize_circle_coord(self, circle):
        raise NotImplementedError


class Circle(Shape2D):
    def __init__(self, x, y, radius):
        super(Circle, self).__init__(x, y)
        self._radius = radius

    def get_radius(self):
        return self._radius

    def set_radius(self, radius):
        self._radius = radius

    def randomize_radius(self, min_radius, max_radius):
        self._radius = random() * (min_radius - max_radius) + min_radius

    def draw(self, canvas, **kwargs):
        if not self._item:
            x, y = self._center
            x1 = x - self._radius
            x2 = x + self._radius
            y1 = y - self._radius
            y2 = y + self._radius
            self._item = canvas.create_oval(x1, y1, x2, y2, **kwargs)

    def redraw(self, canvas, **kwargs):
        if self._item:
            x, y = self._center
            x1 = x - self._radius
            x2 = x + self._radius
            y1 = y - self._radius
            y2 = y + self._radius
            canvas.coords(self._item, x1, y1, x2, y2)
            canvas.itemconfig(self._item, **kwargs)

    def randomize_circle_coord(self, circle):
        # Generate x randomly inside the area.
        radius = circle.get_radius()
        x_range = 2 * (self._radius - radius)
        min_x = self._center[0] - x_range / 2
        x = random() * x_range + min_x

        # Generate y randomly inside the area.
        y_range = 2 * math.sqrt(math.pow(self._radius - radius, 2) - math.pow(x - self._center[0], 2))
        min_y = self._center[1] - y_range / 2
        y = random() * y_range + min_y
        circle.set_center(x, y)
        # return np.array([x, y])

    def contains_point(self, point):
        vector = self._center - point
        distance = np.dot(vector, vector)
        threshold = math.pow(self._radius, 2)
        return distance < threshold or math.isclose(distance, threshold)

    def confine_circle_coord(self, circle, start, end):
        return formula.point_on_circumference(self._radius - circle.get_radius(), self._center, start, end)

    def collides_circle(self, circle):
        vector = self._center - circle.get_center()
        distance = np.dot(vector, vector)
        threshold = math.pow(self._radius + circle.get_radius(), 2)
        return distance < threshold or math.isclose(distance, threshold)

    def overlaps_circle(self, circle):
        vector = self._center - circle.get_center()
        distance = np.dot(vector, vector)
        threshold = math.pow(self._radius + circle.get_radius(), 2)
        return distance < threshold and not math.isclose(distance, threshold)

    def confines_circle(self, circle):
        vector = circle.get_center() - self._center
        distance = np.dot(vector, vector)
        threshold = math.pow(self._radius - circle.get_radius(), 2)
        return distance < threshold or math.isclose(distance, threshold)

    def distance_from_circle(self, circle):
        vector = circle.get_center() - self._center
        distance = math.sqrt(np.dot(vector, vector))
        return distance - self._radius - circle.get_radius()


class Rectangle(Shape2D):
    def __init__(self, x, y, width, height):
        super(Rectangle, self).__init__(x, y)
        self._width = width
        self._height = height

    def get_center(self):
        return self._center

    def get_width(self):
        return self._width

    def set_width(self, width):
        self._width = width

    def get_height(self):
        return self._height

    def set_height(self, height):
        self._height = height

    def draw(self, canvas, **kwargs):
        if not self._item:
            x = self._center[0]
            y = self._center[1]
            width = self._width / 2
            height = self._height / 2
            x1 = x - width
            x2 = x + width
            y1 = y - height
            y2 = y + height
            self._item = canvas.create_rectangle(x1, y1, x2, y2, **kwargs)

    def redraw(self, canvas, **kwargs):
        if self._item:
            canvas.itemconfig(self._item, **kwargs)

    def randomize_circle_coord(self, circle):
        # Generate x randomly inside the area.
        center = self.get_center()
        x_range = self.get_width() - 2 * circle.get_radius()
        min_x = center[0] - x_range / 2
        x = random() * x_range + min_x

        # Generate y randomly inside the area.
        y_range = self.get_height() - 2 * circle.get_radius()
        min_y = center[1] - y_range / 2
        y = random() * y_range + min_y
        circle.set_center(x, y)

    # TODO: Fix this method.
    def confine_circle_coord(self, circle, start, end):
        center = self.get_center()
        radius = circle.get_radius()
        width = self._width
        height = self._height
        point = circle.get_center()

        # canvas.create_oval(point[0] - radius, point[1] - radius, point[0] + radius, point[1] + radius, fill="blue")

        # The boundary of the current quadrant.
        x1 = center[0] - width / 2
        x2 = center[0] + width / 2
        y1 = center[1] - height / 2
        y2 = center[1] + height / 2

        # The four corners of the current quadrant.
        north_west = np.array([x1, y1])
        north_east = np.array([x2, y1])
        south_west = np.array([x1, y2])
        south_east = np.array([x2, y2])

        # The borders of the current quadrant.
        north_border = formula.Segment(north_west, north_east)
        south_border = formula.Segment(south_west, south_east)
        west_border = formula.Segment(north_west, south_west)
        east_border = formula.Segment(north_east, south_east)

        segment = formula.Segment(start, end)

        if not self.contains_point(point):
            if segment.intersects_segment(north_border):
                point = segment.intersection_point(north_border)
                y = point[1] + radius
                x = segment.x(y)
                point[0], point[1] = x, y
            elif segment.intersects_segment(south_border):
                point = segment.intersection_point(south_border)
                y = point[1] - radius
                x = segment.x(y)
                point[0], point[1] = x, y
            elif segment.intersects_segment(west_border):
                point = segment.intersection_point(west_border)
                x = point[0] + radius
                y = segment.y(x)
                point[0], point[1] = x, y
            elif segment.intersects_segment(east_border):
                point = segment.intersection_point(east_border)
                x = point[0] - radius
                y = segment.y(x)
                point[0], point[1] = x, y
        # canvas.create_oval(point[0] - radius, point[1] - radius, point[0] + radius, point[1] + radius, fill="red")
        if np.sign(segment.vector[0]) < 0:
            distance_from_west = west_border.distance_from_point(point)
            if np.sign(segment.vector[1]) < 0:
                distance_from_north = north_border.distance_from_point(point)
                if distance_from_west < distance_from_north:
                    x = point[0] - distance_from_west + radius
                    y = segment.y(x)
                    point[0], point[1] = x, y
                else:
                    y = point[1] - distance_from_north + radius
                    x = segment.x(y)
                    point[0], point[1] = x, y
            elif np.sign(segment.vector[1]) > 0:
                distance_from_south = south_border.distance_from_point(point)
                if distance_from_west < distance_from_south:
                    x = point[0] - distance_from_west + radius
                    y = segment.y(x)
                    point[0], point[1] = x, y
                else:
                    y = point[1] + distance_from_south - radius
                    x = segment.x(y)
                    point[0], point[1] = x, y
        elif np.sign(segment.vector[0]) > 0:
            distance_from_east = east_border.distance_from_point(point)
            if np.sign(segment.vector[1]) < 0:
                distance_from_north = north_border.distance_from_point(point)
                if distance_from_east < distance_from_north:
                    x = point[0] + distance_from_east - radius
                    y = segment.y(x)
                    point[0], point[1] = x, y
                else:
                    y = point[1] - distance_from_north + radius
                    x = segment.x(y)
                    point[0], point[1] = x, y
            elif np.sign(segment.vector[1]) > 0:
                distance_from_south = south_border.distance_from_point(point)
                if distance_from_east < distance_from_south:
                    x = point[0] + distance_from_east - radius
                    y = segment.y(x)
                    point[0], point[1] = x, y
                else:
                    y = point[1] + distance_from_south - radius
                    x = segment.x(y)
                    point[0], point[1] = x, y


        # canvas.create_oval(point[0] - radius, point[1] - radius, point[0] + radius, point[1] + radius, fill="green")

        """distance_from_north = north_border.distance_from_point(point)
            if distance_from_north < radius and not math.isclose(distance_from_north, radius):
                y = point[1] - distance_from_north + radius
                x = segment.x(y)
                point[0], point[1] = x, y

            distance_from_south = south_border.distance_from_point(point)
            if distance_from_south < radius and not math.isclose(distance_from_south, radius):
                y = point[1] + distance_from_south - radius
                x = segment.x(y)
                point[0], point[1] = x, y

            distance_from_west = west_border.distance_from_point(point)
            if distance_from_west < radius and not math.isclose(distance_from_west, radius):
                x = point[0] - distance_from_west + radius
                y = segment.y(x)
                point[0], point[1] = x, y

            distance_from_east = east_border.distance_from_point(point)
            if distance_from_east < radius and not math.isclose(distance_from_east, radius):
                x = point[0] + distance_from_east - radius
                y = segment.y(x)
                point[0], point[1] = x, y
        else:
            distance_from_north = north_border.distance_from_point(point)
            if point[1] < center[1] - height / 2 + radius and not math.isclose(point[1], center[1] - height / 2 + radius):
                y = point[1] + distance_from_north + radius
                x = segment.x(y)
                point[0], point[1] = x, y

            distance_from_south = south_border.distance_from_point(point)
            if point[1] > center[1] + height / 2 - radius and not math.isclose(point[1], center[1] + height / 2 - radius):
                y = point[1] - distance_from_south - radius
                x = segment.x(y)
                point[0], point[1] = x, y

            distance_from_west = west_border.distance_from_point(point)
            if point[0] < center[0] - width / 2 + radius and not math.isclose(point[0], center[0] - width / 2 + radius):
                x = point[0] + distance_from_west + radius
                y = segment.y(x)
                point[0], point[1] = x, y

            distance_from_east = east_border.distance_from_point(point)
            if point[0] > center[0] + width / 2 - radius and not math.isclose(point[0], center[0] + width / 2 - radius):
                x = point[0] - distance_from_east - radius
                y = segment.y(x)
                point[0], point[1] = x, y

            distance_from_north = north_border.distance_from_point(point)
            if distance_from_north < radius and not math.isclose(distance_from_north, radius):
                y = point[1] - distance_from_north + radius
                x = segment.x(y)
                point[0], point[1] = x, y

            distance_from_south = south_border.distance_from_point(point)
            if distance_from_south < radius and not math.isclose(distance_from_south, radius):
                y = point[1] + distance_from_south - radius
                x = segment.x(y)
                point[0], point[1] = x, y

            distance_from_west = west_border.distance_from_point(point)
            if distance_from_west < radius and not math.isclose(distance_from_west, radius):
                x = point[0] - distance_from_west + radius
                y = segment.y(x)
                point[0], point[1] = x, y

            distance_from_east = east_border.distance_from_point(point)
            if distance_from_east < radius and not math.isclose(distance_from_east, radius):
                x = point[0] + distance_from_east - radius
                y = segment.y(x)
                point[0], point[1] = x, y"""
        return point

    def contains_point(self, point):
        x, y = point
        cx, cy = self._center
        x1 = cx - self._width / 2
        x2 = cx + self._width / 2
        y1 = cy - self._height / 2
        y2 = cy + self._height / 2
        return x1 <= x <= x2 and y1 <= y <= y2

    def collides_circle(self, circle):
        cx = self._center[0]
        cy = self._center[1]
        x, y = circle.get_center()
        left = x + circle.get_radius() > cx - self._width / 2 or math.isclose(x + circle.get_radius(),
                                                                              cx - self._width / 2)
        right = x - circle.get_radius() < cx + self._width / 2 or math.isclose(x - circle.get_radius(),
                                                                               cx + self._width / 2)
        top = y + circle.get_radius() > cy - self._height / 2 or math.isclose(y + circle.get_radius(),
                                                                              cy - self._height / 2)
        bottom = y - circle.get_radius() < cy + self._height / 2 or math.isclose(y - circle.get_radius(),
                                                                                 cy + self._height / 2)
        return left and right and top and bottom

    def overlaps_circle(self, circle):
        cx = self._center[0]
        cy = self._center[1]
        x, y = circle.get_center()
        left = x + circle.get_radius() > cx - self._width / 2 and not math.isclose(x + circle.get_radius(),
                                                                                   cx - self._width / 2)
        right = x - circle.get_radius() < cx + self._width / 2 and not math.isclose(x - circle.get_radius(),
                                                                                    cx + self._width / 2)
        top = y + circle.get_radius() > cy - self._height / 2 and not math.isclose(y + circle.get_radius(),
                                                                                   cy - self._height / 2)
        bottom = y - circle.get_radius() < cy + self._height / 2 and not math.isclose(y - circle.get_radius(),
                                                                                      cy + self._height / 2)
        return left and right and top and bottom

    def confines_circle(self, circle):
        cx = self._center[0]
        cy = self._center[1]
        x, y = circle.get_center()
        left = x - circle.get_radius() > cx - self._width / 2 \
               or math.isclose(x - circle.get_radius(), cx - self._width / 2)
        right = x + circle.get_radius() < cx + self._width / 2 \
                or math.isclose(x + circle.get_radius(), cx + self._width / 2)
        top = y - circle.get_radius() > cy - self._height / 2 \
              or math.isclose(y - circle.get_radius(), cy - self._height / 2)
        bottom = y + circle.get_radius() < cy + self._height / 2 \
                 or math.isclose(y + circle.get_radius(), cy + self._height / 2)
        return left and right and top and bottom


class Segment(Shape2D):
    def __init__(self, start, end):
        self.start = start
        self.end = end
        self.vector = self.end - self.start
        center = self.start + self.vector / 2
        rotation = math.acos(self.vector[0] / self.magnitude())
        super(Segment, self).__init__(center[0], center[1], rotation)
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

    def draw(self, canvas, **kwargs):
        if not self._item:
            x1, y1 = self.start[0], self.start[1]
            x2, y2 = self.end[0], self.end[1]
            self._item = canvas.create_line(x1, y1, x2, y2, **kwargs)

    def redraw(self, canvas, **kwargs):
        if self._item:
            x1, y1 = self.start[0], self.start[1]
            x2, y2 = self.end[0], self.end[1]
            canvas.itemconfig(self._item, x1, y1, x2, y2, **kwargs)

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

    def squared_magnitude(self):
        return np.dot(self.vector, self.vector)

    def magnitude(self):
        return math.sqrt(np.dot(self.vector, self.vector))

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
        angle = formula.angle_between(self.vector, segment.vector)

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

    def randomize_circle_coord(self, circle):
        pass
