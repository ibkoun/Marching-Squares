from graphic import Rectangle, Graphic2D
import math
import formula
import numpy as np


class Node:
    def __init__(self, index=None):
        self._index = index
        self._contents = []
        self._leaves = []
    
    def contents(self):
        return self._contents
    
    def leaves(self):
        return self._leaves


class Quadrant(Node, Rectangle):
    def __init__(self, x, y, width, height, index=None):
        Node.__init__(self, index=index)
        Rectangle.__init__(self, x, y, width, height)
        self.canvas = None

    def inside(self, event):
        x, y = event.x, event.y
        if self.contains_point(np.array([x, y])):
            print(self._center)
            self.canvas.create_oval(self._center[0] - 5,
                               self._center[1] - 5,
                               self._center[0] + 5,
                               self._center[1] + 5,
                               outline="blue")

    def draw(self, canvas, fill="", outline="black"):
        super(Quadrant, self).draw(canvas, fill=fill, outline=outline)
        # for i in range(len(self._contents)):
            # self._contents[i].draw(canvas, fill=fill, outline=outline)

    def redraw(self, canvas, fill="", outline="black"):
        super(Quadrant, self).redraw(canvas, fill=fill, outline=outline)
        # for i in range(len(self._contents)):
            # self._contents[i].redraw(canvas, fill=fill, outline=outline)
        
    def partition(self):
        if len(self._leaves) == 0:
            width = self._width / 2
            height = self._height / 2
            for i in range(4):
                coord = self.sub_quadrant_coord(i)
                x, y = coord
                self._leaves.append(Quadrant(x, y, width, height))
                
    def sub_quadrant_coord(self, index):
        assert 0 <= index < 4, "Index is out of bounds."
        x, y = self._center
        width = self._width / 4
        height = self._height / 4
        if index == 0:
            x += width
            y -= width
        elif index == 1:
            x -= width
            y -= height
        elif index == 2:
            x -= width
            y += height
        elif index == 3:
            x += width
            y += height
        return np.array([x, y])
    
    def circle_overlap(self, circle):
        quadrants = []
        for quadrant in self._leaves:
            if quadrant.collides_circle(circle):
                quadrants.append(quadrant)
        return quadrants


class Quadtree(Graphic2D):
    def __init__(self, root=None):
        self._contents = []
        self._root = root
        self.quadtree_lookups = 0
        self.quadtree_comparisons = 0
        self.linear_comparisons = 0
        self.collisions = 0

    def contents(self):
        return self._contents
        
    def get_root(self):
        return self._root
    
    def set_root(self, root):
        self._root = root

    def draw(self, canvas, fill="", outline="black"):
        def depth_first_search(quadrant):
            quadrant.draw(canvas, fill=fill, outline=outline)
            leaves = quadrant.leaves()
            for i in range(len(leaves)):
                depth_first_search(leaves[i])
        depth_first_search(self._root)

    def redraw(self, canvas, fill="", outline="black"):
        def depth_first_search(quadrant):
            quadrant.redraw(canvas, fill=fill, outline=outline)
            leaves = quadrant.leaves()
            for i in range(len(leaves)):
                depth_first_search(leaves[i])
        depth_first_search(self._root)

    def result(self, n):
        count = len(self._contents)
        for i in range(count):
            for j in range(i + 1, count):
                if self._contents[i].collides_circle(self._contents[j]):
                    self.collisions += 1
        print("RESULT")
        print("Collisions: {} ({}%)".format(self.collisions, self.collisions / count * 100))
        print("Count: {}/{} ({}%)".format(count, n, count / n * 100))

    def linear_search(self, circle, overlap=False):
        if not overlap:
            for i in range(len(self._contents)):
                self.linear_comparisons += 1
                if circle.collides_circle(self._contents[i]):
                    return False
        return True

    def linear_search_result(self):
        print("LINEAR SEARCH")
        print("Comparisons: {}".format(self.linear_comparisons))

    def overlapped_by_circle(self, circle, leaves_only=False):
        queue = [self._root]
        quadrants = []
        while len(queue) > 0:
            quadrant = queue.pop(0)
            if len(quadrant.leaves()) > 0:
                if not leaves_only:
                    if quadrant.collides_circle(circle):
                        quadrants.append(quadrant)
                candidates = quadrant.leaves()
                queue += candidates
            else:
                if quadrant.collides_circle(circle):
                    quadrants.append(quadrant)
        return quadrants

    def quadtree_search(self, circle, overlap=False):
        if not overlap:
            queue = [self._root]
            quadrants = []
            while len(queue) > 0:
                self.quadtree_lookups += 1
                quadrant = queue.pop(0)
                if len(quadrant.contents()) == 0:
                    if len(quadrant.leaves()) == 0:
                        quadrants.append(quadrant)
                    else:
                        queue += quadrant.circle_overlap(circle)
                else:
                    circles = quadrant.contents()
                    for i in range(len(circles)):
                        self.quadtree_comparisons += 1
                        if circle.overlaps_circle(circles[i]):
                            quadrants.clear()
                            return quadrants
                    if len(quadrant.leaves()) == 0:
                        width = quadrant.get_width() / 2
                        height = quadrant.get_height() / 2
                        radius = 2 * circle.get_radius()
                        if (width > radius or math.isclose(width, radius))\
                                and (height > radius or math.isclose(height, radius)):
                            quadrant.partition()
                        if len(quadrant.leaves()) == 0:
                            quadrants.append(quadrant)
                        else:
                            queue.append(quadrant)
                    else:
                        queue += quadrant.circle_overlap(circle)
            return quadrants
        return [self._root]

    def rectangle_overlap(self, start, end, margin, canvas):
        queue = [self._root]
        quadrants = []

        # Vectors used to find the four corners of the rectangle.
        v1 = formula.resize_vector(end - start, margin)
        v2 = -v1
        v3 = formula.rotate_vector(v1, 90)
        v4 = -v3

        # The four corners of the rectangle.
        p1 = start + v2 + v3
        p2 = start + v2 + v4
        p3 = end + v1 + v3
        p4 = end + v1 + v4

        # The borders of the rectangle
        border_12 = formula.Segment(p1, p2)
        border_13 = formula.Segment(p1, p3)
        border_24 = formula.Segment(p2, p4)
        border_34 = formula.Segment(p3, p4)

        # Visual representation of the rectangle.
        # canvas.create_line(p1[0], p1[1], p2[0], p2[1], fill="blue", width=2)
        # canvas.create_line(p3[0], p3[1], p4[0], p4[1], fill="blue", width=2)
        # canvas.create_line(p1[0], p1[1], p3[0], p3[1], fill="blue", width=2)
        # canvas.create_line(p2[0], p2[1], p4[0], p4[1], fill="blue", width=2)

        # Iterate through a list of quadrants overlapped by the rectangle.
        while len(queue) > 0:
            quadrant = queue.pop(0)
            center = quadrant.get_center()
            width = quadrant.get_width()
            height = quadrant.get_height()
            centerline = formula.Segment(start, end)
            distance_from_centerline = centerline.distance_from_point(center)

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

            # Check if the quadrant is inside the rectangle.
            if distance_from_centerline < margin or math.isclose(distance_from_centerline, margin):
                quadrants.append(quadrant)
                queue += quadrant.leaves()

            # Check if one of the borders of the rectangle goes through the center of the current quadrant.
            elif border_12.intersects_point(center)\
                    or border_13.intersects_point(center)\
                    or border_24.intersects_point(center)\
                    or border_34.intersects_point(center):
                quadrants.append(quadrant)
                queue += quadrant.leaves()
            
            # Check if one of the corners of the rectangle is inside the current quadrant.
            elif quadrant.contains_point(p1):
                quadrants.append(quadrant)
                queue += quadrant.leaves()
            elif quadrant.contains_point(p2):
                quadrants.append(quadrant)
                queue += quadrant.leaves()
            elif quadrant.contains_point(p3):
                quadrants.append(quadrant)
                queue += quadrant.leaves()
            elif quadrant.contains_point(p4):
                quadrants.append(quadrant)
                queue += quadrant.leaves()
                
            # Check if one of the segments of the rectangle intersects a quadrant.
            elif west_border.intersects_segment(border_12)\
                    or east_border.intersects_segment(border_12)\
                    or north_border.intersects_segment(border_12)\
                    or south_border.intersects_segment(border_12):
                quadrants.append(quadrant)
                queue += quadrant.leaves()
            elif west_border.intersects_segment(border_13)\
                    or east_border.intersects_segment(border_13)\
                    or north_border.intersects_segment(border_13)\
                    or south_border.intersects_segment(border_13):
                quadrants.append(quadrant)
                queue += quadrant.leaves()
            elif west_border.intersects_segment(border_24)\
                    or east_border.intersects_segment(border_24)\
                    or north_border.intersects_segment(border_24)\
                    or south_border.intersects_segment(border_24):
                quadrants.append(quadrant)
                queue += quadrant.leaves()
            elif west_border.intersects_segment(border_34)\
                    or east_border.intersects_segment(border_34)\
                    or north_border.intersects_segment(border_34)\
                    or south_border.intersects_segment(border_34):
                quadrants.append(quadrant)
                queue += quadrant.leaves()
        return quadrants
                    
    def quadtree_search_result(self):
        print("QUADTREE SEARCH")
        print("Lookups: {}".format(self.quadtree_lookups))
        print("Comparisons: {}".format(self.quadtree_comparisons))
        print("Total: {}".format(self.quadtree_lookups + self.quadtree_comparisons))

    def quadrants_count(self):
        count = 0
        queue = [self._root]
        while len(queue) > 0:
            quadrant = queue.pop(0)
            count += 1
            if len(quadrant.leaves()) > 0:
                queue += quadrant.leaves()
        return count
