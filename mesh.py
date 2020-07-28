from node import Quadrant
from graphic import Segment
import numpy as np
import tkinter as tk
import math


class Grid(Quadrant):
    def __init__(self, x, y, width, height, cell_size):
        super(Grid, self).__init__(x, y, width, height)
        origin_x = x - cell_size * width / 2  # Top left corner x.
        origin_y = y - cell_size * height / 2  # Top left corner y.
        for i in range(self._height):
            self._leaves.append([])  # Add a row.
            yi = origin_y + (2 * i + 1) * cell_size / 2
            for j in range(self._width):
                xj = origin_x + (2 * j + 1) * cell_size / 2
                self._leaves[i].append(Quadrant(xj, yi, cell_size, cell_size))  # Add a column.
                self._leaves[i][j].canvas = canvas

    def draw(self, canvas, **kwargs):
        for i in range(self._height):
            for j in range(self._width):
                self._leaves[i][j].draw(canvas, **kwargs)


class Mesh2D:
    def __init__(self, x, y, width, height, cell_size, isovalue=0):
        self.cell_size = cell_size
        self.isovalue = isovalue  # Threshold for the data value of each pixel.
        self.pixels = Grid(x, y, width, height, cell_size)  # Grid of pixels.
        self.contours = Grid(x, y, width - 1, height - 1, cell_size)  # Grid of isolines.
        self.marked_pixel = None  # The vertex of a contour cell is highlighted when the mouse is over it.
        self.bits = []  # Contains the bit value for each pixel.
        self.vertices = []  # Contains the ovals representing the vertices of each contour cell.
        for i in range(height):
            self.bits.append([])
            self.vertices.append([])
            for j in range(width):
                pixels = self.pixels.leaves()
                value = pixels[i][j].contents()
                value.append(self.isovalue - 1)
                self.bits[i].append(self.binary_value(value[0]))
                cx, cy = pixels[i][j].get_center()
                radius = self.cell_size / 10
                self.vertices[i].append(canvas.create_oval(cx - radius, cy - radius, cx + radius, cy + radius))

    # Highlight a vertex of a contour cell.
    def mark_pixel(self, event):
        x, y = event.x, event.y
        center = self.pixels.get_center()
        i = math.floor((y - center[1] + self.cell_size * self.pixels.get_height() / 2) / self.cell_size)
        j = math.floor((x - center[0] + self.cell_size * self.pixels.get_width() / 2) / self.cell_size)

        # Check if the position of the mouse is within the boundary of the grid of pixels.
        if 0 <= i < self.pixels.get_height() and 0 <= j < self.pixels.get_width():
            pixels = self.pixels.leaves()
            cell = pixels[i][j]
            cx, cy = cell.get_center()
            radius = self.cell_size / 5
            if self.marked_pixel:
                canvas.coords(self.marked_pixel, cx - radius, cy - radius, cx + radius, cy + radius)
            else:
                self.marked_pixel = canvas.create_oval(cx - 5, cy - 5, cx + 5, cy + 5, outline="blue")
        else:
            canvas.delete(self.marked_pixel)
            self.marked_pixel = None

    # Return the bit value based on the isovalue.
    def binary_value(self, value):
        return '1' if value > self.isovalue else '0'

    # Invert the bit value of a pixel.
    def invert_bit(self, event):
        # Get the corresponding pixel.
        x, y = event.x, event.y
        center = self.pixels.get_center()
        i = math.floor((y - center[1] + self.cell_size * self.pixels.get_height() / 2) / self.cell_size)
        j = math.floor((x - center[0] + self.cell_size * self.pixels.get_width() / 2) / self.cell_size)

        # Invert the bit value.
        if 0 <= i < self.pixels.get_height() and 0 <= j < self.pixels.get_width():
            self.bits[i][j] = str(int(self.bits[i][j]) ^ 1)
            canvas.itemconfig(self.vertices[i][j], fill="black" if self.bits[i][j] == '1' else "")
            i = min(i, self.contours.get_height() - 1)
            j = min(j, self.contours.get_width() - 1)
            self.define_contour_cell(i, j)
            if i > 0 and j > 0:
                self.define_contour_cell(i - 1, j - 1)
                self.define_contour_cell(i - 1, j)
                self.define_contour_cell(i, j - 1)
            elif i > 0 and j == 0:
                self.define_contour_cell(i - 1, j)
            elif i == 0 and j > 0:
                self.define_contour_cell(i, j - 1)

    # Define the edges to be drawn inside the contour cell.
    def define_contour_cell(self, i, j):
        contours = self.contours.leaves()
        index = int(self.bits[i][j] +  # Top left corner.
                    self.bits[i][j+1] +  # Top right corner.
                    self.bits[i+1][j+1] +  # Bottom right corner.
                    self.bits[i+1][j], 2)  # Bottom left corner.
        cell = contours[i][j]
        edges = cell.contents()

        # Remove the current edges inside the cell.
        for edge in edges:
            canvas.delete(edge.get_item())
        x, y = contours[i][j].get_center()
        edges.clear()

        # Add the new edges inside the cell.
        edges += self.generate_edges(x, y, index)
        for edge in edges:
            edge.draw(canvas, fill="blue", width=2)

    # Return a list of edges based on the index obtained from the contour cell.
    def generate_edges(self, x, y, index):
        edges = []
        if index == 1 or index == 14:
            start = np.array([x, y + self.cell_size / 2])
            end = np.array([x - self.cell_size / 2, y])
            edges.append(Segment(start, end))
        elif index == 2 or index == 13:
            start = np.array([x, y + self.cell_size / 2])
            end = np.array([x + self.cell_size / 2, y])
            edges.append(Segment(start, end))
        elif index == 3 or index == 12:
            start = np.array([x - self.cell_size / 2, y])
            end = np.array([x + self.cell_size / 2, y])
            edges.append(Segment(start, end))
        elif index == 4 or index == 11:
            start = np.array([x, y - self.cell_size / 2])
            end = np.array([x + self.cell_size / 2, y])
            edges.append(Segment(start, end))
        elif index == 5:
            start = np.array([x, y - self.cell_size / 2])
            end = np.array([x - self.cell_size / 2, y])
            edges.append(Segment(start, end))
            start = np.array([x, y + self.cell_size / 2])
            end = np.array([x + self.cell_size / 2, y])
            edges.append(Segment(start, end))
        elif index == 6 or index == 9:
            start = np.array([x, y - self.cell_size / 2])
            end = np.array([x, y + self.cell_size / 2])
            edges.append(Segment(start, end))
        elif index == 7 or index == 8:
            start = np.array([x, y - self.cell_size / 2])
            end = np.array([x - self.cell_size / 2, y])
            edges.append(Segment(start, end))
        elif index == 10:
            start = np.array([x, y + self.cell_size / 2])
            end = np.array([x - self.cell_size / 2, y])
            edges.append(Segment(start, end))
            start = np.array([x, y - self.cell_size / 2])
            end = np.array([x + self.cell_size / 2, y])
            edges.append(Segment(start, end))
        return edges


if __name__ == "__main__":
    root = tk.Tk()
    canvas = tk.Canvas(root, width=1000, height=1000)
    canvas.pack()
    mesh = Mesh2D(500, 500, 10, 10, 50)
    mesh.contours.draw(canvas)
    canvas.bind('<Motion>', mesh.mark_pixel)
    canvas.bind('<Button-1>', mesh.invert_bit)
    root.mainloop()
