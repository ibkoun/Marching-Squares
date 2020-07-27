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

    def draw(self, canvas, fill="", outline="black"):
        for i in range(self._height):
            for j in range(self._width):
                self._leaves[i][j].draw(canvas, fill=fill, outline=outline)


class Mesh2D:
    def __init__(self, x, y, width, height, cell_size):
        self.cell_size = cell_size
        self.pixels = Grid(x, y, width, height, cell_size)  # Grid of pixels.
        self.contours = Grid(x, y, width - 1, height - 1, cell_size)  # Grid of isolines or isobands.
        self.marked_pixel = None  # The corner of a contour cell is highlighted when the mouse is over it.
        self.corners = []  # Contains the ovals representing the corners of each contour cell.
        for i in range(height):
            self.corners.append([])
            for j in range(width):
                pixels = self.pixels.leaves()
                bit = pixels[i][j].contents()
                bit.append(False)
                cx, cy = pixels[i][j].get_center()
                self.corners[i].append(canvas.create_oval(cx - 5, cy - 5, cx + 5, cy + 5))

    # Highlight a corner of a contour cell.
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
            if self.marked_pixel:
                canvas.coords(self.marked_pixel, cx - 5, cy - 5, cx + 5, cy + 5)
            else:
                self.marked_pixel = canvas.create_oval(cx - 5, cy - 5, cx + 5, cy + 5, outline="blue")
        else:
            canvas.delete(self.marked_pixel)
            self.marked_pixel = None

    # Set the bit value of a pixel.
    def set_pixel(self, event):
        x, y = event.x, event.y
        center = self.pixels.get_center()
        i = math.floor((y - center[1] + self.cell_size * self.pixels.get_height() / 2) / self.cell_size)
        j = math.floor((x - center[0] + self.cell_size * self.pixels.get_width() / 2) / self.cell_size)
        if 0 <= i < self.pixels.get_height() and 0 <= j < self.pixels.get_width():
            pixels = self.pixels.leaves()
            contents = pixels[i][j].contents()
            contents[0] = not contents[0]
            bit = contents[0]
            canvas.itemconfig(self.corners[i][j], fill="black" if bit else "")
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
        pixels = self.pixels.leaves()
        index = int(str(int(pixels[i][j].contents()[0])) +  # Top left corner.
                    str(int(pixels[i][j+1].contents()[0])) +  # Top right corner.
                    str(int(pixels[i+1][j+1].contents()[0])) +  # Bottom right corner.
                    str(int(pixels[i+1][j].contents()[0])), 2)  # Bottom left corner.
        cell = contours[i][j]
        edges = cell.contents()
        for edge in edges:
            canvas.delete(edge.item)
        x, y = contours[i][j].get_center()
        edges.clear()
        edges += self.generate_edges(x, y, index)
        for edge in edges:
            edge.draw(canvas, fill="black")

    # Return a list of edges based on the index obtained from the contour cell.
    def generate_edges(self, x, y, case):
        edges = []
        if case == 1 or case == 14:
            start = np.array([x, y + self.cell_size / 2])
            end = np.array([x - self.cell_size / 2, y])
            edges.append(Segment(start, end))
        elif case == 2 or case == 13:
            start = np.array([x, y + self.cell_size / 2])
            end = np.array([x + self.cell_size / 2, y])
            edges.append(Segment(start, end))
        elif case == 3 or case == 12:
            start = np.array([x - self.cell_size / 2, y])
            end = np.array([x + self.cell_size / 2, y])
            edges.append(Segment(start, end))
        elif case == 4 or case == 11:
            start = np.array([x, y - self.cell_size / 2])
            end = np.array([x + self.cell_size / 2, y])
            edges.append(Segment(start, end))
        elif case == 5:
            start = np.array([x, y - self.cell_size / 2])
            end = np.array([x - self.cell_size / 2, y])
            edges.append(Segment(start, end))
            start = np.array([x, y + self.cell_size / 2])
            end = np.array([x + self.cell_size / 2, y])
            edges.append(Segment(start, end))
        elif case == 6 or case == 9:
            start = np.array([x, y - self.cell_size / 2])
            end = np.array([x, y + self.cell_size / 2])
            edges.append(Segment(start, end))
        elif case == 7 or case == 8:
            start = np.array([x, y - self.cell_size / 2])
            end = np.array([x - self.cell_size / 2, y])
            edges.append(Segment(start, end))
        elif case == 10:
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
    mesh = Mesh2D(500, 500, 5, 5, 50)
    mesh.contours.draw(canvas)
    canvas.bind('<Motion>', mesh.mark_pixel)
    canvas.bind('<Button-1>', mesh.set_pixel)
    root.mainloop()
