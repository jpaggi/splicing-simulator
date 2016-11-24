from grid import Grid
from random import random

class Diffuser:
    def __init__(self, x, y, z, constant = 1, num = 1):
        self.x = x
        self.y = y
        self.z = z
        self.grid = Grid(x, y, z)
        self.grid.add_to_block(0, 0, 0, num)
        self.constant = constant

    def get_block(self, x, y, z):
        return self.grid.get_block(x, y, z)

    def add_to_block(self, x, y, z, num = 1):
        self.grid.add_to_block(x, y, z, num)

    def remove_from_block(self, x, y, z, num = 1):
        self.grid.remove_from_block(x, y, z, num)
        
    def update(self):
        new_grid = Grid(self.x, self.y, self.z)
        step = self.constant
        for block in self.grid.block_index():
            for _ in range(self.grid.get_block(*block)):
                delta = [(step, 0, 0), (-step, 0, 0),
                         (0, step, 0), (0, -step, 0),
                         (0, 0, step), (0, 0, -step)][int(random() * 6)]
                new_grid.add_to_block(self.grid.nextX(block[0], delta[0]),
                                      self.grid.nextY(block[1], delta[1]),
                                      self.grid.nextZ(block[2], delta[2]))
        self.grid = new_grid
