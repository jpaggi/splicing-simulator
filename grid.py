
class Grid:
    """
    A 3D grid with objects able to be stored in each grid space.
    The Grid wraps around at the edges
    """
    def __init__(self, y, x, z):
        self.y = y
        self.x = x
        self.z = z
        self.blocks = {}

    def block_index(self):
        return self.blocks.keys()
    
    def nextX(self, x, step):
        return (x+step) % self.x

    def nextY(self, y, step):
        return (y+step) % self.y

    def nextZ(self, z, step):
        return (z+step) % self.z

    def get_block(self, x, y, z):
        if not (x, y, z) in self.blocks: return 0
        return self.blocks[(x, y, z)]

    def add_to_block(self, x, y, z, num=1):
        if (x, y, z) not in self.blocks: self.blocks[(x, y, z)] = 0
        self.blocks[(x, y, z)] += num

    def remove_from_block(self, x, y, z, num = 1):
        self.blocks[(x, y, z)] += -1

    def neighbors(self, x, y, z, dist):
        return [self.blocks[(self.nextX(x, i), self.nextY(y, j), self.nextZ(z, k))]
                for i in range(-dist, dist+1)
                for j in range(-dist, dist+1)
                for k in range(-dist, dist+1)]
        

