from grid import Grid
from random import random
from math import ceil, log, exp

def best_score(seq, pwm):
    best, lor = -1, 0
    for start in range(len(seq) - len(pwm)):
        s = 0
        for i in range(len(pwm)):
            s += log(pwm[i][seq[start+i]] / .25) if pwm[i][seq[start+i]] else - float('inf')
        if s > lor:
            best, lor = start, s
    return best, lor

class RNA:
    def __init__(self, x, y, z,
                 sequence, persistence, u1_pwm, u2_pwm, u1_t, u2_t):
        self.x = x
        self.y = y
        self.z = z
        self.sequence = sequence
        self.persistence = persistence
        self.u1_pwm = u1_pwm
        self.u2_pwm = u2_pwm
        self.u1_t = u1_t
        self.u2_t = u2_t

        self.length = 0
        self.bound_u1 = []
        self.bound_u2 = []
        self.blocks = []
        self.grid = Grid(x, y, z) # just to use next*(...)

    def extend(self):
        self.length = min(len(self.sequence), self.length + 1)

    def segments(self):
        return range(int(ceil(self.length / float(self.persistence))))

    def random_flight(self):
        self.blocks = [(0, 0, 0)]
        for _ in self.segments():
            delta = [(1, 0, 0), (-1, 0, 0),
                     (0, 1, 0), (0, -1, 0),
                     (0, 0, 1), (0, 0, -1)][int(random() * 6)]
            last = self.blocks[-1]
            self.blocks += [(self.grid.nextX(last[0], delta[0]),
                             self.grid.nextY(last[1], delta[1]),
                             self.grid.nextZ(last[2], delta[2]))]

    def bind_snorna(self, u1, u2):
        for i in self.segments():
            if i in self.bound_u1 or i in self.bound_u2: continue
            block = self.blocks[i]
            if u1.get_block(*block):
                seq = self._get_block_seq(i, len(self.u1_pwm))
                best, score = best_score(seq, self.u1_pwm)
                if score > 0:
                    self.bound_u1 += [i]
                    u1.remove_from_block(*block)
            if u2.get_block(*block):
                seq = self._get_block_seq(i, len(self.u2_pwm))
                best, score = best_score(seq, self.u2_pwm)
                if score > 0:
                    self.bound_u2 += [i]
                    u2.remove_from_block(*block)
    
    def release_snorna(self, u1, u2):
        for i in self.segments():
            block = self.blocks[i]
            if i in self.bound_u1:
                seq = self._get_block_seq(i, len(self.u1_pwm))
                best, score = best_score(seq, self.u1_pwm)
                if random() < exp(-score / self.u1_t):
                    self.bound_u1.remove(i)
                    u1.add_to_block(*block)

            if i in self.bound_u2:
                seq = self._get_block_seq(i, len(self.u2_pwm))
                best, score = best_score(seq, self.u2_pwm)
                if random() < exp(-score / self.u2_t):
                    self.bound_u2.remove(i)
                    u2.add_to_block(*block)

    def bind_u3(self, u3):
        for i in self.segments():
            if not i in self.bound_u1: continue
            for j in self.segments():
                if self.blocks[i] != self.blocks[j]: continue
                block = self.blocks[i]
                if i in self.bound_u1 and j in self.bound_u2 and u3.get_block(*block):
                    return (i, j)
        return False

    def _get_block_seq(self, i, motif_length):
        begin = i * self.persistence
        if i: begin += - motif_length + 1
        end = min((i+1) * self.persistence, self.length)
        return self.sequence[begin: end]
