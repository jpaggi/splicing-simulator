from grid import Grid
from random import random
from math import ceil

def best_score(seq, pwm):
    best, score = -1, 0
    for start in range(len(seq) - len(pwm)):
        s = 1
        for i in range(len(pwm)):
            s *= pwm[i][seq[start+i]]
        if s > score:
            best, score = start, s

    best_pos = 1
    worst_pos = 1
    for i in range(len(pwm)):
        best_pos *= max(pwm[i].values())
        worst_pos *= min(pwm[i].values())

    return best, (score - worst_pos) / (best_pos - worst_pos)

class RNA:
    def __init__(self, x, y, z,
                 sequence, persistence, pwm_u1, pwm_u2):
        self.x = x
        self.y = y
        self.z = z
        self.sequence = sequence
        self.persistence = persistence
        self.pwm_u1 = pwm_u1
        self.pwm_u2 = pwm_u2

        self.length = 0
        self.bound_u1 = []
        self.bound_u2 = []
        self.blocks = []
        self.grid = Grid(x, y, z) # just to use next*(...)

    def extend(self):
        self.length += 1

    def random_flight(self):
        self.blocks = [(0, 0, 0)]
        for _ in range(int(ceil(self.length / float(self.persistence)))):
            delta = [(1, 0, 0), (-1, 0, 0),
                     (0, 1, 0), (0, -1, 0),
                     (0, 0, 1), (0, 0, -1)][int(random() * 6)]
            self.blocks += [tuple([i+j for i, j in zip(
                            list(self.blocks[-1]), list(delta))])]

    def bind_snorna(self, u1, u2):
        for i in range(int(ceil(self.length / float(self.persistence)))):
            block = self.blocks[i]
            if u1.get_block(*block) and block not in self.bound_u1:
                seq = self._get_block_seq(i, len(self.pwm_u1))
                best, score = best_score(seq, self.pwm_u1)
                if score:
                    self.bound_u1 += [i]
                    u1.remove_from_block(*block)
            if u2.get_block(*block) and block not in self.bound_u2:
                seq = self._get_block_seq(i, len(self.pwm_u2))
                best, score = best_score(seq, self.pwm_u2)
                if score:
                    self.bound_u2 += [i]
                    u2.remove_from_block(*block)
    
    def release_snorna(self, u1, u2):
        for i in range(int(ceil(self.length / float(self.persistence)))):
            block = self.blocks[i]
            if block in self.bound_u1:
                seq = self._get_block_seq(i, len(u1_pwm))
                best, score = best_score(seq, u1_pwm)
                if random() > score:
                    self.bound_u1.remove(i)
                    u1.add_to_block(*block)

            if block in self.bound_u2:
                seq = self._get_block_seq(i, len(u2_pwm))
                best, score = best_score(seq, u2_pwm)
                if random() > score:
                    self.bound_u2.remove(i)
                    u2.add_to_block(*block)

    def bind_u3(self, u3):
        for i in range(int(ceil(self.length / float(self.persistence)))):
            for j in range(int(ceil(self.length / float(self.persistence)))):
                if self.blocks[i] != self.blocks[j]: continue
                block = self.blocks[i]
                if i in self.bound_u1 and j in self.bound_u2 and u3.get_block(*block):
                    return i, j
        return False

    def _get_block_seq(self, i, motif_length):
        begin = i * self.persistence
        if i: begin += - motif_length - 1
        end = min((i+1) * self.persistence, self.length)
        if end != self.length: end = min(end + motif_length -1, self.length)

        return self.sequence[begin: end]
