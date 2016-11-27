from diffuser import Diffuser
from rna import RNA
from process_seq import pwms, get_genes

def init_diffusers(x, y, z, n_u1, n_u2, n_u3, d_u1, d_u2, d_u3, num_steps):
    u1 = Diffuser(x, y, z, d_u1, n_u1)
    u2 = Diffuser(x, y, z, d_u2, n_u2)
    u3 = Diffuser(x, y, z, d_u3, n_u3)
    for _ in range(num_steps):
        u1.update()
        u2.update()
        u3.update()
    return u1, u2, u3



def main(x, y, z, num_steps, u1, u2, u3,
         pwm_u1, pwm_u2,
         u1_t, u2_t,
         persistence, sequence):
    rna = RNA(x, y, z, sequence, persistence, pwm_u1, pwm_u2, u1_t, u2_t)
    i = 0
    while i < num_steps:
        i += 1
        rna.extend()
        rna.random_flight()
        u1.update()
        u2.update()
        u3.update()
        rna.bind_snorna(u1, u2)
        rna.release_snorna(u1, u2)
        splice = rna.bind_u3(u3)
        if splice: return  i, splice
    return i, False

u1_pwm, u2_pwm = pwms()
persistence = 15
u1, u2, u3 = init_diffusers(50, 50, 50, 1000, 1000, 5000, 1, 1, 1, 500)
for seq, five, three in get_genes():
    f_block, t_block = (five+len(u1_pwm)-1) / persistence, (three+len(u2_pwm)-1) / persistence
    iterations, splice = main(50, 50, 50, int(len(seq) * 1.75), u1, u2, u3, u1_pwm, u2_pwm, 1.0, 1.0, persistence, seq)
    if not splice:
        print '\t'.join(map(str, [f_block, t_block, -1, -1]))
    else:
        print '\t'.join(map(str, [f_block, t_block, splice[0], splice[1]]))
