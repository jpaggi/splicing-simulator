def revcomp(seq):
    comps = {'A': 'T', 'G':'C', 'T':'A', 'C':'G'}
    return ''.join(comps[char] for char in seq[::-1])

def get_seq(gff):
    gff = open(gff, 'r')
    line = gff.readline()
    while line[0] != '>':
        line = gff.readline()

    chrom = line[1:].strip()
    seq = {chrom:''}
    for line in gff:
        line = line.strip()
        if line[0] == '>':
            chrom = line[1:]
            seq[chrom] = ''
        else:
            seq[chrom] += line.upper()
    gff.close()
    return seq

def get_splice_sites(gff):
    gff = open(gff, 'r')
    line = gff.readline()
    while line[0] == '#': line = gff.readline()

    line = line.strip().split()
    fives, threes = [], []
    while len(line) > 7:
        chrom, source, feature, start, end, a, strand = line[:7]
        if feature == 'intron':
            if strand == '+':
                five, three = int(start)-1, int(end)
            else:
                three, five = int(start)-1, int(end)
            fives += [(chrom, five, strand)]
            threes += [(chrom, three, strand)]
        line = gff.readline().strip().split()
    gff.close()
    return fives, threes

def get_splice_site_seqs(ss, up_len, down_len, seq):
    seqs = []
    for chrom, site, strand in ss:
        if strand == '+':
            seqs += [seq[chrom][site - up_len: site + down_len]]
        else:
            seqs += [revcomp(seq[chrom][site - down_len: site + up_len])]
    return seqs

def get_pwm(seqs):
    counts = [{'A':0, 'T':0, 'C':0, 'G':0} for _ in range(len(seqs[0]))]
    for seq in seqs:
        assert len(seq) == len(counts)
        for i, char in enumerate(seq):
            counts[i][char] += 1
    for d in range(len(counts)):
        total = float(sum(counts[d].values()))
        for key, val in counts[d].items():
            counts[d][key] = val / total
    return counts

def pwms():
    gff = 'saccharomyces_cerevisiae.gff'
    seq = get_seq(gff)
    fives, threes = get_splice_sites(gff)
    five_seqs = get_splice_site_seqs(fives, 3, 7, seq)
    three_seqs = get_splice_site_seqs(threes, 15, 3, seq)
    five_seqs = filter(lambda x: x[3:5] == 'GT', five_seqs)
    three_seqs = filter(lambda x: x[13:15] == 'AG', three_seqs)
    five_pwm = get_pwm(five_seqs)
    three_pwm = get_pwm(three_seqs)
    return five_pwm, three_pwm
