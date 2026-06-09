import numpy as np
from itertools import combinations

def fft_similarity(seq1, seq2):
    n = len(seq1) + len(seq2)
    s1 = np.zeros(n)
    s2 = np.zeros(n)
    for i, c in enumerate(seq1.upper()):
        s1[i] = ord(c)
    for i, c in enumerate(seq2.upper()):
        s2[i] = ord(c)
    f1 = np.fft.fft(s1)
    f2 = np.fft.fft(s2)
    corr = np.fft.ifft(f1 * np.conj(f2))
    return float(np.max(np.abs(corr))) / (len(seq1) * len(seq2) + 1)

def nw_align(seq1, seq2, match=1, mismatch=-1, gap=-2):
    n, m = len(seq1), len(seq2)
    dp = np.zeros((n+1, m+1))
    for i in range(n+1): dp[i][0] = i * gap
    for j in range(m+1): dp[0][j] = j * gap
    for i in range(1, n+1):
        for j in range(1, m+1):
            s = match if seq1[i-1] == seq2[j-1] else mismatch
            dp[i][j] = max(dp[i-1][j-1]+s, dp[i-1][j]+gap, dp[i][j-1]+gap)
    a1, a2, i, j = [], [], n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0:
            s = match if seq1[i-1] == seq2[j-1] else mismatch
            if dp[i][j] == dp[i-1][j-1]+s:
                a1.append(seq1[i-1]); a2.append(seq2[j-1]); i -= 1; j -= 1; continue
        if i > 0 and (j == 0 or dp[i][j] == dp[i-1][j]+gap):
            a1.append(seq1[i-1]); a2.append('-'); i -= 1
        else:
            a1.append('-'); a2.append(seq2[j-1]); j -= 1
    return ''.join(reversed(a1)), ''.join(reversed(a2))

def build_guide_tree(sequences):
    names = list(sequences.keys())
    n = len(names)
    dist = {}
    for i, j in combinations(range(n), 2):
        sim = fft_similarity(sequences[names[i]], sequences[names[j]])
        dist[(i,j)] = 1.0 / (sim + 1e-9)
    clusters = [[i] for i in range(n)]
    order = []
    while len(clusters) > 1:
        best = min(
            [(ci, cj) for ci in range(len(clusters)) for cj in range(ci+1, len(clusters))],
            key=lambda p: np.mean([dist.get((min(a,b),max(a,b)),0)
                                   for a in clusters[p[0]] for b in clusters[p[1]]])
        )
        ci, cj = best
        order.append((list(clusters[ci]), list(clusters[cj])))
        clusters[ci] = clusters[ci] + clusters[cj]
        del clusters[cj]
    return order, names

def consensus_seq(profile):
    length = len(profile[0])
    result = []
    for i in range(length):
        col = [s[i] for s in profile if s[i] != '-']
        result.append(max(set(col), key=col.count) if col else '-')
    return ''.join(result)

def reinsert_gaps(seqs, aligned_consensus):
    result = []
    for seq in seqs:
        new_seq, si = [], 0
        for c in aligned_consensus:
            if c == '-':
                new_seq.append('-')
            else:
                new_seq.append(seq[si] if si < len(seq) else '-')
                si += 1
        result.append(''.join(new_seq))
    return result

def mafft(sequences):
    if not sequences:
        return {}
    if len(sequences) == 1:
        return dict(sequences)
    order, names = build_guide_tree(sequences)
    profiles = {i: [sequences[names[i]]] for i in range(len(names))}
    for group1, group2 in order:
        p1 = [profiles[i][0] for i in group1]
        p2 = [profiles[i][0] for i in group2]
        con1 = consensus_seq(p1)
        con2 = consensus_seq(p2)
        ac1, ac2 = nw_align(con1, con2)
        new_p1 = reinsert_gaps(p1, ac1)
        new_p2 = reinsert_gaps(p2, ac2)
        for idx, ni in enumerate(group1):
            profiles[ni] = [new_p1[idx]]
        for idx, ni in enumerate(group2):
            profiles[ni] = [new_p2[idx]]
    return {names[i]: profiles[i][0] for i in range(len(names))}
