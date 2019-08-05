#!usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 15:44:48 2019

@author: xytan
"""

'''Download a fasta file'''
import requests
url = 'https://d28rh4a8wq0iu5.cloudfront.net/ads1/data/chr1.GRCh38.excerpt.fasta'
r = requests.get(url)
with open('c4w3.fasta', 'wb') as f:
    f.write(r.content)
    f.close()

'''Read fasta file'''
def readFasta(filename):
    sequences = ''
    with open(filename) as f:
        for line in f:
            if line[0] != '>':
                sequences+=line.rstrip()
    f.close()
    return sequences

def approximateMatch_ed(x, y):
    # Create distance matrix
    D = []
    for i in range(len(x)+1):
        D.append([0]*(len(y)+1))
    # Initialize first row and column of matrix
    for i in range(len(x)+1):
        D[i][0] = i
    for i in range(len(y)+1):
        D[0][i] = 0
    # Fill in the rest of the matrix
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            distHor = D[i][j-1] + 1
            distVer = D[i-1][j] + 1
            if x[i-1] == y[j-1]:
                distDiag = D[i-1][j-1]
            else:
                distDiag = D[i-1][j-1] + 1
            D[i][j] = min(distHor, distVer, distDiag)
    # Edit distance is the value in the bottom right corner of the matrix
    return min(D[-1])

# Q1&2
p = 'GCTGATCGATCGTACG'
t = readFasta('c4w3.fasta')
p = 'GATTTACCAGATTGAG'
approximateMatch_ed(p,t)

'''Download a fastq file'''
import requests
url = 'https://d28rh4a8wq0iu5.cloudfront.net/ads1/data/ERR266411_1.for_asm.fastq'
r = requests.get(url)
with open('c4w3.fastq', 'wb') as f:
    f.write(r.content)
    f.close()

'''Read fasta file'''
def readFastq(filename):
    sequences = []
    with open(filename) as f:
        while True:
            f.readline()
            seq = f.readline().rstrip()
            f.readline()
            f.readline()
            if len(seq) == 0:
                break
            sequences.append(seq)
    return sequences

def overlap(a, b, min_length=3):
    """ Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's prefix in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match
        
def kmerDict(reads, k):
    '''create a dict {k-mer:set(reads)}'''
    d = {}
    for r in reads:
        for i in range(len(r) - k +1):
            if r[i:i+k] not in list(d.keys()):
                d[r[i:i+k]] = {r}
            d[r[i:i+k]].add(r)
    return d

def overlap_all_pairs(reads, min_length = 30):
    '''find all overlap pairs'''
    d = kmerDict(reads, min_length)
    olaps = []
    for a in reads:
        suffix = a[-min_length:]
        reads_b = d[suffix]
        for b in reads_b:
            if b != a and overlap(a, b, min_length) > 0:   
                olaps.append((a, b))
    return olaps

#Q3
reads = readFastq('c4w3.fastq')
pairs = overlap_all_pairs(reads, 30)

#Q4
nodes = set()
for p in pairs:
    nodes.add(p[0])
