#!usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 26 16:07:26 2019

@author: xytan
"""

'''
Here are 4 basic functions that used to solve 
exam questions:
'''

# F1 build a dictionary of reads
def buildRead(path):
    '''
    This function read the records in certain .fasta file and 
    returns a dictionary in the form id:seq. 
    The argument is the path of target .fasta file.
    '''
    try:
        f = open(path)
    except IOError:
        print('The file does not exist')
    read = {}
    for line in f:
        # let's discard the newline at the end (if any!)
        line = line.rstrip()
        if line[0] == ">": #line.startswith('>')
            words = line.split()
            key = words[0][1:]
            read[key] = ''
        else:
            read[key] = read[key] + line
    f.close()
    return read

# F2 length of sequence
def seqlength(read):
    '''
    This function compute the sequence length of each record.
    It returns a dictionary in the form id:seqence length.
    The input argument should be a dict in the form id:seq.
    '''
    length = {}
    for i in read.keys():
        length[i] = len(read[i])
    return length

# F3 read frame
def ORF(read, position = 0):
    '''
    This function identifies all potential reading frames
    in the DNA sequences based on a start codon and 
    a stop codon.
    It returns a tuple where the 1st 
    element is a list of record ids (may have repeats), 
    the 2nd element is a list of all reading frames 
    with given reading frame.
    
    The arguments includes 1) a dict in the form id:seq, 
    2) a starting reading position(0,1,2 for reading frame 1,2,3
    respectively).
    '''
    name = []
    orf = []
    for k in read.keys():
        seq = read[k]
        start = position
        for i in range(start, len(seq), 3):
            frame = ''
            if seq[i:i+3] != 'ATG':
                continue
            name.append(k)
            frame = frame + seq[i:i+3]
            for j in range(i+3, len(seq), 3):
                frame = frame + seq[j:j+3]
                if seq[j:j+3] in ("TAA", "TAG", "TGA"):
                    start = j+3
                    orf.append(frame)
                    break
            else:
                name.pop()
    return name, orf

# F4 repeart region
def repeat(n, string):
    '''
    This function figures out all repeats with a given length 
    and count the times of repeats for each repeat.
    
    It returns a tuple where the 1st element is a list of
    the nucleotide sequences of all repeats and the 2nd 
    element is a list of the times of repeats for each repeat.   
    
    Specially, this function captures the repeats overlapping 
    themselves. 
    
        i.e.: repeat(3, 'ACACACACA') 
              -> (['ACA', 'CAC'], [4, 3])
              
    The input arguments include 1) the length of repeats,
    n, 2) a single DNA seq.
    
    It can also compute the position of the repeats
    and return the starting index, but it is not involved in
    the exam.
    '''
#    index_all = []
    repeats = []
    seq = []
    # initiate the search (loop1):
    for i in range(len(string)):
        if string[i:i+n] not in seq:
        # if the seq appeared before, skip to avoid redundent.
        # Does the seq repeat at leat once? 
            find = string.find(string[i:i+n], i+1)
            if find != -1:
#                    index = [i, find]
                num = 2 # at least 2 bcz the seq is repeating
                # Does the seq repeats 2 or more (loop2)?
                j = find+1
                while j < len(string):
                    find = string.find(string[i:i+n], j)
                    if find != -1:
#                            index.append(find)
                        j = find + 1
                        num +=1
                    else:
                        j+=1 # an infinite loop if forget this.
                # record times of repeats and the seq after searching
                repeats.append(num)
                seq.append(string[i:i+n])
#                    index_all.append(index)
#    print('repeat sequence, index, times of repeats')
    return (seq, repeats)

# F5 merge the times of repeats of identical seqs
def merge(repeat_seq, times_of_repeat):
    '''
    This function merges the times of repeats corresponding to 
    identical repeat regions by adding up thes times of repeats.
    
    It returns a tuple where the 1st element is a list of the 
    DNA seqs of unique repaet regions and the 2nd element is a 
    list of corresponding times of repeats.
    
    The input arguments includes 1) a list of all repeat regions 
    in DNA seqs from all records, 2) a list of corresponding
    times of repeats.    
    '''
    repeat_seq_copy = repeat_seq
    uni_seq = list(set(repeat_seq))
    merge_repeat = []
    for seq in uni_seq:
        num = 0
        for i in repeat_seq_copy:
            if i == seq:
                index = repeat_seq_copy.index(i) # always return the first one
                num += times_of_repeat[index]
                repeat_seq_copy[index] = ''
        merge_repeat.append(num)
    return uni_seq, merge_repeat

'''
Here are the answers to the questions:
'''
# Q1 number of records
path=".\Downloads\dna2.fasta"
read = buildRead(path)
len(read)

#Q2, Q3 the longest/shortest sequence
length = seqlength(read)
max(list(length.values()))
min(list(length.values()))

# Q4 length of longest OFR in reading frame 2 ?
frame2 = ORF(read, 1)
frame2[1]
ORFlength = [len(i) for i in frame2[1]]
max(ORFlength)

# Q5 ORF starting position
frame3 = ORF(read, 2)
ORFlength = [len(i) for i in frame3[1]]
index_long = ORFlength.index(max(ORFlength))
key = frame3[0][index_long]
read[key].find(frame3[1][index_long])

# Q6 length of longest ORF in any fw reading frame
for i in range(0,3):
    print('frame', i+1)
    frame = ORF(read, i)
    ORFlength = [len(i) for i in frame[1]]
    print(max(ORFlength))

# Q7 length of longest ORF with specific id
ID = 'gi|142022655|gb|EQ086233.1|16'
for i in range(0,3):
    print('frame', i+1)
    frame = ORF(read, i)
    index = [i for i in range(0, len(frame[0])) if frame[0][i] == ID]
    ORFlength = [len(frame[1][i]) for i in index]
    print(max(ORFlength))

# Q8 number of most freq repeats, n=6
all_seq = [];all_repeat = []
for s in list(read.values()):
    result = repeat(6, s)
    seq = result[0]
    num_repeat = result[1]
    all_seq.extend(seq)
    all_repeat.extend(num_repeat)

result = merge(all_seq, all_repeat)
max(result[1])

# Q9 n =12
all_seq = [];all_repeat = []
for s in list(read.values()):
    result = repeat(12, s)
    seq = result[0]
    num_repeat = result[1]
    all_seq.extend(seq)
    all_repeat.extend(num_repeat)
    
result = merge(all_seq, all_repeat)
maximum = max(result[1])    
result[1].count(maximum)

# Q10 n = 7, seq = ?
all_seq = [];all_repeat = []
for s in list(read.values()):
    result = repeat(7, s)
    seq = result[0]
    num_repeat = result[1]
    all_seq.extend(seq)
    all_repeat.extend(num_repeat)

result = merge(all_seq, all_repeat)
maximum = max(result[1])  
index = result[1].index(maximum)  
result[0][index]
