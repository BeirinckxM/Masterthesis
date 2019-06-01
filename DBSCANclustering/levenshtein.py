import numpy as np

def minimumEditDistance(s1,s2):
    if len(s1) > len(s2):
        s1,s2 = s2,s1
    distances = range(len(s1) + 1)
    for index2,char2 in enumerate(s2):
        newDistances = [index2+1]
        for index1,char1 in enumerate(s1):
            if char1 == char2:
                newDistances.append(distances[index1])
            else:
                newDistances.append(1 + min((distances[index1],
                                             distances[index1+1],
                                             newDistances[-1])))
        distances = newDistances
    return distances[-1]
 
print(minimumEditDistance("kitten","sitting"))
print(minimumEditDistance("rosettacode","raisethysword"))

def create_dist(column):
    x = column
    y = column
    
    matrix = np.zeros((len(x),len(y)),dtype=np.int)

    for i in range(0,len(x)):
      for j in range(0,len(y)):
        matrix[i,j] = minimumEditDistance(x[i],y[j])
        
    return matrix

def add_CW(seq):
    seq=str(seq)
    if seq[0] != 'C' and seq[-1] != 'W':
        seq = 'C' + seq + 'W'
    if seq[0] != 'C' and seq[-1] == 'W':
        seq = 'C' + seq 
    if seq[0] == 'C' and seq[-1] != 'W':
        seq = seq + 'W'
    return seq

def add_CW_strict(seq):
    seq=str(seq)
    newseq = 'C' + seq + 'W'
    return newseq
