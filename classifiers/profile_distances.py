#!/usr/bin/env python3

import numpy as np

from scipy.stats import entropy
from scipy.stats import zscore



#Z-normalization of dictionaries
#Using zscore from scipy.stats. d stands for dictionary? ddof = degress of freedom.
#zscore = (x - mu)/sigma
def znorm_scipy(d):
    keys, vals = zip(*d.items())
    #variadic argument, creates list with (key, value)
    return dict(zip(keys, zscore(vals, ddof=1)))
    #Zip to a dictionary again.
    #ddof = 1, does this mean n-1?

#Weighted euclidean
def weightedL2(a,b,w):
    q = a-b
    return np.sqrt((w*q*q).sum())


# chemical properties (dictionaries)
basicity = {'A': 206.4, 'B': 210.7, 'C': 206.2, 'D': 208.6, 'E': 215.6, 'F': 212.1, 'G': 202.7,
            'H': 223.7, 'I': 210.8, 'K': 221.8, 'L': 209.6, 'M': 213.3, 'N': 212.8, 'P': 214.4,
            'Q': 214.2, 'R': 237.0, 'S': 207.6, 'T': 211.7, 'V': 208.7, 'W': 216.1, 'X': 210.2,
            'Y': 213.1, 'Z': 214.9}

hydrophobicity = {'A': 0.16, 'B': -3.14, 'C': 2.50, 'D': -2.49, 'E': -1.50, 'F': 5.00, 'G': -3.31,
                  'H': -4.63, 'I': 4.41, 'K': -5.00, 'L': 4.76, 'M': 3.23, 'N': -3.79, 'P': -4.92,
                  'Q': -2.76, 'R': -2.77, 'S': -2.85, 'T': -1.08, 'V': 3.02, 'W': 4.88, 'X': 4.59,
                  'Y': 2.00, 'Z': -2.13}

helicity = {'A': 1.24, 'B': 0.92, 'C': 0.79, 'D': 0.89, 'E': 0.85, 'F': 1.26, 'G': 1.15, 'H': 0.97,
            'I': 1.29, 'K': 0.88, 'L': 1.28, 'M': 1.22, 'N': 0.94, 'P': 0.57, 'Q': 0.96, 'R': 0.95,
            'S': 1.00, 'T': 1.09, 'V': 1.27, 'W': 1.07, 'X': 1.29, 'Y': 1.11, 'Z': 0.91}

molecular_weight ={'A': 89.1, 'R': 174.2, 'N': 132.1, 'D': 133.1, 'C': 121.2, 'E': 147.1, 'Q': 146.2, 
                   'G': 75.1, 'H': 155.2, 'I': 131.2, 'L': 131.2, 'K': 146.2, 'M': 149.2, 'F': 165.2, 
                   'P': 115.1, 'S': 105.1, 'T': 119.1, 'W': 204.2, 'Y': 181.2, 'V': 117.1, 'Z': 147.0, 
                   'B': 133.0,'X': 140.0}

mutation_stability = {'A': 13, 'C': 52, 'D': 11, 'E': 12, 'F': 32, 'G': 27, 'H': 15, 'I': 10,
                      'K': 24, 'L': 34, 'M':  6, 'N':  6, 'P': 20, 'Q': 10, 'R': 17, 'S': 10,
                      'T': 11, 'V': 17, 'W': 55, 'Y': 31, 'X': 21}

elektrochem = {'molecular_weight': znorm_scipy(molecular_weight),'basicity': znorm_scipy(basicity),
               'hydrophobicity': znorm_scipy(hydrophobicity), 'helicity': znorm_scipy(helicity),
               'mutation_stability': znorm_scipy(mutation_stability)}

def make_profile(sequence, prop = 'basicity'):
    #Make the respective profile for the amino acids in a given sequence.
    #prop is the dictionary key with elektrochem as dictionary.
    #x is the amino acid in the respective sequence,
    #which to refers to the value for that amino acid in the normalized dictionary.
    return [elektrochem[prop][x] for x in sequence]

def distribution_similarity_realign(distriA, distriB):

    """Calculate the similarity between two distributions.
    Similarity is based on the distance measure used in the calc_distance function.
    If the distributions are of unequal length, slide the shorter one over the longer one
    and return the best similarity match.
    Return the best similarity score and the start position of the shorter distribution
    on the longer one that resulted in the best similarity score.
    """
    # First of all what are these distributions? Just sequences?
    # if distributions are of equal length, return similarity
    if len(distriA)==len(distriB):
        return calc_distance(distriA, distriB), 0 #What's the O here? Return False?

    # if distributions are of unequal length
    else:
        # determine which is the longer distribution
        if len(distriA) > len(distriB):
            long_distri = distriA
            short_distri= distriB
        else:
            long_distri = distriB
            short_distri= distriA

        # slide shorter distribution over longer one, calculating the similarity between them
        # initialize first similarity value
        #print('Calculating initial similarity for position 0:\n{}\n{}.'.format(short_distri, long_distri[0:len(short_distri)]))
        similarity = calc_distance(short_distri, long_distri[0:len(short_distri)])
        pos = 0 #Is this necessary?
        for i in range(1, len(long_distri)-len(short_distri)+1): #aligning and starting at pos = 1 and ending at + 1, as initial similarity was pos = 0
            #len(long_distri)-len(short_distri)+1) overlap between both distributions, only need to advance over this difference.
            #However, is it possible to advance in the opposite direction?
            #print('Calculating similarity for position {}:\n{}\n{}.'.format(i,short_distri, long_distri[i:i+len(short_distri)]))
            new_similarity = calc_distance(short_distri, long_distri[i:i+len(short_distri)])
            if new_similarity < similarity:
                similarity = new_similarity
                pos = i
        return similarity, pos #The best similarity score and it's starting position in the longest distibution


def distribution_similarity(distriA, distriB):

    # if distributions are of equal length, return similarity
    if len(distriA)==len(distriB):
        return calc_distance(distriA, distriB), 0 #What's the O?

    else:
        # determine which is the longer distribution
        if len(distriA) > len(distriB):
            long_distri = distriA
            short_distri= distriB
        else:
            long_distri = distriB
            short_distri= distriA

        diff = len(long_distri) - len(short_distri)

        if diff % 2 == 0: #If difference is even, amount of AA 'unmatched'
            shift = int(diff/2) #integer!
            return calc_distance(short_distri, long_distri[shift:len(long_distri)-shift]), diff/2

            '''For example diff = 19 - 14 = 5, shift = int(2,5) = 2.
            Starting at index 2 in the long distribution and ending at -2'''
            #centering of short sequence in the long sequence

        else: #If difference is odd

            #print("Odd distance : ",diff)
            shift = int(diff/2)
            #print("Shift : ",shift)

            #print("Short length : ", len(short_distri))
            #print("Old length : ", len(long_distri))
            #print("New length : ", len(long_distri[shift+1:len(long_distri)-shift]))

            #Two possibilities? Left or right centered?
            first = calc_distance(short_distri, long_distri[shift+1:len(long_distri)-shift])
            second = calc_distance(short_distri, long_distri[shift:len(long_distri)-shift-1])

            if(first > second):
                return first,int(diff/2)+1
            else:
                return second,int(diff/2)


def calc_distance(distriA, distriB):

    """ Take two distributions of equal length and return their similarity.
    Similarity is calculated as the Euclidean distance
    between the points in the distribution. """

    w = distriA.copy() #This is quicker than making an empty matrix
    middle = len(distriA)/2 #what if odd number?

    for i in range(len(distriA)):
        w[i] = (middle - abs(middle-i))+1 #Distance from center for i, but why +1 ?

    tot = w.sum()

    wn = [weight/tot for weight in w] #weight normalized?

    return weightedL2(distriA, distriB,wn) #Function defined earlier?


def profile_distance(seqA,seqB,prop='hydrophobicity',st = 0,l = 0):
    #Why is this seperate?
    #Create the two distributions with a numpy array, uses make_profile function
    # st = 0 stdev?
    # l = 0 ?
    distance, position = distribution_similarity(np.asarray(make_profile(seqA[start:len(seqA)-l], prop)),
                                                 np.asarray(make_profile(seqB[start:len(seqB)-l], prop)))
    return distance

def profile_distance_allprop(seqA,seqB,st = 0,l = 0):
    distanceALL = 0
    for prop in elektrochem:
        distance, position = distribution_similarity(np.asarray(make_profile(seqA[st:len(seqA)-l], prop)),
                                                     np.asarray(make_profile(seqB[st:len(seqB)-l], prop)))
        distanceALL = distanceALL + distance #Distance over all properties between two sequences?
    return distanceALL

def profile_distance_allprop_realign(seqA,seqB,st = 0,l = 0):
    distanceALL = 0
    for prop in elektrochem:
        distance, position = distribution_similarity_realign(np.asarray(make_profile(seqA[st:len(seqA)-l], prop)),
                                                     np.asarray(make_profile(seqB[st:len(seqB)-l], prop)))
        distanceALL = distanceALL + distance
    return distanceALL

if __name__ == '__main__':

    # Calculate distance between two CDR3 sequences as follows:
    CDR3_A = 'CASSLWTGSHEQYF'
    CDR3_B = 'CASSLWTGSHEQYF'
    #CDR3_B = 'CSARDRTGNGYTF'
    distance = profile_distance_allprop(CDR3_A,CDR3_B) #Distance is summed
    print(distance)
