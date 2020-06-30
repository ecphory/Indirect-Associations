#!/usr/bin/env python
# -*- coding: utf-8 -*-
# some code written by Saranya Venkatraman
# the rest written by M. A. Kelly
import numpy as np
import math
from numpy.fft import fft,ifft

"""
Function to load pre-trained GloVe embeddings
Input : Path to pre-trained GloVe embeddings file
Returns: A dictionary mapping word -> embedding

"""

N = 1024 # the dimensions
embeddings_file_path = ['L1W5v2order.csv',
                        'L2W5v2order.csv',
                        'L3W5v2order.csv',
                        'L4W5v2order.csv']
pos_file_path = 'words.txt'
out_path = 'L1toL4W5.csv'
shortNames = ['L1','L2','L3','L4']
thresholdFreq = 5

def loadBeagleModel(vectorFile):
    print("Loading BEAGLE Model from " + vectorFile)
    f = open(vectorFile,'r')
    model = {}
    for line in f:
        splitLine = line.strip().split(",")
        word   = splitLine[0].strip()
        vector = splitLine[1:(N+1)]
        embedding = np.array([float(val.strip()) for val in vector])
        model[word] = embedding
    print("Done.",len(model)," words loaded!")
    f.close()
    return model

def magnitude(a):
    return math.sqrt(np.dot(a,a))

def normalize(a):
    return a / magnitude(a)

def cosine(a,b):
    return np.dot(a,b) / (magnitude(a) * magnitude(b))

def cconv(a,b):
    return ifft(fft(a) * fft(b)).real

count = [0 for i in range(len(embeddings_file_path))]
score = [0 for i in range(len(embeddings_file_path))]

file_num = 0
for file_name in embeddings_file_path:
    embeddings = loadBeagleModel(file_name)
    word_list = list(embeddings.keys()) # BEAGLE's vocabulary
    print("Finished loading BEAGLE Model")

    freq_list  = {}
    index_list = {}
    pos_list   = {}
    g = open(pos_file_path,'r')
    for line in g:
        splitLine = line.strip().split(" ")
        itemList  = []
        for item in splitLine:
            if item:
                itemList = itemList + [item]
        word  = itemList[0]
        index = int(itemList[1])
        freq  = int(itemList[2])
        pos   = ''.join(sorted(itemList[3]))
        if freq >= thresholdFreq:
            index_list[word] = index
            freq_list[word]  = freq
            pos_list[word]   = pos
    g.close()
    print("Finished loading Part of Speech and Word Frequencies")

    pos_vec  = {}
    pos_count= {}
    for word in pos_list:
        vector  = embeddings[word]
        freq    = freq_list[word]
        pos     = pos_list[word]
        if pos in pos_vec:
            prototype= pos_vec[pos]
        else:
            prototype = np.zeros(N)
        if thresholdFreq > 0:
            prototype = prototype + vector
        else:
            prototype = prototype + np.log(freq)*vector
        pos_vec[pos] = prototype
        if pos in pos_count:
            if freq > pos_count[pos]:
                pos_count[pos] = freq
        else:
            pos_count[pos] = freq

    numPOS = len(pos_vec)

    print("Finished constructing POS vectors")
    print("Detected "+ str(numPOS) + " unique part of speech tag combinations")
    for pos in pos_vec:
        print(pos + " " + str(pos_count[pos]))


    # count all instances where the closes POS tag is the right one
    for word in pos_list:
        wordVector = embeddings[word]
        wordPOS    = pos_list[word]
        wordFreq   = freq_list[word]
        cos_list   = {}
        if wordFreq < thresholdFreq:
            print("Error: " + word + " below threshold frequency")
        elif magnitude(wordVector) == 0:
            print("Empty vector for word '" + word + "'. Skipping to next word...")
        else:
            for pos, vec in pos_vec.items():
                if pos_count[pos] < thresholdFreq:
                    print("Error: " + pos + " below threshold frequency")
                else:
                    if magnitude(vec) == 0:
                        print("Empty vector for " + pos + " at threshold " + str(thresholdFreq))
                        cos_list[pos] = 0
                    else:
                        cos_list[pos] = cosine(wordVector,vec)
            # find the maximum cosine for each frequency threshold
            pos_tag = max(cos_list, key=lambda key: cos_list[key])
            if pos_tag == wordPOS:
                score[file_num] = score[file_num] + 1
            count[file_num] = count[file_num] + 1
    file_num = file_num + 1

print(score)
print(count)

h = open(out_path, "w")
for short in shortNames:
    h.write("%s," % short)
h.write("\n")
for value in score:
    h.write("%s," % value)
h.write("\n")
for value in count:
    h.write("%s," % value)
h.write("\n")
for numerator, denominator in zip(score,count):
    percent = numerator/denominator
    h.write("%s," % percent)
h.close()
