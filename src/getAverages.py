import h5py
import statistics
import itertools
# import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
import argparse, sys

def get_dict_of_permutations(kmerLength, alphabet = "ACTG"):
    seperator = ''
    perm = list(itertools.product(alphabet, repeat=5))
    # return dict.fromkeys(map(lambda x: seperator.join(x), perm),[])
    return dict.fromkeys(map(lambda x: seperator.join(x), perm),0)

def make_avg_signal_dict_from_data(attrs, events, data, kmerLength = 5):
    dict = get_dict_of_permutations(kmerLength)
    nuc_dict = {"A":[], "C":[], "T":[], "G":[]} 
    start = attrs['mapped_start']
    # end = attrs['mapped_end']
    for eventIndex in range(0, len(events)-5):
        kmer = ""
        avg = 0
        for i in range (kmerLength):
            event = events[eventIndex+i]
            eventData = data[start+event[2]:start+event[2]+event[3]]
            avgForEvent = statistics.mean(eventData)
            kmer = kmer+str(event[4])[2:3]
            #avg = avg + data[event[2]]
            avg = avg + avgForEvent
        avg = avg / kmerLength
        # dict[kmer].append(avg)
        if dict[kmer] == 0:
            dict[kmer] = avg
        else:
            dict[kmer] = (dict[kmer] + avg)/2
        
        nuc_legnth = events[eventIndex][3]
        nuc = str(events[eventIndex][4])[2:3]
        nuc_dict[nuc].append(nuc_legnth)
        # if nuc_dict[nuc] == 0:
        #     nuc_dict[nuc] = nuc_legnth
        # else:
        #     nuc_dict[nuc] = (nuc_dict[nuc] + nuc_legnth)/2
    # for key in dict: 
    #     if len(dict[key]) == 0:
    #         dict[key] = 0
    #     else:
    #         dict[key] = statistics.mean(dict[key])
    return (dict,nuc_dict)

parser=argparse.ArgumentParser()

# parser.add_argument('--fasta-filename', help='Path to fasta file')
parser.add_argument('--fast5-filename', help='Path to fast5 file')
# parser.add_argument('--plot', help='Whether to plot the results or not')
args=parser.parse_args()
print(args)

#FAST5
filename = 'db6b45aa-5d21-45cf-a435-05fb8f12e839.fast5'
fast5_filename = args.fast5_filename
if fast5_filename is not None:
    print(fast5_filename)
    filename = fast5_filename
f = h5py.File(filename, 'r')
keys = list(f.keys())

#TOMBO
analyses =  f[keys[0]]
tombo_out = analyses[list(analyses.keys())[1]]
baseCalledTemplate = tombo_out[list(tombo_out.keys())[0]]
alignment = baseCalledTemplate[list(baseCalledTemplate.keys())[0]]
alignmentAttrs = dict(alignment.attrs.items())
events = baseCalledTemplate[list(baseCalledTemplate.keys())[1]]
events = list(events)

#RAW
raw = f[keys[1]]
reads = raw[list(raw.keys())[0]]
read = reads[list(reads.keys())[0]]
signal = read[list(read.keys())[0]]

avgsDict, nucDict = make_avg_signal_dict_from_data(alignmentAttrs, events, signal)
# print(nucDict)
for key in nucDict:
    std = np.std(nucDict[key])
    mean = np.mean(nucDict[key])
    median = np.median(nucDict[key])
    print((key, mean, median, std))
# plt.bar(range(len(avgsDict)), list(avgsDict.values()), align='center')
# plt.xticks(range(len(avgsDict)), list(avgsDict.keys()))
# plt.show()
