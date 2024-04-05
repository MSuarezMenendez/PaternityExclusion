#! /usr/bin/env python3
from itertools import combinations
import operator
import sys
import multiprocessing
import types
import pandas as pd
import random
Microsatellites= "./Microsatellites.tsv"
micro = pd.read_csv(Microsatellites, sep="\t")
Microheader = list(micro)

Extra = 20 #Number of extra loci to be added

#List of loci
Loci = enumerate(Microheader)
Locuslist = list()
Locus = {}
for number, locus in Loci:
    if number % 2 == 0 and number != 0:
        Locuslist.append(locus)
for number, locus in enumerate(Locuslist):
    number = number + 1
    Locus[number] = locus[:-3]


for i in range(1, Extra + 1): 
    Dup = random.choice(list(Locus.values()))
    Dup1 = Dup + "(1)"
    Dup2 = Dup + "(2)"
    micro[Dup1 + str(i)] = micro.loc[:, Dup1]
    micro[Dup2 + str(i)] = micro.loc[:, Dup2]
micro.to_csv("Extra_loci.tsv", index = False, sep='\t')
