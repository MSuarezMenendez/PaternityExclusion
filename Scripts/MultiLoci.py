

#! /usr/bin/env python3
from itertools import combinations
import operator
import sys
import multiprocessing
import types
import pandas as pd
import random
Microsatellites= "Micro_32_loci_GM_nodups.txt"
micro = pd.read_csv(Microsatellites, sep=" ")
Microheader = list(micro)

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
print(Locus)


for i in range(1,19):
    Dup = random.choice(list(Locus.values()))
    Dup1 = Dup + "(1)"
    Dup2 = Dup + "(2)"
    micro[Dup1 + str(i)] = micro.loc[:, Dup1]
    micro[Dup2 + str(i)] = micro.loc[:, Dup2]
micro.to_csv("GM_dupli50loci_from32.csv", index = False, sep=',')
