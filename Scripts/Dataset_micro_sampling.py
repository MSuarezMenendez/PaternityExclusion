#! /usr/bin/env python3
from itertools import combinations
import operator
import sys
import multiprocessing
import types
import pandas as pd
import random
Microsatellites= "Mn_32loci_Paternity.txt"
micro = pd.read_csv(Microsatellites, sep="\t")
Microheader = list(micro)

Number_loci= 32
Repe = Number_loci - 1

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
Locus = list(Locus.values())

def subset(Loci, Set, rep):
    global Locus
    Dup = random.sample(Locus, Loci)
    print(Dup)
    Locus = [x for x in Locus if (x not in Dup)]
    print(Locus)
    for locus in Dup:
        Dup1 = locus + "(1)"
        Dup2 = locus + "(2)"
        df[Dup1] = micro[Dup1]
        df[Dup2] = micro[Dup2]
    df.to_csv(str(Set) + "_" + str(rep) +".csv", index = False, sep='\t')

for rep in range(1,51): #Replicates
    Loci = enumerate(Microheader)
    Locuslist = list()
    Locus = {}
    for number, locus in Loci:
        if number % 2 == 0 and number != 0:
            Locuslist.append(locus)
    for number, locus in enumerate(Locuslist):
        number = number + 1
        Locus[number] = locus[:-3]
    Locus = list(Locus.values())
    df = pd.DataFrame()
    df['Category'] = micro['Category']
    subset(8, 8, rep)
    for i in range(9,Repe):
        subset(1, i, rep)


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
Locus = list(Locus.values())



comb = combinations(Locus, Number_loci - 1)
counter = 0
for set in comb:
    counter += 1
    df = pd.DataFrame()
    df['Category'] = micro['Category']
    for locus in set:
        Dup1 = locus + "(1)"
        Dup2 = locus + "(2)"
        df[Dup1] = micro[Dup1]
        df[Dup2] = micro[Dup2]
    df.to_csv(str(Number_loci - 1) + "_" + str(counter) +".csv", index = False, sep='\t')


df = pd.DataFrame()
df['Category'] = micro['Category']
for locus in Locus:
    Dup1 = locus + "(1)"
    Dup2 = locus + "(2)"
    df[Dup1] = micro[Dup1]
    df[Dup2] = micro[Dup2]
df.to_csv(str(Number_loci) + "_" + str(1) +".csv", index = False, sep='\t')

