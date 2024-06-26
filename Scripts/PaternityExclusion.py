#! /usr/bin/env python3
from itertools import combinations
import operator
import sys
import multiprocessing
import types
import pandas as pd

from argparse import ArgumentParser
parser = ArgumentParser()
flag = parser.add_argument_group('Arguments')
flag.add_argument("-m", action="store", dest="MinimumML", help="Match threshold")
flag.add_argument("-i", action="store", dest="Microsatellites", help="Input file with microsatellite data")
flag.add_argument("-d", action="store", dest="Delimiter", help="Microsatellite\
 file column delimiter", default ="\t")
flag.add_argument("-s", action="store", dest="Sires", help="Text file with list of males (potential sires)")
flag.add_argument("-c", action="store", dest="CalfCow", help="Text file with potential calf-cow pairs")
flag.add_argument("-t", action="store", dest="Cores", help="Number of cores to use", default = 1)
if len(sys.argv) < 2: #If no arguments are provided, help is printed
	sys.stderr.write("Paternity exclusion analysis from microsatellite markers.\n\
Marine Evolution and Conservation\
\nUniversity of Groningen\nGitHub: https://github.com/MSuarezMenendez/PaternityExclusion\nContact: marcos.sume@gmail.com\n\n")
	parser.print_help()
	sys.exit(1)

args = parser.parse_args()
MinimumML = int(args.MinimumML)
Microsatellites = args.Microsatellites
Sires = args.Sires
Delimiter = args.Delimiter
CalfCow = args.CalfCow
Cores = int(args.Cores)


def missingdata(allele):
    if allele[0] == "-9":
        allele0 = "NA"
    else:
        allele0 = allele[0]
    if allele[1] == "-9":
        allele1 = "NA"
    else:
        allele1 = allele[1]
    allelef = [allele0, allele1]
    return allelef


with open(CalfCow) as f: # Calf mother pairs file. Calfs, first column. Mothers, second column. Tab separated.
    content = f.readlines()
content = [x.strip().split('\t') for x in content]
Todo = content
with open(Sires) as f: #List of male sample IDs
    Males = f.readlines()
Males = [x.strip().split('\t') for x in Males]

micro = pd.read_csv(Microsatellites, sep=Delimiter, engine='python')
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
Num_loci = len(Locuslist)
def mergealleles(info):
    lista = list() #Merge alleles from same locus
    micros = list()
    count = 0
    for allele in info:
        allele = str(allele)
        count += 1
        if count <= Num_loci*2:
            if len(micros) < 2:
                micros.append(allele)
            else:
                lista.append(micros)
                micros = list()
                micros.append(allele)
            if count == Num_loci * 2 -1:
                temp = list()
                temp.append(allele)
            if count == Num_loci * 2:
                temp.append(allele)
                lista.append(temp)
        else:
            continue
    return lista

def findinlist(what,lista):
    list(filter(lambda f:what in f, lista))

print("Calf\tDam\tDam-Calf_MS\tDam-Calf_MSL\tDam-Calf_MM\tDam-Calf_MML\tSire\tMatching\tSire_MS\tSire_MSL\tMismatch\tMismatch_Loci")


def search(Parejas):
    Allelemis = {}
    Mismatchloci = {}
    Resultados = list()
    microcalf = micro.loc[micro['IDs'] == Parejas[0]].values.tolist()[0][1::]
    micromother = micro.loc[micro['IDs'] == Parejas[1]].values.tolist()[0][1::]
    allelesmother = mergealleles(micromother)
    allelescalf = mergealleles(microcalf)
    Sirecheck = list()
    Counter = 0
    for allelem, allelec, in zip(allelesmother,allelescalf): # Check the behaviour of this changing test data
        Counter += 1
        Mismatch = 0
        Missing = 0
        allelem = missingdata(allelem)
        if allelec[0] in allelem and allelec[1] in allelem: # Scenarios: Mother 153 153 Calf 153 153, Mother 153 175 Calf 153 175,  Mother 153 153 Calf 153 173
            templist = list()
            templist.append(allelec[0])
            templist.append(allelec[1])
            Sirecheck.append(templist)
            continue
        if allelec[0] not in allelem and allelec[1] not in allelem:
                if "NA" not in allelem: #Scenario: Mother 153 171 Calf 150 -9
                    if allelec[0] != "-9" and allelec[1] == "-9":
                        Sirecheck.append(allelec[0])
                        Allelemis[Counter] = Missing
                        continue
                    if allelec[1] != "-9" and allelec[0] == "-9":
                        Sirecheck.append(allelec[1])
                        Allelemis[Counter] = Missing
                        continue
        if allelec[0] in allelem and "NA" not in allelem: #Scenarion: Mother 153 -9 Calf 153 171, Mother 153 187 Calf 153 -9
                if allelec[1] != "-9":
                    Sirecheck.append(allelec[1])
                    Allelemis[Counter] = Missing
                    continue
                else:
                    Sirecheck.append("NA")
                    Missing+=1
                    Allelemis[Counter] = Missing
                    continue
        elif allelec[1] in allelem and "NA" not in allelem: #Scenarion: Mother 171 -9 Calf 153 171
                if allelec[0] != "-9":
                    Sirecheck.append(allelec[0])
                    Allelemis[Counter] = Missing
                    continue
                else:
                    Sirecheck.append("NA")
                    Missing+=1
                    Allelemis[Counter] = Missing
                    continue
        elif "NA" not in allelem and "-9" not in allelec:
            Sirecheck.append("NA")
            Mismatch += 1
            Mismatchloci[Counter] = Mismatch
            continue
        else:
            Sirecheck.append("NA")
            Missing+=1
            Allelemis[Counter] = Missing
            continue
    if Missing > Num_loci - MinimumML or Mismatch > Num_loci - MinimumML:
        return
    Locusmissing = list()
    for key, value in Allelemis.items():
        if value == 1:
            Locusmissing.append(Locus[key])
    Locusmismatch = list()
    for key, value in Mismatchloci.items():
        if value == 1:
            Locusmismatch.append(Locus[key])
    Potentials = list()
    Locusinfo ={}
    Malemismatchloci = {}
    Malemissingloci = {}
    Missing = {}
    Mismatch = {}
    for male in Males:
        if male[0] != Parejas[0]:
            count = 0
            micromale = micro.loc[micro['IDs'] == male[0]].values.tolist()[0][1::]
            allelesmale = mergealleles(micromale)
            Locusn = 0
            Missalelle = 0
            Mismatchmale = 0
            Locilistmissing = list()
            Locilistmismatch = list()
            for malleles, calfallele in zip(allelesmale,Sirecheck):
                Locusn += 1
                if isinstance(calfallele, list):
                    check = 0
                    for alelo in calfallele:
                        check += 1
                        if alelo in malleles:
                            count += 1
                            break
                        else:
                            if check == 2:
                                if "-9" in malleles:
                                    Missalelle += 1
                                    Locilistmissing.append(Locus[Locusn])
                                else:
                                    Mismatchmale += 1
                                    Locilistmismatch.append(Locus[Locusn])
                else:
                    if calfallele in malleles:
                        count += 1
                    else:
                        if "-9" in malleles:
                            Locilistmissing.append(Locus[Locusn])
                            Missalelle += 1
                        elif "NA" not in calfallele:
                            Mismatchmale += 1
                            Locilistmismatch.append(Locus[Locusn])
            if count >= MinimumML:
                Missing[male[0]] = Missalelle
                Mismatch[male[0]] = Mismatchmale
                Malemismatchloci[male[0]] = str(" ".join(Locilistmismatch))
                Malemissingloci[male[0]] = str(" ".join(Locilistmissing))
                Potentials.append(male)
                Locusinfo[male[0]] = count
                if len(Locusmissing) != 0:
                    Locmissing = " ".join(Locusmissing)
                else:
                    Locmissing = "None"
                if len(Locusmismatch) != 0:
                    Locmissm = " ".join(Locusmismatch)
                else:
                    Locmissm = "None"
                if Missing[male[0]] != 0:
                    Locmissing_M = Malemissingloci[male[0]]
                else:
                    Locmissing_M = "None"
                if Mismatch[male[0]] != 0:
                    Locmissm_M = Malemismatchloci[male[0]]
                else:
                    Locmissm_M = "None"
                print(Parejas[0],Parejas[1],len(Locusmissing),Locmissing,len(Locusmismatch),Locmissm,male[0],Locusinfo.get(male[0]),Missing[male[0]],Locmissing_M,Mismatch[male[0]],Locmissm_M, sep="\t", flush=True)


with multiprocessing.Pool(Cores) as pool:
    pool.map(search, Todo)
