
import random
import argparse
from argparse import ArgumentParser
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
import seaborn as sns
from kneed import KneeLocator
from scipy.signal import savgol_filter #Smothing

parser = ArgumentParser()
flag = parser.add_argument_group('Arguments')
flag.add_argument("-r", action="store", dest="Repetitions", help="Repetitions")
flag.add_argument("-i", action="store", dest="Input", help="Input")
flag.add_argument("-o", action="store", dest="Output", help="Output")
args = parser.parse_args()
#Num_loci = int(args.Num_loci)
Repetitions = int(args.Repetitions)
Input = args.Input
Output = args.Output

with open("Males.txt") as f: #List of male sample IDs
    Males = f.readlines()
Males = [x.strip().split('\t') for x in Males]

with open("Females.txt") as f: #List of male sample IDs
    Females = f.readlines()
Females = [x.strip().split('\t') for x in Females]

with open(Input) as f:
    micro = f.readlines()
micro = [x.strip() for x in micro]

def counts(lst, x):
    return lst.count(x)

def mergealleles(info):
    lista = list() #Merge alleles from same locus
    micros = list()
    count = 0
    for allele in info:
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

Loci = enumerate(micro[0].split(" "))
Locuslist = list()
for number, locus in Loci:
    if number % 2 == 0 and number != 0:
        Locuslist.append(locus[:-1])

#print("Number of loci: {}".format(Num_loci))
print("Replicates: {}".format(Repetitions))
#Halfsib of calf (father's side)
Repetition = list()
Represult_all = list()
Represult_1ms = list()
Represult_2ms = list()
Represult_3ms = list()
Locinum = list()
for Num_loci in range(8,51):
    print("Loci number: {}".format(Num_loci))
    for i in range(1,51):
        RepID = "Rep{}".format(i)
        print(RepID)
        Result = list()
        for i in range(Repetitions):
            #Father
            male = random.choice(Males)
            micromale = list(filter(lambda f:male[0] in f, micro))
            micromale = micromale[0].split(" ")[1::]
            allelesmale = mergealleles(micromale)
            #Calf
            Fatherinherit = list()
            for Locusmale in allelesmale:
                Maleallele = random.choice(Locusmale)
                Fatherinherit.append(Maleallele)
            #HalfSib
            female = random.choice(Females)
            microfemale = list(filter(lambda f:female[0] in f, micro))
            microfemale = microfemale[0].split(" ")[1::]
            allelesfemale = mergealleles(microfemale)
            Halfsiballeles = list()
            for Locusmale, Locusfemale in zip(allelesmale, allelesfemale):
                Locus = list()
                Calfallele = random.choice(Locusmale)
                Parent2allele = random.choice(Locusfemale)
                Locus.append(Calfallele)
                Locus.append(Parent2allele)
                Halfsiballeles.append(Locus)
            Sex = ["Male", "Female"] # Sex of calf of calf
            Matchnum = 0
            if random.choice(Sex) == "Male":
                for Locusfather, Locushalfsib in zip(Fatherinherit, Halfsiballeles):
                    if Locusfather in Locushalfsib and Locusfather != "-9":
                        Matchnum += 1
            Result.append(Matchnum)
        All = counts(Result,Num_loci) * 100 / Repetitions
        Unms = counts(Result, Num_loci - 1) * 100 / Repetitions
        Dosms = counts(Result, Num_loci - 2) * 100 / Repetitions
        #Tresms = counts(Result, Num_loci - 3) * 100 / Repetitions
        Represult_all.append(All)
        Represult_1ms.append(Unms)
        Represult_2ms.append(Dosms)
        #Represult_3ms.append(Tresms)
        Repetition.append(RepID)
        Locinum.append(Num_loci)
df = pd.DataFrame(list(zip(Locinum,Represult_all)), columns =['Loci','FP(%)'])
#sns.relplot(data=df, x = "Loci",y = "FP(%)", ci=95, estimator=np.median, kind="line")
#plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}')) # No decimal places
#plt.savefig("myImagePDF.pdf", format="pdf")

df1 = pd.DataFrame(list(zip(Locinum,Represult_1ms)), columns =['Loci','FP(%)'])
df2 = pd.DataFrame(list(zip(Locinum,Represult_2ms)), columns =['Loci','FP(%)'])


X = list()
Y = list()
for Num_loci in range(8,51):
    Temp = df.loc[df['Loci'] == Num_loci]
    X.append(Temp['Loci'].median())
    Y.append(Temp['FP(%)'].median())
yhat = savgol_filter(Y, 26, 4)
Finalnom = pd.DataFrame(list(zip(X,yhat)), columns =['Loci','FP(%)'])
kn = KneeLocator(X, yhat, curve='convex', direction='decreasing', online = True, S=1)
print(kn.knee)



X = list()
Y = list()
for Num_loci in range(8,51):
    Temp = df1.loc[df1['Loci'] == Num_loci]
    X.append(Temp['Loci'].median())
    Y.append(Temp['FP(%)'].median())
yhat = savgol_filter(Y, 26, 4)
Final1m = pd.DataFrame(list(zip(X,yhat)), columns =['Loci','FP(%)'])
kn = KneeLocator(X, yhat, curve='convex', direction='decreasing', online = True, S=1)
print(kn.knee)




X = list()
Y = list()
for Num_loci in range(8,51):
    Temp = df2.loc[df2['Loci'] == Num_loci]
    X.append(Temp['Loci'].median())
    Y.append(Temp['FP(%)'].median())
yhat = savgol_filter(Y, 26, 4)
Final2m = pd.DataFrame(list(zip(X,yhat)), columns =['Loci','FP(%)'])
kn = KneeLocator(X, yhat, curve='convex', direction='decreasing', online = True, S=1)
print(kn.knee)




concatenated = pd.concat([Finalnom.assign(Mismatches='No mismatch'), Final1m.assign(Mismatches='One mismatch'), Final2m.assign(Mismatches='Two mismatches')],ignore_index=True)
print(concatenated)
concatenated.to_csv("Result_simu.txt", sep='\t')
ax = sns.relplot(x = "Loci",y = "FP(%)", hue = 'Mismatches',kind ="line", data=concatenated)
sns.move_legend(ax, "center right")
plt.ylim(0, 18)
x = np.arange(8, 56, 6)
plt.xticks(x)
plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}')) # No decimal places
plt.savefig(Output, format="pdf")
