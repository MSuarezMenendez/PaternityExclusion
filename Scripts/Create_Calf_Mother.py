#! /usr/bin/env python3
import csv
import itertools


with open("Samples.txt") as f:
    po = f.readlines()
po = [x.strip().split('\t')[0] for x in po]
Sex_dict = {}
with open("Sexes.txt") as f:
    Sex = f.readlines()
    Sex = [x.strip().split('\t') for x in Sex]
    for x in Sex:
        Sex_dict[x[0]] = x[1]

Prueba= itertools.combinations(po, 2)
for line in Prueba:
    indv1 = line [0]
    indv2 = line[1]
    Sex1 = Sex_dict.get(indv1)
    Sex2= Sex_dict.get(indv2)
    if Sex1 == "male" and Sex2 == "male":
        continue
    if Sex1 == "female" and Sex2 == "male":
        print(indv2,indv1, sep='\t')
        continue
    if Sex1 == "male" and Sex2 == "female":
        print(indv1,indv2, sep='\t')
        continue
    if Sex1 == "female" and Sex2 == "female":
        print(indv1,indv2, sep='\t')
        print(indv2,indv1, sep='\t')
        continue
