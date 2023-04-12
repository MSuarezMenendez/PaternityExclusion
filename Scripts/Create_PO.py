import csv
import itertools


with open("Samples") as f:
    po = f.readlines()
po = [x.strip().split('\t')[0] for x in po]
Sex_dict = {}


Prueba= itertools.combinations(po, 2)
for line in Prueba:
    indv1 = line [0]
    indv2 = line[1]
    print(indv1, indv2, sep='\t')
