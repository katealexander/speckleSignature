#!/usr/bin/python

## this is for three groups

import sys, re, os
from numpy import median
from scipy.stats import ttest_ind

def main(args):
    if len(args) != 4: sys.exit("USAGE: python getMedianExpressionOfPatientGroup.py hg19_genes.txt KIRC_specklePatientGroups.txt ~/Desktop/Alexander2023_filesTooBigForGithub/speckleSignature/KIRC_FPKM_UQ.txt > medianGeneExpression_KIRC_specklepatientGroups.txt")
    
    genes = []
    f = open(args[1])
    line = f.readline()[:-1]
    while line != "":
        genes.append(line)
        line = f.readline()[:-1]
    f.close

    IDdict = {}
    f = open(args[2])
    line = f.readline()[:-1]
    while line != "":
        items = line.split('\t')
        patientID = items[0]
        patientGroup = items[1]
        IDdict[patientID] = patientGroup
        line = f.readline()[:-1]
    f.close
    
    group1i = []
    group2i = []
    group3i = []
    f = open(args[3])
    line = f.readline()[:-1]
    
    # this will make a list of index numbers that each group has
    headerItems = line.split('\t')
    i = 0
    for x in headerItems:
        if x in IDdict.keys():
            if IDdict[x] == "1":
                group1i.append(i)
            elif IDdict[x] == "2":
                group2i.append(i)
            elif IDdict[x] == "3":
                group3i.append(i)
        i += 1
    
    header = "Ensembl\tGeneSymbol\tSigI\tNormal\tSigII\tpIvsN\tpIIvsN\tpIvsII"
    print header
    line = f.readline()[:-1]
    while line != "":
        items = line.split('\t')
        ensID = items[0]
        gene = items[1]
        if gene in genes:
            medianGroup1 = median([float(items[i]) for i in group1i])
            medianGroup2 = median([float(items[i]) for i in group2i])
            medianGroup3 = median([float(items[i]) for i in group3i])
            ratio1to3 = medianGroup1/medianGroup2
            tstat, p1vs2 = ttest_ind( [float(items[i]) for i in group1i], [float(items[i]) for i in group2i])
            tstat, p3vs2 = ttest_ind( [float(items[i]) for i in group3i], [float(items[i]) for i in group2i])
            tstat, p1vs3 = ttest_ind( [float(items[i]) for i in group1i], [float(items[i]) for i in group3i])
            if medianGroup1 > 0 and medianGroup2 > 0 and medianGroup3 > 0:
                print ensID + "\t" + gene + "\t" + str(medianGroup1) + "\t" + str(medianGroup2) + "\t" + str(medianGroup3) + "\t" + str(p1vs2) + "\t" + str(p3vs2) + "\t" + str(p1vs3)
        line = f.readline()[:-1]
    f.close
    
    
    

if __name__ == "__main__": main(sys.argv)
