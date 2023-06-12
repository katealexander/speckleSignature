#!/usr/bin/python

## this is for three groups

import sys, re, os
from numpy import median

def main(args):
    if len(args) != 4: sys.exit("USAGE: python getExpressionForPatientGroupsOfGenes.py KIRC_VHLmutant_specklePatientGroups.txt TCGA-KIRC.STAR_NormalizedCounts_VHLlof.txt geneName > outfile")
    
    geneName = args[3]
    
    IDdict = {}
    f = open(args[1])
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
    f = open(args[2])
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
    
    line = f.readline()[:-1]
    line = f.readline()[:-1]
    while line != "":
        items = line.split('\t')
        gene = items[1]
        if gene == geneName:
            group1 = [items[i] for i in group1i]
            group2 = [items[i] for i in group2i]
            group3 = [items[i] for i in group3i]
            for sample in group1:
                print str(sample) + '\tgroup1'
            for sample in group2:
                print str(sample) + '\tgroup2'
            for sample in group3:
                print str(sample) + '\tgroup3'
                
            break
        line = f.readline()[:-1]
    f.close
    
    
    
if __name__ == "__main__": main(sys.argv)
