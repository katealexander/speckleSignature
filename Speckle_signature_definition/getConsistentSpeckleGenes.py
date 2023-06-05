#!/usr/bin/python

import sys, re

def main(args):
    if len(args) != 2: sys.exit("USAGE: python getConsitentSpeckleGenes.py speckleProteinGenes_rotations_allStudiesWithRNA_noNA.txt > outfile")
    
    ## these are the cancer types with consistent rotations (22 cancer types)
    cancerTypes = ['blca_tcga_pan_can_atlas_2018', 'cesc_tcga_pan_can_atlas_2018', 'esca_tcga_pan_can_atlas_2018', 'hnsc_tcga_pan_can_atlas_2018', 'kich_tcga_pan_can_atlas_2018', 'kirc_tcga_pan_can_atlas_2018', 'kirp_tcga_pan_can_atlas_2018', 'lgg_tcga_pan_can_atlas_2018', 'lihc_tcga_pan_can_atlas_2018', 'meso_tcga_pan_can_atlas_2018', 'pcpg_tcga_pan_can_atlas_2018', 'paad_tcga_pan_can_atlas_2018', 'prad_tcga_pan_can_atlas_2018', 'sarc_tcga_pan_can_atlas_2018', 'skcm_tcga_pan_can_atlas_2018', 'laml_tcga_pan_can_atlas_2018', 'stad_tcga_pan_can_atlas_2018', 'thca_tcga_pan_can_atlas_2018', 'thym_tcga_pan_can_atlas_2018', 'ucec_tcga_pan_can_atlas_2018', 'ucs_tcga_pan_can_atlas_2018', 'brca_tcga_pan_can_atlas_2018']
    
    ## count up how many of cancer types have positive (Sig I) and negative (Sig II) rotations
    geneDict = {}
    f = open(args[1])
    line = f.readline()[:-1]
    header = line.split("\t")
    indexList = [header.index(i) for i in cancerTypes]
    line = f.readline()[:-1]
    while line != "":
        line = line.strip()
        items = line.split("\t")
        gene = items[1]
        geneDict[gene] = [0,0] ## number of positive rotations and number of negative rotations
        for i in indexList:
            if float(items[i]) > 0:
                geneDict[gene][0] += 1
            elif float(items[i]) < 0:
                geneDict[gene][1] += 1
        line = f.readline()[:-1]
    print "Gene\tSigI\tSigII"
    for gene in geneDict.keys():
        print gene + "\t" + "\t".join(str(x) for x in geneDict[gene])
        
        
if __name__ == "__main__": main(sys.argv)
