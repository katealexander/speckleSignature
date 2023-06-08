#!/usr/bin/python

import sys, re, numpy, os

def main(args):
	if len(args) <= 3: sys.exit("USAGE: python makeExpressionTable_GDC.py gdc_sample_sheet_KIRC_VHLmut.2021-09-20.tsv ENS_geneName.txt pathToExpression/*/*UQ.txt > KIRC_VHLmutant_FPKM_UQ.txt")

	### Get sample ID link from sample sheet
	fileNameDict = {}
	f = open(args[1])
	line = f.readline()[:-1]
	while line != "":
		items = line.split('\t')
		fileName = items[1] 
		sampleID = items[6]
		fileNameDict[fileName] = sampleID
		line = f.readline()[:-1]
	f.close
	
	### Create dictionary linking ensemble id to gene name (downloaded from UCSC table browser)
	f = open(args[2])
	geneNameDict = {}
	line = f.readline()[:-1]
	while line != "":
		items = line.split('\t')
		eID = items[1]
		name = items[3]
		geneNameDict[eID] = name
		line = f.readline()[:-1]
	f.close
		
	header = ["ENS_ID", "GENE_SYM"]
	expressionDict = {}
	for i in range(3,len(args)):
		## add gene symbol to header
		header.append(fileNameDict[os.path.basename(args[i]) + '.gz'])
		f = open(args[i])
		line = f.readline()[:-1]
		while line != "":
			ensID = line.split('\t')[0].split('.')[0]
			exp = line.split('\t')[1]
			if ensID in expressionDict.keys():	
				## append expression of sample to the dictionary
				expressionDict[ensID].append(exp)
			else:
				if ensID in geneNameDict.keys():
					expressionDict[ensID] = [geneNameDict[ensID], exp]
			line = f.readline()[:-1]
		f.close
		
	print '\t'.join(header)
	for key in expressionDict.keys():
		print key + '\t' + '\t'.join(str(i) for i in expressionDict[key])
	

	
if __name__ == "__main__": main(sys.argv)
