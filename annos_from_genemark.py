#!/usr/bin/python
# Callum Le Lay 01/03/2019 This script takes the output files from genemark.hmm (v3.25 - heuristics as in 1999 pub) (GFF), batch CD-Search and a fasta file of virus sequences and creates tables suited for littlegenomes
# Usage: make_anno_tables.py [GENEMARK.OUT] [GENOMES.FASTA] [OUTPUT_NAME]

import sys
import re
from random import randint

# Load tables
def load(file):
	print('Loading: {}'.format(file))
	f = open(file,'r')
	data = []
	for line in f.readlines():
		if not '#' == line[0]: 
			line = line.strip().split()
			if not line == []:
				data.append(line)
	f.close()
	return data
gm_t = load(sys.argv[1]) # Loads genemark output table data
gf_t = load(sys.argv[2]) # Loads genome fasta file data

# Extract relevant data
data = {}
gene_d = {}
i = 0
while i < len(gf_t):
	data[gf_t[i][0][1:]] = [len(gf_t[i+1][0])] 
	i += 2
for i in gm_t:
	gene_id = re.search('(\d+)',i[8]).group(0)
	gene_d[gene_id] = (i[3],i[4],i[6])
	data[i[0]].append((i[0],i[3],i[4],'0','0','white','','m')) # First 8 columns are directly for the resulting table - the 9th and 10th are for sorting out the start finish of domain hits

# Write annotation table
ids = data.keys()
f = open(sys.argv[3],'w')
for i in ids:
	j = 1
	while j < len(data[i]):
		anno = "\t".join(data[i][j])
        	f.write(anno+"\n")
		j += 1
f.close()


