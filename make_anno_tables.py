#!/usr/bin/python
# Callum Le Lay 01/03/2019 This script takes the output files from genemark.hmm (v3.25 - heuristics as in 1999 pub) (GFF), batch CD-Search and a fasta file of virus sequences and creates tables suited for littlegenomes
# Usage: make_anno_tables.py [GENEMARK.OUT] [CDSEARCH.OUT] [GENOMES.FASTA] [OUTPUT_NAME]

import sys
import re
from random import randint

# Load tables
def load(file):
	print 'Loading: {}'.format(file)
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
cs_t = load(sys.argv[2]) # Loads cd-search output table data
cs_t.pop(0) # Removes header
gf_t = load(sys.argv[3]) # Loads genome fasta file data

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
for i in cs_t:
	m = re.match('(>gene_\d*_)(.*)',i[2])
	id = m.group(2)
	gene_id = re.search('\d+',m.group(1)).group(0)
	data[id].append((id,i[5],i[6],'0','1','red',i[10],'s',gene_id))		

# Check for input errors

# Make set of domain names and assign colours 
domains = []
for i in data.keys():
	[domains.append(j[6]) for j in data[i][1:]] 
	data[i] = [data[i][0]]+list(set(data[i][1:])) # Removes duplicates
domains = list(set(domains))
dom2col = {}
cols = ('orange','blue','red','green','light green','brown','pink','yellow')
j = 0
for i in domains:
	dom2col[i] = cols[j]
	j = 0 if j == len(cols)-1 else j+1

# Convert index values of domain hits to relative to genome nucleotides
for i in data.keys():
	j = 1
	while j<len(data[i]):
		if data[i][j][7] == 's':
			strt = int(data[i][j][1])
			end = int(data[i][j][2])
			gid = data[i][j][8]
			if gene_d[gid][2] == '+':
				nstrt = strt*3+int(gene_d[gid][0])-3
				nend = end*3+int(gene_d[gid][0])
			else:
				nend = int(gene_d[gid][1])-3*strt+3
                                nstrt = int(gene_d[gid][1])-3*end
			data[i][j] = list(data[i][j])
			data[i][j][1] = str(nstrt)
			data[i][j][2] = str(nend)
			data[i][j][5] = str(dom2col[data[i][j][6]])
			data[i][j].pop(8)
		j += 1
			
# Print out tables
# Write genome info table
ids = data.keys()
f = open(sys.argv[4]+"_genome_info.txt",'w')
for i in ids:
	f.write("{0}\t{1}\n".format(i,data[i][0]))
f.close()
# Write annotation table
f = open(sys.argv[4]+"_annos.txt",'w')
for i in ids:
	j = 1
	while j < len(data[i]):
		anno = "\t".join(data[i][j])
        	f.write(anno+"\n")
		j += 1
f.close()


