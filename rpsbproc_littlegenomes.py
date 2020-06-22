#!/usr/bin/python
# Callum Le Lay This script converts the output of rpsbproc (which reformats rps-blast output to similar to CD-Search) to littlegenomes annotation (and genome) table

# Load in libraries
import re
import argparse

# Take arguments
parser = argparse.ArgumentParser(description="Convert rpsbproc output to littlegenomes input.")
parser.add_argument("input_file",help="rpsbproc output to be converted")
parser.add_argument("-g","--genome_file", help="genome data will be retrieved, formatted for littlegenomes and stored here",default=None)
parser.add_argument("-o","--anno_file",help='destination for littlegenomes format annotations table, STDOUT by default')
parser.add_argument("-s","--separator",help='character used to separate fields in output, default is tab',default='\t')
args = parser.parse_args()

# Containers for annotation and genome data
names = {}
annos = []

# Load in the rpsbproc output file
f = open(args.input_file,'r')
lines = [i.strip('\n').split('\t') for i in f.readlines()]
f.close()

# Derive data
# Parse through lines
qid=None
in_domains = False
for line in lines:
	if len(line) < 2:
		in_domains = True if line[0] == "DOMAINS" else in_domains
		in_domains = False if line[0] == "ENDDOMAINS" else in_domains
	# When a 'QUERY' line is hit, add name to dict with query id as key
	elif line[0] == "QUERY":
                qid = line[1] # Keep current query id
		names[qid] = (line[3],line[4])
        # Check if line corresponds to current query id
	elif in_domains and re.match(qid, line[1]):
		e = float(line[6])
		if e < 0.005: 	
        		# Store data for annotation, in littlegenomes format
			name = names[qid][1]
			strt = line[4]
			end = line[5]
			label = line[9]+' ('+line[8]+')' 
			annos.append((name,strt,end,'NA','NA','NA',label,'m'))

# Output annotation table
for line in annos:
	print(args.separator.join(line))	

# Output genome table
if args.genome_file:
	f = open(args.genome_file,'w')
	for key in names.keys():
		f.write(names[key][1]+'\t'+names[key][0]+'\n')
	f.close()
