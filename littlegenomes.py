#!/usr/bin/python3

# Callum Le Lay 210811
# Usage: littlegenomes.py GENOMES_TSV ANNOTATIONS_TSV OUPUT_NAME

# Import
import svgwrite
import sys
from math import ceil
import argparse

# Argument parsing
parser = argparse.ArgumentParser(description="Takes tables with descriptions of small genomes and their annotations and produces diagramtic representations of these viruses (.svg).")
parser.add_argument("geno_f", metavar="<GENOME_TABLE>", help="tab-delimited table with genome names corresponding to annotation table and sizes (in bp/nt) of the genomes")
parser.add_argument("anno_f", metavar="<ANNOTATION_TABLE>", help="tab-delimited annotation table in littlegenomes format")
parser.add_argument("out_f", metavar="<OUTPUT_FILE>", help=".svg output filename")
parser.add_argument("--squeeze", choices=['x','y',"both"], default=None, help="genomes are scaled to fit in a A4 sized page")
parser.add_argument("--multi", default=False, action="store_true", help="genomes are grouped by size and displayed with individual scaling and scalebars to aid readability in the case of large numbers of genomes")
parser.add_argument("--max_height",default=950,type=float,help="max size of figure. If y-axis squeeze is used, this will scale all genomes to fit within this height - defaults to A4 height")
parser.add_argument("--incre_height",default=100,type=int, help="space between individual sequences on the y-axis")
parser.add_argument("--offset_shift",default=False,choices=['s','m','l'],help="size of which offset is a proportion of. Is the width of annotations by default")

args = parser.parse_args()

class littlegenomes():
	def __init__(self, geno_f, anno_f, out_f, squeeze, multi, maxheight, incre, shift):
		# Set up of object variables
		self.names =[] # List of names for genomes to produce graphics for
		self.data = {} # Dictionary of annotation data, with genome names as keys
		self.grouping = {} # Dictionary of which genomes belong to which groups
		self.x = 50 # x and y are the coordinates for the top left corner of the diagram
		self.y = -1000 
		self.yincre = incre # Distance objects i.e. the genome diagrams
		self.maxHeight = maxheight # Entire figure cannot exceed this number of px for height and width
		self.maxWidth = 550.0 
		self.xscale = [1] # This guy is a list, because if we group genomes each group gets its own xscaling
		self.yscale = 1 
		# x and y scale are applied to diagrams to fit them within the specifed max width and height
		self.shift = shift
		self.nameGap = 10 # Distance between name label and the end of the genome diagram
		self.tick = 500.0 # Number of nucleotides at which the scalebar shows a tick
		
		# These styles are specifically for getting text to display in the correct position when 
		# opened in inkscape
		self.textStyle = "font-variant:normal;font-weight:bold;font-stretch:normal;font-size:11.19999981px;font-family:Arial;-inkscape-font-specification:Arial-BoldMT;writing-mode:lr-tb;fill:#231f20;fill-opacity:1;fill-rule:nonzero;stroke:none"
		self.annoTextStyle = "font-variant:normal;font-weight:bold;font-stretch:normal;font-size:10px;font-family:Arial;-inkscape-font-specification:Arial-BoldMT;writing-mode:lr-tb;fill:#231f20;fill-opacity:1;fill-rule:nonzero;stroke:none;"+"text-anchor:middle;dominant-baseline:hanging;"
		self.sbTextStyle = "font-variant:normal;font-weight:bold;font-stretch:normal;font-size:10px;font-family:Arial;-inkscape-font-specification:Arial-BoldMT;writing-mode:lr-tb;fill:#231f20;fill-opacity:1;fill-rule:nonzero;stroke:none;text-anchor:middle;"
		
		# The genome diagrams represent rectangles to represent annotations. The following are some constants
		# for use when creating these annotation representaitons
		self.annoWidth = {'s':15, 'm':22, 'l':30}
		self.annoOffset = lambda x:0.5*x
		self.annoColour ={'orange':svgwrite.rgb(246,146,30),
				'blue':svgwrite.rgb(35,192,240),
				'red':svgwrite.rgb(236,28,36),
				'green':svgwrite.rgb(62,144,78),
				'light green':svgwrite.rgb(139,197,63),
				'brown':svgwrite.rgb(195,153,107),
				'pink':svgwrite.rgb(224,117,173),
				'gray':svgwrite.rgb(146,148,151),
				'yellow':svgwrite.rgb(249,236,49),
				'purple':svgwrite.rgb(228,116,255),
				'black':'black',
				'white':'white'}
		
		# Open given files, extract the data and format for use by object
		self.loadData(geno_f,anno_f)

		# After we know how many genomes will be displayed and their lengths, we can set scaling constants
		lengths = [self.data[i][0] for i in self.names]
		if squeeze:
			scales = self.scaleData(squeeze,max(lengths),len(self.names)) 
			self.xscale[0] = scales[0]
			self.yscale = scales[1]
		if multi:
                        self.groupGenomes(squeeze,lengths)
		
		# Sort annotations in order of layer (makes it easy to render them in the correct order)
		for i in self.names:
			size = self.data[i].pop(0)
			self.data[i].sort(key=lambda x:x[3])
			self.data[i].insert(0,size)
		
		# Set up the svgwrite object - we give it commands, it creates the svg
		self.dwg = svgwrite.Drawing(out_f, debug=True)

	def main(self):
		loc = (self.x,self.y)
		
		# For each group of genomes (without --multi this only runs once)
		group = 0
		while group < len(self.grouping):
			group_names = self.grouping[group]
			xscale = self.xscale[group]
			
			sb_length = ceil(max([self.data[i][0] for i in group_names])/self.tick)*self.tick
			self.drawScaleBar(loc, sb_length, self.xscale[group])			
			loc = (loc[0],loc[1]+(self.yincre*self.yscale))

			# For each genome to be displayed
			i = 0
			while i < len(group_names):
				self.drawGenome(group_names[i],loc,self.xscale[group])

				loc = (loc[0],loc[1]+(self.yincre*self.yscale))
				i += 1
			
			group += 1
		
		# After scalebar and genome diagrams have been designed, svgwrite saves to file
		self.dwg.save()		

	def loadData(self, geno_f, anno_f):
		# Get the names and lengths of the genomes
		f = open(geno_f,'r')
		for line in f.readlines():
			line = line.strip().split('\t')
			self.names.append(line[0]) # Names go here
			self.data[line[0]] = [int(line[1])] # Name is dict key to list, list starts with genome length
		f.close()
		self.grouping[0] = self.names # Sets grouping to default, where all genomes belong to one group
		
		# Get the annotation data from the table and store in list (self.data[name])
		f = open(anno_f,'r')
		for line in f.readlines():
			line = line.strip().split('\t')
			if len(line) != 8:
				raise Exception(
"Error: Incorrect number of fields in annotation table input.\n Detected: {} instead of 8.".format(len(line)))
			# Need a test for if line[3] is a float here
			if not line[5] in self.annoColour.keys():
				raise Exception("Error: Colour option incorrect:\n",line[5],line)
			elif not line[7] in self.annoWidth.keys():
				raise Exception("Error: Width option incorrect:\n",line[7],line)

			try:
				anno = (int(line[1]), int(line[2]), line[3], int(line[4]), line[5], line[6], line[7])
				
				if not anno[4] in self.annoColour.keys():
					raise Exception("Colour option incorrect") 
				elif not anno[6] in self.annoWidth.keys():
					raise Exception("Width option incorrect")
			except:
				print("Error: annotation is in incorrect format:\n",line,
					"\nCorrect format should fit:\n",
					"<genome name>\\t<start position in genome>\\t<end position in genome>\\t<offset from line>{-1 being below, 0 on line, 1 above}\\t<layer {0,1,2,3}>\\t<colour>\\t<annotation label>\\t<width of annotation graphic {s,m,l}>" )
				sys.exit()
			if not line[0] in self.data.keys():
				raise Exception("This annotation's name does not correspond to a genome:\n",
					line)
			self.data[line[0]].append(anno) # Probably should check value is in self.names first
		f.close()
	
	def scaleData(self, squeeze, scale_to, num_genomes):
		# Work out the height of this figure - how many genomes spaced how far apart?
		height = num_genomes*self.yincre
		
		# Won't scale the height of the genomes unless it looks like we will need more than a A4 page worth
		if height > self.maxHeight and squeeze != "x":
			yscale = self.maxHeight/height
		else:
			yscale = self.yscale
		
		# Scale the x axis
		if squeeze != "y":
			xscale = self.maxWidth/scale_to
		else:
			xscale = self.xscale
		
		return((xscale,yscale))
	
	def groupGenomes(self, squeeze, lengths):
		# Get indexes for sorted list - will correspond to self.names
		sorted_lengths = sorted(lengths,reverse=True)
		sorted_indexs = sorted(range(len(lengths)), key=lambda x: lengths[x], reverse=True)		
		
		# Go through and choose where the groups start and end
		# Will make groups where members are at least 80% of the scalebar and assigns x-scales for each group
		longest = max(lengths)
		threshold = 0.8 # 80% of the scalebar
		group = 0
		i = 0
		self.grouping[group] = []
		while i < len(sorted_lengths):
			length = sorted_lengths[i]
			if length < threshold*longest:
				longest = length
				group += 1
				self.grouping[group] = []
				self.xscale.append(self.scaleData(squeeze,longest,1)[0])
			self.grouping[group].append(self.names[sorted_indexs[i]]) # Change self.grouping to assign genomes to groups 
			i += 1

	def drawGenome(self, name, start, xscale):
		dwg = self.dwg

		# Draw base line of genome
		length = self.data[name][0] * xscale
		end = (start[0]+length, start[1])
		dwg.add(dwg.line(start, end, stroke='black')) # Have to work out the rgb bit svgwrite.rgb(10,10,16,'%'))

		# Add genome title
		dwg.add(dwg.text(name, insert=(end[0]+self.nameGap,end[1]), style=self.textStyle))
		
		# For each annotation 
		for a in self.data[name][1:]:
			wid = self.annoWidth[a[6]]*self.yscale
			shift = self.annoWidth[self.shift]*self.yscale if self.shift else wid
			offset = float(a[2])
			pos = (start[0]+(a[0]*xscale), start[1]-(wid+shift*offset)/2)
			dim = ((a[1]-a[0])*xscale, wid*self.yscale)
			# Add rectangle
			dwg.add(dwg.rect(insert=pos,size=dim,stroke='black',fill=self.annoColour[a[4]]))
			# Add text
			dwg.add(dwg.text(a[5],insert=(pos[0]+dim[0]/2,pos[1]+(wid-6)/2*self.yscale),
				style=self.annoTextStyle)) # Kind of a messy hack to get the text to sit in the right place
	
	def drawScaleBar(self, start, length, xscale):
		# I want the length of the scalebar (in pixels) to be the largest size rounded up to the nearest 500bp
		scale_len = length*xscale 
		scale_incre = self.tick*xscale

		dwg = self.dwg

		# Creates a group - the scalebar is made up of a bunch of line objects
		sb = dwg.add(dwg.g())
		# Creates the 'spine' of the scalebar
		end = (start[0]+scale_len, start[1])
		sb.add(dwg.line(start, end, stroke='black'))
		
		# Ticks are lines on the scale bar that represent a distance in nucleotides. Every second tick 
		# is emphasized and given a label to show the distance in nucleotides
		tick = 0 # Current tick
		emph_size = (8,5) # Dimensions for a larger than normal line
		emph = 0 # Are we emphasizing the current tick
		nt = 0 # Distance in nucleotides along the scalebar
		while tick <= scale_len+(scale_len*0.001): # There is an extra thousandth there to handle the float error
			s = (start[0]+tick, start[1])
			e = (start[0]+tick, start[1]-emph_size[emph])
			if not nt == 0: # We never emphasize tick 0
				sb.add(dwg.line(s, e,stroke='black'))
				if tick >= scale_len:
					dwg.add(dwg.text("{} nt".format(int(nt)), e, style=self.sbTextStyle))
				elif emph == 0:
					dwg.add(dwg.text(int(nt), e, style=self.sbTextStyle))
			nt += self.tick
			tick += scale_incre
			emph = 1 if emph == 0 else 0
		
lg = littlegenomes(args.geno_f, args.anno_f, args.out_f, args.squeeze, args.multi, args.max_height, args.incre_height,args.offset_shift)
lg.main()
