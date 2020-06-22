#!/usr/bin/python3

# Callum Le Lay 190124
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
parser.add_argument("--squeeze", choices=['x','y',"both"], default=None, help="genomes are scaled by height to fit in a A4 sized page")
parser.add_argument("--max_height",default=950,type=float,help="max size of figure, if y-axis squeeze is used will scale all genomes to fit within this height - defaults to A4 height")
parser.add_argument("--incre_height",default=100,type=int, help="space between individual sequences on the y-axis")

args = parser.parse_args()

class littlegenomes():
	def __init__(self, geno_f, anno_f, out_f, squeeze,maxheight,incre):
		# Set up of object variables
		self.squeeze = squeeze
		self.names =[] # List of names for genomes to produce graphics for
		self.data = {} # Dictionary of annotation data, with genome names as keys
		self.x = 50 # x and y are the coordinates for the top left corner of the diagram
		self.y = -1000 
		self.yincre = incre # Distance objects i.e. the genome diagrams
		self.maxHeight = maxheight # Entire figure cannot exceed this number of px for height and width
		self.maxWidth = 550.0 
		self.xscale = 1 
		self.yscale = 1 
		# x and y scale are applied to diagrams to fit them within the specifed max width and height
		self.nameGap = 10 # Distance between name label and the end of the genome diagram
		self.textLength = 100 # Obsolete DELETE
		self.tick = 500.0 # Number of nucleotides at which the scalebar shows a tick
		
		# These styles are specifically for getting text to display in the correct position when 
		# opened in inkscape
		#self.lineStyle = DELETE
		self.textStyle = "font-variant:normal;font-weight:bold;font-stretch:normal;font-size:11.19999981px;font-family:Arial;-inkscape-font-specification:Arial-BoldMT;writing-mode:lr-tb;fill:#231f20;fill-opacity:1;fill-rule:nonzero;stroke:none"
		self.annoTextStyle = "font-variant:normal;font-weight:bold;font-stretch:normal;font-size:10px;font-family:Arial;-inkscape-font-specification:Arial-BoldMT;writing-mode:lr-tb;fill:#231f20;fill-opacity:1;fill-rule:nonzero;stroke:none;"+"text-anchor:middle;dominant-baseline:hanging;"
		self.sbTextStyle = "font-variant:normal;font-weight:bold;font-stretch:normal;font-size:10px;font-family:Arial;-inkscape-font-specification:Arial-BoldMT;writing-mode:lr-tb;fill:#231f20;fill-opacity:1;fill-rule:nonzero;stroke:none;text-anchor:middle;"
		
		# The genome diagrams represent rectangles to represent annotations. The following are some constants 
		# for use when creating these annotation representaitons
		self.annoWidth = {'s':15, 'm':22, 'l':30}
		self.annoOffset = {'-1':0, '0':0.5, '1':1}
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
		if self.squeeze:
			self.scaleData(self.squeeze) 

		# Sort annotations in order of layer (makes it easy to render them in the correct order)
		for i in self.names:
			size = self.data[i].pop(0)
			self.data[i].sort(key=lambda x:x[3])
			self.data[i].insert(0,size)
		
		# Set up the svgwrite object - we give it commands, it creates the svg
		self.dwg = svgwrite.Drawing(out_f, debug=True)

	def main(self):
		# Creates scalebar from top left corner
		i = 0
		loc = (self.x,self.y)
		self.drawScaleBar(loc)

		# For each genome to be displayed
		while i < len(self.names):
			loc = (loc[0],loc[1]+(self.yincre*self.yscale))
			name = self.names[i]
			self.drawGenome(name,loc)
			i += 1
		
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
		
		# Get the annotation data from the table and store in list (self.data[name])
		f = open(anno_f,'r')
		for line in f.readlines():
			line = line.strip().split('\t')
			if not line[3] in self.annoOffset.keys():
				raise Exception("Error: Readframe option incorrect:\n",line[3],line)
			elif not line[5] in self.annoColour.keys():
				raise Exception("Error: Colour option incorrect:\n",line[5],line)
			elif not line[7] in self.annoWidth.keys():
				raise Exception("Error: Width option incorrect:\n",line[7],line)

			try:
				anno = (int(line[1]), int(line[2]), line[3], int(line[4]), line[5], line[6], line[7])
				
				if not anno[2] in self.annoOffset.keys():
					raise Exception("Readframe option incorrect")
				elif not anno[4] in self.annoColour.keys():
                                        raise Exception("Colour option incorrect") 
				elif not anno[6] in self.annoWidth.keys():
                                        raise Exception("Width option incorrect")
			except:
				print("Error: annotation is in incorrect format:\n",line,
					"\nCorrect format should fit:\n",
					"<genome name>\\t<start position in genome>\\t<end position in genome>\\t<readframe{-3,-2,-1,0,1,2,3}>\\t<layer {0,1,2,3}>\\t<colour>\\t<annotation label>\\t<width of annotation graphic {s,m,l}>" )
				sys.exit()
			if not line[0] in self.data.keys():
				raise Exception("This annotation's name does not correspond to a genome:\n",
					line)
			self.data[line[0]].append(anno) # Probably should check value is in self.names first
		f.close()
	
	def scaleData(self, which):
		# Work out the height of this figure - how many genomes spaced how far apart?
		height = len(self.names)*self.yincre
		
		# Won't scale the height of the genomes unless it looks like we will need more than a A4 page worth
		if height > self.maxHeight and which != "x":
			self.yscale = self.maxHeight/height
		
		# Scale the x axis
		if which != "y":
			width = max([self.data[i][0] for i in self.names])
			self.xscale = self.maxWidth/width

	def drawGenome(self, name, start):
		dwg = self.dwg

		# Draw base line of genome
		length = self.data[name][0]*self.xscale
		end = (start[0]+length,start[1])
		dwg.add(dwg.line(start, end, stroke='black')) # Have to work out the rgb bit svgwrite.rgb(10,10,16,'%'))

		# Add genome title
		dwg.add(dwg.text(name, insert=(end[0]+self.nameGap,end[1]), style=self.textStyle))

		# For each annotation 
		for a in self.data[name][1:]:
			wid = self.annoWidth[a[6]]
			pos = (start[0]+(a[0]*self.xscale),start[1]-self.annoOffset[a[2]]*self.yscale*wid)
			dim = ((a[1]-a[0])*self.xscale,wid*self.yscale)
			#anno = dwg.add(dwg.g()) <- got sick of ungrouping when reformatting
			# Add rectangle
			dwg.add(dwg.rect(insert=pos,size=dim,stroke='black',fill=self.annoColour[a[4]]))
			# Add text
			dwg.add(dwg.text(a[5],insert=(pos[0]+dim[0]/2,pos[1]+(wid-6)/2*self.yscale),
				style=self.annoTextStyle)) # Kind of a messy hack to get the text to sit in the right place
	
	def drawScaleBar(self, start):
		# I want the length of the scalebar (in pixels) to be the largest size rounded up to the nearest 500bp
		scale_len = ceil((self.maxWidth/self.xscale)/self.tick)*self.tick*self.xscale # ceil rounds decimals
		scale_incre = self.tick*self.xscale
		
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
		while(tick < scale_len):
			s = (start[0]+tick, start[1])
			e = (start[0]+tick, start[1]-emph_size[emph])
			if not nt == 0: # We never emphasize tick 0
				sb.add(dwg.line(s, e,stroke='black'))
			if not nt == 0 and emph == 0:	
				dwg.add(dwg.text(int(nt), e, style=self.sbTextStyle))
			nt += self.tick
			tick += scale_incre
			emph = 1 if emph == 0 else 0

		# Lastly, the units for the scalebar is labelled at the end
		dwg.add(dwg.text('nt', (end[0],end[1]-emph_size[0]), style=self.sbTextStyle))
		
lg = littlegenomes(args.geno_f, args.anno_f, args.out_f, args.squeeze, args.max_height, args.incre_height)
lg.main()

## Error handling
## Better text scaling -> wrap?
## Labels offset when they do not fit
## Add ability to process fasta files
## Individual scalebars?
