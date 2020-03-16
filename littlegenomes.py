#!/usr/bin/python

# Callum Le Lay 190124
# Usage: littlegenomes.py GENOMES_TSV ANNOTATIONS_TSV OUPUT_NAME

# Import
import svgwrite
import sys
from math import ceil

class littlegenomes():
	def __init__(self, geno_f, anno_f, out_f):
		self.names =[]
		self.data = {}
		self.x = 50
                self.y = -1000
                self.yincre = 100
                self.maxHeight = 950.0
                self.maxWidth = 550.0
		self.xscale = 1
		self.yscale = 1
                self.nameGap = 10
		self.textLength = 100
		self.tick = 500.0

		#self.lineStyle =
		self.textStyle = "font-variant:normal;font-weight:bold;font-stretch:normal;font-size:11.19999981px;font-family:Arial;-inkscape-font-specification:Arial-BoldMT;writing-mode:lr-tb;fill:#231f20;fill-opacity:1;fill-rule:nonzero;stroke:none"
		self.annoTextStyle = "font-variant:normal;font-weight:bold;font-stretch:normal;font-size:10px;font-family:Arial;-inkscape-font-specification:Arial-BoldMT;writing-mode:lr-tb;fill:#231f20;fill-opacity:1;fill-rule:nonzero;stroke:none;"+"text-anchor:middle;dominant-baseline:hanging;"
		self.sbTextStyle = "font-variant:normal;font-weight:bold;font-stretch:normal;font-size:10px;font-family:Arial;-inkscape-font-specification:Arial-BoldMT;writing-mode:lr-tb;fill:#231f20;fill-opacity:1;fill-rule:nonzero;stroke:none;text-anchor:middle;"
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
				'black':'black',
				'white':'white'}
		
		self.loadData(geno_f,anno_f)
                self.scaleData()

		# Sort annotations in order of layer (makes it easy to render them in the correct order)
		for i in self.names:
			size = self.data[i].pop(0)
			self.data[i].sort(key=lambda x:x[3])
			self.data[i].insert(0,size)
		
		self.dwg = svgwrite.Drawing(out_f, debug=True) # Full??	

	def main(self):
		i = 0
		loc = (self.x,self.y)
		self.drawScaleBar(loc)		
		while i < len(self.names):
			loc = (loc[0],loc[1]+(self.yincre*self.yscale))
			name = self.names[i]
			self.drawGenome(name,loc)
			i += 1

		self.dwg.save()
			
	def loadData(self, geno_f, anno_f):
		f = open(geno_f,'r')
		for line in f.readlines():
			line = line.strip().split('\t')
			self.names.append(line[0])
			self.data[line[0]] = [int(line[1])]
		f.close()

		f = open(anno_f,'r')
		for line in f.readlines():
			line = line.strip().split('\t')
			anno = (int(line[1]), int(line[2]), line[3], int(line[4]), line[5], line[6], line[7]) # start,end,readframe,layer,colour,label,width - also want some error handling here
			self.data[line[0]].append(anno) # Probably should check value is in self.names first
		f.close()
	
	def scaleData(self):
		# If y exceeds 950 (off the bottom of the A4 template)
		height = len(self.names)*self.yincre
		if height > self.maxHeight:
			self.yscale = self.maxHeight/height
		
		# Always scale the x axis
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
			dwg.add(dwg.text(a[5],insert=(pos[0]+dim[0]/2,pos[1]+dim[1]-(wid-9)/2),
				style=self.annoTextStyle)) # Kind of a messy hack to get the text to sit in the right place
	
	def drawScaleBar(self, start):
		# I want the length of the scalebar (in pixels) to be the largest size rounded up to the nearest 500bp
		scale_len = ceil((self.maxWidth/self.xscale)/self.tick)*self.tick*self.xscale
		scale_incre = self.tick*self.xscale

		dwg = self.dwg
		sb = dwg.add(dwg.g())
		end = (start[0]+scale_len, start[1])
		sb.add(dwg.line(start, end, stroke='black'))
		
		tick = 0
		emph_size = (8,5)
		emph = 0
		nt = 0
		while(tick < scale_len):
			s = (start[0]+tick, start[1])
			e = (start[0]+tick, start[1]-emph_size[emph])
			if not nt == 0:
				sb.add(dwg.line(s, e,stroke='black'))
			if not nt == 0 and emph == 0:	
				dwg.add(dwg.text(int(nt), e, style=self.sbTextStyle))
			nt += self.tick
			tick += scale_incre
			emph = 1 if emph == 0 else 0
		dwg.add(dwg.text('nt', (end[0],end[1]-emph_size[0]), style=self.sbTextStyle))
		
lg = littlegenomes(sys.argv[1],sys.argv[2], sys.argv[3])
lg.main()

## Error handling
## Add ability to process fasta files
## Better text scaling -> wrap?
