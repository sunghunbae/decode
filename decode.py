import sys
import csv
import os
import re

from itertools import combinations
from optparse import OptionParser

""" return a regular expression representing the encoding scheme """

def read_encoding ( filename ) :
	f= open(filename, 'r')
	scheme= f.read()
	scheme= re.sub('[\s\t\n\r]','',scheme)
	f.close()

	regex= ""
	nextpos= idx= rnd= 0
	for m in re.finditer('[\{\(]\d+[\}\)]', scheme):
		n= re.search('\d+', m.group(0)).group(0) # string
		if m.group(0).startswith('{') :
			idx+= 1
			var= "(?P<bb%d>.{%s})" % (idx, n)
		else:
			rnd= 1
			var= "(?P<rnd>.{%s})" % n
		regex += scheme[nextpos:m.start()] + var
		nextpos= m.end()
	regex += scheme[nextpos:]
	return regex, idx, rnd


""" read a building block sequence file and return its map object """

def read_building_block_sequence ( filename ) :
	mapobj= {}
  	f= open(filename, 'r')
	for line in f :
		col= line.strip().split()
		(seq,mol,cyc) = (col[0],col[1],col[2])
		if not cyc in mapobj :
			mapobj[cyc]= {}
		if seq in mapobj[cyc] :
			print "Error (check duplicate)", line.strip()
		else:
			mapobj[cyc][seq]= mol
	return mapobj


""" return reverse complementary sequence """

Comp= {'A':'T','T':'A','G':'C','C':'G','N':'N'}
def revcomp (s) :
	l= len(s)
	revcomp= ''
	for i in range(l-1,-1,-1) :
		revcomp += Comp[s[i]]
	return revcomp
	
"""
handling command line options
	prefix <- options.prefix
	args <- filename(s)
"""

parser = OptionParser()

parser.add_option(
	"-p", 
	"--prefix", 
	dest= "prefix", 
	default= "decode",
  	help= "output filename prefix", 
	metavar= "PREFIX")

parser.add_option(
	"-e", 
	"--encoding",
	dest= "encoding", 
  	default= "encoding.txt",
	help= "input encoding scheme file",
	metavar= "ENCODING")

parser.add_option(
	"-b", 
	"--bbseq",
	dest= "bbseq", 
  	default= "BBS.txt",
	help= "input building block sequence file",
	metavar= "BBS")

parser.add_option(
	"-r",
	"--revcomp",
	dest= "revcomp",
	action= "store_true",
	default= False, 
	help="option for reverse complementary trial")

(options, args) = parser.parse_args()

if options.revcomp : 
	options.revcomp= True


# read encoding scheme
(scheme, numbb, numrs)= read_encoding ( options.encoding )
print "read encoding scheme"
print "regex=", scheme
print "number of building blocks=", numbb
print "number of random sequence blocks=", numrs


# read building block sequence info
BBTag= read_building_block_sequence ( options.bbseq )
print "read building block sequence info."

""" setup Count for combined building blocks

  	Count[1] : single building block count
		Count[1]['1']['bb1'] --> count

	Count[2] : combined building blocks(2) count
		Count[2]['12']['bb1 bb2'] --> count
		Count[2]['23']['bb2 bb3'] --> count

	Count[3] : combined building blocks(3) count
		Count[3]['123']['bb1 bb2 bb3'] --> count
"""

Unique= {}
Count= {}	
for i in range(1, numbb+1) :
	Count[i]= {}
	for j in combinations(range(1,numbb+1),i) :
		Count[i][''.join(map(str,j))]= {}

logfile= open(options.prefix+'.log','wb')

def lookup(seq) :
	blockId= []
	try:
		m= re.search(scheme, seq)
		for i in range(1,numbb+1) :
			var= m.group('bb%d' % i)
			blockId.append( BBTag['%d' % i][var] )
		if numrs > 0 :
			blockId.append( m.group('rnd') )
	except:
		return []
	return blockId

# filenames: .fastq files
for filename in args :

	print "reading", filename

	linenum= 0
	reading= 0
	success= 0

	fastaq= open(filename,'r')
	for line in fastaq :
		linenum += 1 # linenum is 1 at the first line
		if linenum % 4 == 2 : # sequence info resides at the second line	
			reading += 1
			raw= line.strip()
			blockId= lookup(raw)
			if not blockId :
				if options.revcomp: # try again
					blockId= lookup(revcomp(raw))
					if not blockId:
						continue
				else:
					continue
					
			success += 1

			# store blockId
			for i in range(1,numbb+1) :
				for j in combinations(range(1,numbb+1), i) :
					classId= ''.join(map(str,j))
					combiId= []
					for k in j :
						combiId.append(blockId[k-1])
					combiId= ' '.join(combiId)
					if combiId in Count[i][classId] :
						Count[i][classId][combiId] += 1
					else:
						Count[i][classId][combiId]  = 1

			# store random sequence
			combiId= ' '.join(blockId[:numbb])
			randomSeq= blockId[-1]
			if combiId in Unique:
				Unique[combiId]['total'] += 1
				Unique[combiId]['set'].add(randomSeq)
			else:
				Unique[combiId]= {}
				Unique[combiId]['total'] = 1
				Unique[combiId]['set']= set([randomSeq])

	print "read", reading, "sequences"
	print "success", success

	logfile.write(filename+"\n")
	logfile.write("reading %d\n" % reading)
	logfile.write("success %d\n" % success)

# write out csv output

for i in range(1,numbb+1) :
	for j in combinations(range(1,numbb+1), i) :
		classId= ''.join(map(str,j))
		csvfile= open(options.prefix+'_'+classId+'.csv','wb')
 		csvout= csv.writer(csvfile, delimiter=',', quotechar='"',)
		header= ['Class']
		for k in j :
			header.append('block%d' % k)
		header.append('count')
		csvout.writerow( header )
		for k in sorted(Count[i][classId], 
			key= Count[i][classId].get, reverse= True) :
			row= [ classId ]
			row.extend ( k.split() )
			row.append ( Count[i][classId][k] )
			csvout.writerow( row )

# write out amplification factor

csvfile= open(options.prefix+'_unique.csv','wb')
csvout= csv.writer(csvfile, delimiter=',', quotechar='"',)
header= []
for i in range(1,numbb+1) :
 	header.append ('block%d' % i )
header.extend (['unique','total','unique/total'])
csvout.writerow( header )
for k in Unique:
	total= Unique[k]['total']
	unique= len(Unique[k]['set'])
	row= k.split()
	row.extend( [unique, total, float(unique)/float(total)] )
	csvout.writerow ( row )
