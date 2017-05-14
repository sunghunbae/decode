import sys
import csv
import os
import re

from optparse import OptionParser

""" modify the following line for a different DEL library """

design= ''.join([
	'AAATCGATGTG' , 	# opening
	'(?P<bb1>.{6})GAG' , 	# building block 1
	'(?P<bb2>.{6})AGT' , 	# building block 2
	'(?P<bb3>.{6})CGA' , 	# building block 3
	'ACTGAATCTACT' , 	# closing
	'(?P<rnd>.{12})' , 	# random sequence for uniqueness test
	'TCAGACAAGCTTCACCTGC', 	# ending
	])
		
BBTag= {}
BBcombi= {} 	# building block combination
BBcount= {}	# single building block count
BBcount12= {}	# two building block combination count
BBcount13= {}	# two building block combination count
BBcount23= {}	# two building block combination count

Comp= {'A':'T','T':'A','G':'C','C':'G','N':'N'}

""" first column should be sequence
    second column should be building block name
    last column should be cycle number
"""

def read_bbs ( filename ) :
	mapobj= {}
  	f= open(filename, 'r')
	for line in f :
		col= line.strip().split()
		(seq, mol, cyc) = (col[0], col[1], col[-1])
		key= seq + str(cyc) # ex. ACATCA2
		if key in mapobj :
			print "Error", line.strip()
		else:
			mapobj[ key ]= mol
	return mapobj


def reverse_complement (s) :
	l= len(s)
	revcomp= ''
	for i in range(l-1,-1,-1) :
		revcomp += Comp[s[i]]
	return revcomp
	
"""
command line options
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
	"-b", 
	"--bbseq",
	dest= "bbseq", 
  	default= "10k-DEL-BBS.txt",
	help= "input building block sequence file",
	metavar= "BBS")

parser.add_option(
	"-r",
	"--revcomp",
	dest= "revcomp",
	action= "store_true",
	default= False, 
	help="search reverse complementary sequence")

(options, args) = parser.parse_args()

if options.revcomp : 
	options.revcomp= True


# output files
countcsvfile= open(options.prefix+'_count.csv','wb')
countcsvout= csv.writer(countcsvfile, delimiter=',', quotechar='"',)
countcsvout.writerow( ['Building block','count'] )

count12csvfile= open(options.prefix+'_count12.csv','wb')
count12csvout= csv.writer(count12csvfile, delimiter=',', quotechar='"',)
count12csvout.writerow( ['BB1','BB2','count'] )

count13csvfile= open(options.prefix+'_count13.csv','wb')
count13csvout= csv.writer(count13csvfile, delimiter=',', quotechar='"',)
count13csvout.writerow( ['BB1','BB3','count'] )

count23csvfile= open(options.prefix+'_count23.csv','wb')
count23csvout= csv.writer(count23csvfile, delimiter=',', quotechar='"',)
count23csvout.writerow( ['BB2','BB3','count'] )

combicsvfile= open(options.prefix+'_combi.csv','wb')
combicsvout= csv.writer(combicsvfile, delimiter=',', quotechar='"',)
combicsvout.writerow( ['BB1','BB2','BB3','count','unique','ratio'] )

logfile= open(options.prefix+'.log','wb')

# read building block sequence info
BBTag= read_bbs ( options.bbseq )

# filenames: .fastq files
for filename in args :

	linenum= 0
	reading= 0
	success= 0

	fastaq= open(filename,'r')
	for line in fastaq :
		linenum += 1 # linenum is 1 at the first line
		if linenum % 4 == 2 : # sequence info resides at the second line	
			reading += 1
			raw= line.strip()

			"""
			bb1 & bb2 & bb3 & rnd should be defined
			otherwise skip the processing and continue to
			the next sequence
			"""

			try:
				m= re.search(design, raw)
				assert BBTag[m.group('bb1')+'1'] 
				assert BBTag[m.group('bb2')+'2']
				assert BBTag[m.group('bb3')+'3']
				assert m.group('rnd')
			except:
				if options.revcomp: # try again
					try:
						m= re.search(design, \
							reverse_complement(raw))
						assert BBTag[m.group('bb1')+'1']
						assert BBTag[m.group('bb2')+'2']
						assert BBTag[m.group('bb3')+'3']
						assert m.group('rnd')
					except:
						continue
				else:
					continue
					
			success += 1

			# prepare output
			bb= []
			for block,cyc in [
				('bb1','1'),
				('bb2','2'),
				('bb3','3'),
				] :
				k= BBTag[m.group(block)+cyc]
				if k in BBcount:
					BBcount[k] += 1
				else:
					BBcount[k] = 1
				bb.append(k)

			bb12= ' '.join([bb[0],bb[1]])
			if bb12 in BBcount12 :
				BBcount12[bb12] += 1
			else:
				BBcount12[bb12]= 1

			bb13= ' '.join([bb[0],bb[2]])
			if bb13 in BBcount13 :
				BBcount13[bb13] += 1
			else:
				BBcount13[bb13]= 1

			bb23= ' '.join([bb[1],bb[2]])
			if bb23 in BBcount23 :
				BBcount23[bb23] += 1
			else:
				BBcount23[bb23]= 1
	
			rnd= m.group('rnd')
			k= ' '.join ( bb )
			if k in BBcombi :
				BBcombi[k].append ( rnd )
			else :
				BBcombi[k]= [ rnd ]


	print filename
	print "reading", reading
	print "success", success
	logfile.write(filename+"\n")
	logfile.write("reading %d\n" % reading)
	logfile.write("success %d\n" % success)

cmin= None
cmax= None
for k in sorted(BBcount, key= BBcount.get, reverse=True) :
	if not cmin or cmin > BBcount[k] : 
		cmin= BBcount[k]
	if not cmax or cmax < BBcount[k] :
		cmax= BBcount[k]
	countcsvout.writerow ( [k, BBcount[k]] )
print "building block count min %d max %d" % (cmin, cmax)
logfile.write("building block count min %d max %d\n" % (cmin,cmax))

for k in sorted(BBcount12, key= BBcount12.get, reverse=True) :
	count12csvout.writerow ( [k, BBcount12[k]] )
for k in sorted(BBcount13, key= BBcount13.get, reverse=True) :
	count13csvout.writerow ( [k, BBcount13[k]] )
for k in sorted(BBcount23, key= BBcount23.get, reverse=True) :
	count23csvout.writerow ( [k, BBcount23[k]] )

for k in BBcombi :
	count= len(BBcombi[k])
	unique= len(set(BBcombi[k]))
	row= k.split()
	row.append ( count )
	row.append ( unique )
	row.append ( float(unique)/float(count) )
	combicsvout.writerow ( row )
