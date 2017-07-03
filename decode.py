import sys
import csv
import os
import re

from itertools import combinations
from optparse import OptionParser

""" data structure for 'Count', combined building blocks

  	Count[1] : single building block count
		Count[1]['1']['bb1'] --> count

	Count[2] : combined building blocks(2) count
		Count[2]['12']['bb1 bb2'] --> count
		Count[2]['23']['bb2 bb3'] --> count

	Count[3] : combined building blocks(3) count
		Count[3]['123']['bb1 bb2 bb3'] --> count
"""

Unique= {}
Count = {}
CountSum= {}

""" return a regular expression representing the encoding scheme """
def read_encoding(filename):
    f = open(filename, 'r')
    scheme = f.read()
    scheme = re.sub('[\s\t\n\r]', '', scheme)
    f.close()

    regex = ""
    nextpos = idx = rnd = 0
    for m in re.finditer('[\{\(]\d+[\}\)]', scheme):
        n = re.search('\d+', m.group(0)).group(0)  # string
        if m.group(0).startswith('{'):
            idx += 1
            var = "(?P<bb%d>.{%s})" % (idx, n)
        else:
            rnd += 1
            var = "(?P<rnd>.{%s})" % n
        regex += scheme[nextpos:m.start()] + var
        nextpos = m.end()
    regex += scheme[nextpos:]
    return regex, idx, rnd


""" read a building block sequence file and return its map object """
def read_building_block_sequence(filename):
    mapobj = {}
    f = open(filename, 'r')
    linenumber = 0
    for line in f:
        linenumber += 1
        line = line.strip()
        if line.startswith("#"):
            continue

        col = line.split()

        try:
            (seq, mol, cyc) = (col[0], col[1], col[2])
        except:
            print "[Error] expecting 3 columns for each line: <sequence> <building-block-id> <cycle-number>"
            print "[Error] line %d: %s" % (linenumber, line)
            sys.exit(1)

        if not cyc in mapobj:
            mapobj[cyc] = {}

        if seq in mapobj[cyc]:
            print "[Error] identical sequence already assigned to a building block %s for cycle %s" % (
                mapobj[cyc][seq], cyc)
            print "[Error] line %d: %s" % (linenumber, line)
            sys.exit(2)
        else:
            mapobj[cyc][seq] = mol

    return mapobj


""" return reverse complementary sequence """
Comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}

def revcomp(s):
    l = len(s)
    revcomp = ''
    for i in range(l - 1, -1, -1):
        revcomp += Comp[s[i]]
    return revcomp


""" look up a matching building block sequence """
def lookup(seq):
    blockId = []
    try:
        m = re.search(scheme, seq)
        for i in range(1, numbb + 1):
            var = m.group('bb%d' % i)
            blockId.append(BBTag['%d' % i][var])
        if numrs > 0:
            blockId.append(m.group('rnd'))
    except:
        return []
    return blockId


""" parsing command line """
parser= OptionParser()
parser.add_option(
    "-p",
    "--prefix",
    dest="prefix",
    default="decode",
    help="output filename prefix",
    metavar="PREFIX")

parser.add_option(
    "-e",
    "--encoding",
    dest="encoding",
    default="encoding.txt",
    help="input encoding scheme file",
    metavar="ENCODING")

parser.add_option(
    "-b",
    "--bbseq",
    dest="bbseq",
    default="BBS.txt",
    help="input building block sequence file",
    metavar="BBS")

parser.add_option(
    "-r",
    "--revcomp",
    dest="revcomp",
    action="store_true",
    default=False,
    help="option to try reverse complementary seq.")

(options, args) = parser.parse_args()
logfile = open(options.prefix + '.log', 'w')
message = ' '.join(sys.argv)
logfile.write(message + "\n\n")


""" read encoding scheme """
(scheme, numbb, numrs) = read_encoding(options.encoding)
message = "Reading encoding scheme from '%s'\n" % options.encoding
message += "regex= %s\n" % scheme
message += "number of building blocks= %d\n" % numbb
message += "number of random seq. blocks= %d\n" % numrs
print message
logfile.write(message + "\n")


""" read building block sequence info. """
BBTag = read_building_block_sequence(options.bbseq)
message = "Reading building blocks from '%s'\n" % options.bbseq
print message
logfile.write(message + "\n")


""" initialize Count """
for fileIdx in range(len(args)): # .fastq files
    Count[fileIdx]= {}
    Unique[fileIdx]= {}
    CountSum[fileIdx]= {}
    for i in range(1, numbb + 1):
        Count[fileIdx][i]= {}
        CountSum[fileIdx][i]= {}
        for j in combinations(range(1, numbb + 1), i):
            classId= ''.join(map(str, j))
            Count[fileIdx][i][classId]= {}
            CountSum[fileIdx][i][classId]= 0


""" multiple .fastq files for comparison """
for fileIdx, filename in enumerate(args):
    message = "Reading '%s' as S%d\n" % (filename, fileIdx)

    linenum = 0
    reading = 0
    success = 0

    fastaq = open(filename, 'r')
    for line in fastaq:
        linenum += 1  # linenum is 1 at the first line
        if linenum % 4 == 2:  # sequence info resides at the second line
            reading += 1
            raw = line.strip()
            blockId = lookup(raw)
            if not blockId and options.revcomp : # try again
                blockId = lookup(revcomp(raw))
                if not blockId:
                    continue
            else:
                continue

            success += 1

            # store blockId
            for i in range(1, numbb + 1):
                for j in combinations(range(1, numbb + 1), i):
                    classId = ''.join(map(str, j))
                    combiId = []
                    for k in j:
                        combiId.append(blockId[k - 1])
                    combiId = ' '.join(combiId)
                    if combiId in Count[fileIdx][i][classId]:
                        Count[fileIdx][i][classId][combiId] +=1
                    else:
                        Count[fileIdx][i][classId][combiId] = 1
                    CountSum[fileIdx][i][classId] += 1

            # store random sequence
            combiId = ' '.join(blockId[:numbb])
            randomSeq = blockId[-1]
            if combiId in Unique[fileIdx]:
                Unique[fileIdx][combiId]['total'] += 1
                Unique[fileIdx][combiId]['set'].add(randomSeq)
            else:
                Unique[fileIdx][combiId] = {}
                Unique[fileIdx][combiId]['total'] = 1
                Unique[fileIdx][combiId]['set'] = {randomSeq}

    message += "number of seq.= %d (success= %d)\n" % (reading, success)
    print message
    logfile.write(message + "\n")

""" write out csv output """
message = "Writing sequence counts"
print message
logfile.write(message + "\n")

for i in range(1, numbb + 1):
    for j in combinations(range(1, numbb + 1), i):
        classId = ''.join(map(str, j))
        filename = options.prefix + '_' + classId + '.csv'
        csvfile = open(filename, 'w')
        print filename
        logfile.write(filename + "\n")

        # For a given classId,
        # populate building block combinations appeared at least once
        remaining_k= []
        for h in range(len(args)) :
            remaining_k.extend(list(Count[h][i][classId].keys()))

        # header
        csvout = csv.writer(csvfile, delimiter=',', quotechar='"', )
        header = ['Class']
        for k in j:
            header.append('block%d' % k)
        for l in range(len(args)) :
            header.append('count(S%d)' % l)
            header.append('pop.(S%d)' % l)
        csvout.writerow(header)

        # rows of data - sorted by the count of the 1st file data
        for k in sorted(Count[0][i][classId], key=Count[0][i][classId].get, reverse=True):
            row = [classId]
            row.extend(k.split())
            for h in range(len(args)):
                if k in Count[h][i][classId] :
                    row.append(Count[h][i][classId][k])
                    try:
                        row.append(100.*Count[h][i][classId][k]/float(CountSum[h][i][classId]))
                    except:
                        row.append(0.0)

                else:
                    row.append(0)
                    row.append(0.0)
            remaining_k.remove(k)
            csvout.writerow(row)

        # rows of data - unsorted, not in the 1st file data
        for k in remaining_k :
            row = [classId]
            row.extend(k.split())
            for h in range(len(args)):
                if k in Count[h][i][classId]:
                    row.append(Count[h][i][classId][k])
                    try:
                        row.append(100.*Count[h][i][classId][k] / float(CountSum[h][i][classId]))
                    except:
                        row.append(0.0)
                else:
                    row.append(0)
                    row.append(0.0)
            csvout.writerow(row)


""" write out amplification factor """
filename = options.prefix + '_unique.csv'
csvfile = open(filename, 'w')
message = "\nWriting sequence uniqueness to '%s'" % filename
print message
logfile.write(message + "\n")
csvout = csv.writer(csvfile, delimiter=',', quotechar='"', )

# header
header = []
for i in range(1, numbb + 1):
    header.append('block%d' % i)
for h in range(len(args)) :
    header.extend(['unique(S%d)' % h, 'total(S%d)' % h, 'unique/total(S%d)' % h])
csvout.writerow(header)

# rows of data
for k in Unique[h]:
    row = k.split()
    for h in range(len(args)):
        if k in Unique[h] :
            total = Unique[h][k]['total']
            unique = len(Unique[h][k]['set'])
            row.extend([unique, total, float(unique) / float(total)])
        else:
            row.extend([0, 0, 0])
    csvout.writerow(row)