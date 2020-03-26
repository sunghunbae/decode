from __future__ import print_function

from itertools import combinations
from optparse import OptionParser

import sys
import csv
import os
import re

class Transcript(object):
    def __init__(self, filename):
        self.terminal = sys.stdout
        self.logfile = open(filename, "w")
    def write(self, message):
        self.terminal.write(message)
        self.logfile.write(message)
    def flush(self):
        # python 3 compatibility
        pass


class DNAEncoded(object):
    def __init__(self, options, silent=False):
        """
        Internal variables

            comp: complementary sequence definitions
            options: command options
            silent: verbose or not
            scheme: DNA encoding scheme
            regex: regular expression for re.search
            regex_build: progressive build up of regular expressions
                for testing partial match
            numbb: number of building blocks
            numrs: number of random sequence (should be 0 or 1)
            bbtag: map of sequence to building block name
            U: unique random sequences if numrs > 0
            C: hit counts
            S: sum of hit counts
            Cp: partly matched hit counts
            Sp: partly matched hit counts sum
        
        Regular expression search

            additional sequences at 5' or 3' end does not matter
            regular expression search is strict, 
            and any mismatch inside the regular expression results in 'no match'

        Data structure of 'C', combined building block count
        
            C[1] : single building block count
                C[0][1]['1']['bb1'] --> count

            C[2] : combined building blocks(2) count
                C[0][2]['12']['bb1 bb2'] --> count
                C[0][2]['23']['bb2 bb3'] --> count

            C[3] : combined building blocks(3) count
                C[0][3]['123']['bb1 bb2 bb3'] --> count
        """

        self.comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        self.options = options
        self.silent = silent
        self.scheme = ''
        self.regex = ''
        self.regex_build = [] 
        self.numbb = 0
        self.numrs = 0
        self.bbtag = {}
        self.U = {}
        self.C = {}
        self.S = {}
 
        self.set_regular_expression(options.encoding)
        self.map_seq_to_building_block(options.bbseq)


    def set_regular_expression(self, filename):
        with open (filename, 'r') as f:
            self.scheme = f.read()
            self.scheme = re.sub('[\s\t\n\r]', '', self.scheme)
        nextpos = 0
        for m in re.finditer('[\{\(]\d+[\}\)]', self.scheme):
            n = re.search('\d+', m.group(0)).group(0)  # string
            if m.group(0).startswith('{'):
                self.numbb += 1
                var = "(?P<bb%d>.{%s})" % (self.numbb, n)
            else:
                self.numrs += 1
                var = "(?P<rnd>.{%s})" % n
            self.regex_build.append(self.scheme[nextpos:m.start()] + var)
            nextpos = m.end()
        self.regex_build.append(self.scheme[nextpos:])
        self.regex = ''.join(self.regex_build)
        if not self.silent:
            print("Encoding scheme")
            print("  filename= {}".format(filename))
            print("  scheme= {}".format(self.scheme))
            # print("  regex_build= {}".format(self.regex_build))
            # print("  regex= {}".format(self.regex))
            print("  number of building blocks= {}".format(self.numbb))
            print("  number of random sequence blocks= {}".format(self.numrs))
            print()


    def map_seq_to_building_block(self, filename):
        with open(filename, 'r') as f:
            l = 0
            for line in f:
                l += 1
                line = line.strip()
                if line.startswith("#"): continue
                col = line.split()
                try:
                    (seq, mol, cyc) = col[:3]
                except:
                    print("[{}] {}".format(l, line))
                    print("[{}] valid format: <sequence> <building-block-id> <cycle-number>".format(l))
                    sys.exit(1)
                if not cyc in self.bbtag:
                    self.bbtag[cyc] = {}
                if seq in self.bbtag[cyc]:
                    print("[{}] {}".format(l, line))
                    print("[{}] sequence assigned to a building block {} for cycle {}".format(
                        l, self.bbtag[cyc][seq], cyc))
                    sys.exit(2)
                else:
                    self.bbtag[cyc][seq] = mol
            if not self.silent:
                print("Building block sequences")
                print("  filename= {}".format(filename))
                for cyc in self.bbtag:
                    print("  cycle {} number of building block sequences= {}".format(
                        cyc,len(self.bbtag[cyc])))
                print()


    """ return reverse complementary sequence """
    def revcomp(self,s):
        l = len(s)
        revcomp = ''
        for i in range(l-1, -1, -1):
            revcomp += self.comp[s[i]]
        return revcomp


    """ look up a matching building block sequences """
    def lookup(self, args):
        # initialize count
        self.C = {}
        self.U = {}
        self.S = {}
        for f in range(len(args)): # number of .fastq files
            self.C[f]= {}
            self.U[f]= {}
            self.S[f]= {}
            for i in range(1, self.numbb + 1):
                self.C[f][i]= {}
                self.S[f][i]= {}
                for j in combinations(range(1, self.numbb + 1), i):
                    classId= ''.join(map(str, j))
                    self.C[f][i][classId]= {}
                    self.S[f][i][classId]= 0
        if not self.silent:
            print("Reading FASTQ file(s)......")
        
        total_success = 0
        for fileIdx, filename in enumerate(args):
            linenum = 0
            reading = 0
            success = 0
            partly  = 0
            with open(filename, 'rt') as fastq:
                if not self.silent:
                    print("  S{}= {}".format(fileIdx, filename))
                for line in fastq:
                    linenum += 1
                    if linenum % 4 == 2:  # sequence info resides at the second line
                        s = line.strip()
                        reading += 1
                        # note: regular expression is strict
                        blockId = self.lookup_subprocess(s,self.regex)
                        if not blockId and self.options.revcomp: # try again
                            blockId = self.lookup_subprocess(self.revcomp(s),self.regex)
                        if blockId:
                            success += 1
                            self.count(fileIdx, blockId)
                        elif options.partly: # why doesn't it strictly match?
                            partial_blockId = []
                            for i in range(len(self.regex_build)-1,0,-1):
                                partial_regex = ''.join(self.regex_build[:i])
                                blockId = self.lookup_subprocess(s,partial_regex)
                                if not blockId and self.options.revcomp: # try again
                                    blockId = self.lookup_subprocess(self.revcomp(s),partial_regex)
                                if blockId:
                                    partial_blockId = blockId
                                    break
                            if partial_blockId:
                                partly += 1
                                if not self.silent:
                                    print("  [{:6d} / {:6d}] partly matched {}".format(
                                        partly, reading, ','.join(partial_blockId)))
                total_success += success
                if not self.silent:
                    if options.partly:
                        print("  number of sequences= {}, success= {}, partly= {}".format(
                            reading, success, partly))
                    else:
                        print("  number of sequences= {}, success= {}".format(reading, success))

        if not self.silent:
            print()
        if total_success > 0:
            self.writecsv(args)
            print()
        elif not self.silent:
            print("Skipping .csv files since there is no matched sequence")
            print()



    def lookup_subprocess(self,s,regex):
        blockId = []
        m = re.search(regex, s)
        if not m: return blockId

        for i in range(1, self.numbb + 1):
            try:
                seq = m.group('bb%d' % i)
                bb = self.bbtag[str(i)][seq]
                blockId.append(bb)
            except IndexError:
                continue
            except KeyError:
                continue
        if len(blockId) == self.numbb:
            if self.numrs > 0:
                try:
                    rndseq = m.group('rnd')
                    blockId.append(rndseq)
                except IndexError:
                    pass
            return blockId
        return []


    def count (self,fileIdx,blockId):
        for i in range(1, self.numbb + 1):
            for j in combinations(range(1, self.numbb + 1), i):
                classId = ''.join(map(str, j))
                combiId = []
                for k in j:
                    combiId.append(blockId[k - 1])
                combiId = ' '.join(combiId)
                if combiId in self.C[fileIdx][i][classId]:
                    self.C[fileIdx][i][classId][combiId] +=1
                else:
                    self.C[fileIdx][i][classId][combiId] = 1
                self.S[fileIdx][i][classId] += 1
        # random sequence
        if self.numrs > 0:
            combiId = ' '.join(blockId[:self.numbb])
            randomSeq = blockId[-1]
            if combiId in self.U[fileIdx]:
                self.U[fileIdx][combiId]['total'] += 1
                self.U[fileIdx][combiId]['set'].add(randomSeq)
            else:
                self.U[fileIdx][combiId] = {}
                self.U[fileIdx][combiId]['total'] = 1
                self.U[fileIdx][combiId]['set'] = {randomSeq}


    def writecsv(self, args):
        if not self.silent:
            print("Writing hit counts for building block combinations")
        for i in range(1, self.numbb + 1):
            for j in combinations(range(1, self.numbb + 1), i):
                classId = ''.join(map(str, j))
                outfilename = options.prefix + '_' + classId + '.csv'
                with open(outfilename, 'w') as csvfile:
                    csvout = csv.writer(csvfile, delimiter=',', quotechar='"', )
                    header = ['Class'] + ['block%d' % k for k in j]
                    for f in self.C :
                        header.append('cnt_S%d' % f)
                        header.append('pop_S%d' % f)
                    csvout.writerow(header)
                    # For a given classId,
                    # populate building block combinations appeared at least once
                    remaining_k= []
                    for f in self.C :
                        remaining_k.extend(list(self.C[f][i][classId].keys()))
                    # rows of data - sorted by the count of the 1st file data
                    for k in sorted(self.C[0][i][classId], key=self.C[0][i][classId].get, reverse=True):
                        row = [classId] + k.split()
                        for f in self.C:
                            if k in self.C[f][i][classId] :
                                try:
                                    p = 100.*self.C[f][i][classId][k]/float(self.S[f][i][classId])
                                except:
                                    p = 0.0
                                row += [self.C[f][i][classId][k], p]
                            else:
                                row += [0, 0.0]
                        remaining_k.remove(k)
                        csvout.writerow(row)

                    # rows of data - unsorted, not in the 1st file data
                    for k in remaining_k :
                        row = [classId] + k.split()
                        for f in self.C:
                            if k in self.C[f][i][classId]:
                                try:
                                    p = 100.*self.C[f][i][classId][k]/float(self.S[f][i][classId])
                                except:
                                    p = 0.0
                                row += [self.C[f][i][classId][k], p]
                            else:
                                row += [0, 0.0]
                        csvout.writerow(row)

        """ write out amplification factor or uniqueness """
        if not self.silent:
            print("Writing sequence uniqueness")
        for f in self.U:
            outfilename = '{}_S{}_uniqueness.csv'.format(options.prefix,f)
            with open(outfilename, 'w') as csvfile:
                csvout = csv.writer(csvfile, delimiter=',', quotechar='"', )
                header = ['block%d' % i for i in range(1, self.numbb + 1) ]
                header += ['S%d_unique' % f, 'S%d_total' % f, 'S%d_unique/total' % f]
                csvout.writerow(header)
                for k in self.U[f]:
                    row = k.split()
                    if k in self.U[f] :
                        total = self.U[f][k]['total']
                        unique = len(self.U[f][k]['set'])
                        row += [unique, total, float(unique) / float(total)]
                    else:
                        row += [0, 0, 0]
                    csvout.writerow(row)


if __name__ == "__main__" :

    parser= OptionParser()

    parser.add_option("-p", "--prefix", dest="prefix", default="decoded",
        help="output filename prefix", metavar="PREFIX")
    parser.add_option("-e", "--encoding", dest="encoding", default="encoding.txt",
        help="input encoding scheme file", metavar="ENCODING")
    parser.add_option("-b", "--bbseq", dest="bbseq", default="BBS.txt",
        help="input building block sequence file", metavar="BBS")
    parser.add_option("-r", "--revcomp", dest="revcomp", action="store_true",
        default=False, help="option to try reverse complementary seq.")
    parser.add_option("-y", "--partly", dest="partly", action="store_true",
        default=False, help="option to report partly matched sequences")

    (options, args) = parser.parse_args()
    if len(args) == 0 :
        parser.print_help()
        sys.exit(0)

    # logging
    sys.stdout = Transcript("{}.log".format(options.prefix))
    echo_commands = "$ "+' '.join(sys.argv)
    print(echo_commands)

    Lib = DNAEncoded(options)
    Lib.lookup(args)