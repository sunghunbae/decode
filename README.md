# decode
Decode and analyze NGS(next generation sequencing) data from DNA-Encoded Library screening 

## usage
<pre>
Usage: decode.py [options]

Options:
  -h, --help            show this help message and exit
  -p PREFIX, --prefix=PREFIX
                        output filename prefix
  -e ENCODING, --encoding=ENCODING
                        input encoding scheme file
  -b BBS, --bbseq=BBS   input building block sequence file
  -r, --revcomp         option for reverse complementary trial
</pre>

## encoding file format
A text file with single or multiples lines explain sequence design for encoding building block information and random sequence.

* {#} - sequence of # bases containing building block information. building block 1, 2, ... are assigned from the 5' end of the sequence.
* (#) - sequence of # bases containing random sequence to check out whether hit count is biased by PCR amplification or not.

<pre>
AAATCGATGTG
{6}GAG
{6}AGT
{6}CGA
ACTGAATCTACT
(12)
TCAGACAAGCTTCACCTGC
</pre>

In this example encoding, there are three building blocks with six bases and a random sequence with twelve bases.  Sequences surounding the {} or () blocks should match with the NGS data, but not necessarily from the 5'-end. For example, the following sequence satisfies the encoding scheme.

<pre>
AGTTGACTCCC AAATCGATGTG TGTATG GAG GCTATG AGT GCTGGCCGAACTGAATCTACTAGGGAGAGTGCGTCAGACAAGCTTCACCTGCAATAGATCG
</pre>

## building block data file format
<pre>
GCTGCC FAa-001 1 10K DEL
GCTTGC FAa-002 1 10K DEL
GCTGAC FAa-003 1 10K DEL
GCTGAG FAa-004 1 10K DEL
GCTTTC FAa-005 1 10K DEL
......
</pre>

- column 1 : sequence tag for a building block
- column 2 : name or molecule code for a building block
- column 3 : building cycle number (a same building block can be used with different sequence tags for different building cycles)
- column 4-: extra information

## fastq file format

<pre>
@MG00HS13:1108:H7TNJBCXY:1:1101:1744:2114 1:N:0:CGATGT
AGTTGACTCCCAAATCGATGTGTGTATGGAGGCTATGAGTGCTGGCCGAACTGAATCTACTAGGGAGAGTGCGTCAGACAAGCTTCACCTGCAATAGATCG
+
DDDDDIIIIIIIHIIIHI?HHHHHHHIHIIIIHIIIFGHIHHIHGHHHIIIHIIIIH1CFHIHHI?CCHCEHHHHHHHIIGEHHFHIIIHEHGHIHHIGHH
@MG00HS13:1108:H7TNJBCXY:1:1101:1908:2146 1:N:0:CGATGT
AGTTGACTCCCAAATCGATGTGTGTGGCGAGGCTGGCAGTGCTATCCGAACTGAATCGACTAGAGGCATCCGGTCAGACAAGCTTCACGTGCAATAGATCG
+
B0@<DH@1<11D<FHFH1C<CCGEEEHHH//C/<F=0FG@FEH1FC?C0D/DC<D<@1<<CE<GH1GHGF?C/<<C1<FF@FEECCCG11<@H@EHIH<<E

......
</pre>

## test run
<pre>
cd test
python ../decode.py -b BBS.txt -e encoding.txt -p bae ../data/testlg.fastq
</pre>

## comparison to Zhang et al. (2017)
A related C++ source code was published in the supplementary material to the Zhang et al. (2017) paper on the DNA-encoded library

### credit
<pre>
 Written by Yixin Zhang
 Modified by Hannes Röst, November 2009
 Modified by Fabian Buller, 14.11.2009
 ESACH VERSION 12.2.2010
 Modified by Michael Stravs, 19. Jan. 2012
 Modified by Moreno Wichert & Adrian Rabenseifner, 23. Apr. 2012
</pre>

### structure.txt

It is a main input file describing DNA encoding structure
<pre>
../../data/testlg.fasta
50
23 	24 	28	0
2BB
0
1	12	22	AAATCGATGTG
x	23	28	bb1.txt
2	29	31	GAG
y	32	37	bb2.txt
3	38	40	AGT
z	41	46	bb3.txt
4	47	49	CGA
5	50	61	ACTGAATCTACT
</pre>

- Line 1: FASTA format data file containing NGS data
- Line 2: minimum allowed length of FASTA sequence
Sequences shorter than this length are ignored.
- Line 3: number of building blocks for cycle 1, 2, 3, and 4
- Line 4: output format. Counts are reported for a combination of two building blocks. Other options are 3BB and 4BB. 
- Line 5: maximum allowed mismatches. Sequences with more mismatches are not counted.
- Line 6-: DNA encoding scheme. 
-- Column 1: x, y, z, and $ are reserved for the building block sequences. 1, 2, 3, 4, and 5 are reserved for constant sequence block.
-- Column 2: start position(number, note: the number starts from 1)
-- Column 3: end position(number)
-- Column 4: filename if 1st column is x, y, z or $ and constant sequence if that is 1-5. 

### bb1.txt - building block information

bb2.txt and bb3.txt are similar files for the 2nd and 3rd building blocks.

<pre>
GCTGCC
GCTTGC
GCTGAC
GCTGAG
GCTTTC
.....
</pre>
