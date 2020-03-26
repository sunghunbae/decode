# decode
Decode and analyze NGS(next generation sequencing) data from DNA-Encoded Library screening 

## How to use (-h or --help)
```
Usage: decode.py [options]

Options:
  -h, --help            show this help message and exit
  -p PREFIX, --prefix=PREFIX
                        output filename prefix
  -e ENCODING, --encoding=ENCODING
                        input encoding scheme file
  -b BBS, --bbseq=BBS   input building block sequence file
  -r, --revcomp         option to try reverse complementary seq.
  -y, --partly          option to report partly matched sequences
```

## How to run the test data

```
cd test
python ../decode.py -e encoding.txt -b BBS.txt -r ../data/testlg.fastq.gz
```

### Example CSV output files

`decoded_1.csv`
```
Class,block1,cnt_S0,pop_S0
1,FAa-004,3519,5.38625197067332
1,FAa-013,3440,5.265332986392788
1,FAa-006,3431,5.251557405905132
1,FAa-001,3387,5.184210123521038
......
```

`decoded_12.csv`
```
Class,block1,block2,cnt_S0,pop_S0
12,FAa-006,Am-010,219,0.335205791866285
12,FAa-009,Am-010,199,0.30459339078260605
12,FAa-013,Am-022,198,0.3030627707284221
12,FAa-009,Am-004,193,0.2954096704575023
......
```

`decoded_123.csv`
```
Class,block1,block2,block3,cnt_S0,pop_S0
123,FAa-006,Am-054,Am-063,16,0.0244899208669432
123,FAa-009,Am-005,Am-007,16,0.0244899208669432
123,FAa-004,Am-006,Am-022,16,0.0244899208669432
...
```

## DNA encoding file format (-e filename or --encoding=filename)

A text file with single or multiple lines describes a DNA sequence design encoding building blocks and a random sequence.

### Syntax

- {#}: Sequence of # bases encoding a building block. The building block numbers 1, 2, ... are assigned sequentially from the 5' end.
- (#): Sequence of # bases encoding a random sequence to test if hit count is biased by PCR amplification.

### Example (encoding.txt)

```
AAATCGATGTG
{6}GAG
{6}AGT
{6}CGA
ACTGAATCTACT
(12)
TCAGACAAGCTTCACCTGC
```

is equivalent to

```
AAATCGATGTG{6}GAG{6}AGT{6}CGAACTGAATCTACT(12)TCAGACAAGCTTCACCTGC
```

Both of above example files describe an identical DNA encoding scheme in which 3 building blocks of 6 bases and a random sequence of 12 bases are placed between defined constant DNA sequence blocks(i.e. opening, cycle, closing, terminal tag sequences). NGS sequences are to be aligned with the encoding scheme. For instance, the following sequence satisfies the above encoding scheme. Note that additional
sequences before 5' end or after 3' end also satisfy the above encoding scheme. But a single nucleotide 
difference inside the scheme results in mismatch.

```AGTTGACTCCCAAATCGATGTGTGTATGGAGGCTATGAGTGCTGGCCGAACTGAATCTACTAGGGAGAGTGCGTCAGACAAGCTTCACCTGCAATAGATCG```

## Building block data file format (-b filename or --bbseq=filename)

### Example (BBS.txt)

```
GCTGCC FAa-001 1 10K DEL
GCTTGC FAa-002 1 10K DEL
GCTGAC FAa-003 1 10K DEL
GCTGAG FAa-004 1 10K DEL
GCTTTC FAa-005 1 10K DEL
......
```

- Column 1: Sequence tag for a building block
- Column 2: Name or molecule code for a building block
- Column 3: Building cycle number. A same building block can be used with different sequence tags for different building cycles.
- Column 4 and the following columns may contain extra information not used in the program.

## FASTQ file format

In FASTQ format, each record has four lines and the second line contains the DNA sequence.

### Example (testlg.fastq)
```
@MG00HS13:1108:H7TNJBCXY:1:1101:1744:2114 1:N:0:CGATGT
AGTTGACTCCCAAATCGATGTGTGTATGGAGGCTATGAGTGCTGGCCGAACTGAATCTACTAGGGAGAGTGCGTCAGACAAGCTTCACCTGCAATAGATCG
+
DDDDDIIIIIIIHIIIHI?HHHHHHHIHIIIIHIIIFGHIHHIHGHHHIIIHIIIIH1CFHIHHI?CCHCEHHHHHHHIIGEHHFHIIIHEHGHIHHIGHH
@MG00HS13:1108:H7TNJBCXY:1:1101:1908:2146 1:N:0:CGATGT
AGTTGACTCCCAAATCGATGTGTGTGGCGAGGCTGGCAGTGCTATCCGAACTGAATCGACTAGAGGCATCCGGTCAGACAAGCTTCACGTGCAATAGATCG
+
B0@<DH@1<11D<FHFH1C<CCGEEEHHH//C/<F=0FG@FEH1FC?C0D/DC<D<@1<<CE<GH1GHGF?C/<<C1<FF@FEECCCG11<@H@EHIH<<E
......
```


## Troubleshooting or rescuing unsuccessful matches (`--partly` or `-y` option)

During building a DNA encoded library, DNA sequence tag is appended as new building blocks are added to a starting scaffold.
These stepwise DNA ligations may introduce a sequence mismatch at certain step and DNA sequence up to that step may be 
still usable by modifying the encoding file (truncating the 3' part from that step). 
Using `-y` or `--partly` option turns on testing and reporting of stepwise matching.
For example, if the encoding has three building blocks and one random sequence block, it would test partial matching of sequences
in the order of bb1-bb2-bb3-rnd, bb1-bb2-bb3, bb1-bb2, bb1 until it finds the longest match.

```
python ../decode.py -e encoding.txt -b BBS.txt -r -p partly -y ../data/testlg.fastq.gz
```

`partly.log`
```
$ ../decode.py -e encoding.txt -b BBS.txt -r -p partly -y ../data/testlg.fastq.gz
Encoding scheme
  filename= encoding.txt
  scheme= AAATCGATGTG{6}GAG{6}AGT{6}CGAACTGAATCTACT(12)TCAGACAAGCTTCACCTGC
  number of building blocks= 3
  number of random sequence blocks= 1

Building block sequences
  filename= BBS.txt
  cycle 1 number of building block sequences= 23
  cycle 2 number of building block sequences= 24
  cycle 3 number of building block sequences= 28

Reading FASTQ file(s)......
  S0= ../data/testlg.fastq.gz
  [     1 /      2] partly matched FAa-018,Am-010,Am-011
  [     2 /      3] partly matched FAa-023,Am-002,Am-001,TCATGCTCACAA
  [     3 /     11] partly matched FAa-018,Am-022,Am-032
  [     4 /     14] partly matched FAa-004,Am-030,Am-062,CCCATCTTTAGT
  [     5 /     23] partly matched FAa-004,Am-062,Am-036
  [     6 /     25] partly matched FAa-018,Am-041,Am-062
  [     7 /     26] partly matched FAa-017,Am-023,Am-023,ACGAATTGCACG
  [     8 /     27] partly matched FAa-002,Am-001,Am-011,CCACTCGGTCGT
  [     9 /     30] partly matched FAa-012,Am-022,Am-041,ATAAATGGGTAT
  [    10 /     31] partly matched FAa-004,Am-063,Am-022,CCCACGTGTCTG
  [    11 /     38] partly matched FAa-010,Am-023,Am-010,ACATTGCAACTG
  [    12 /     43] partly matched FAa-014,Am-005,Am-054,ACGCGAATGGTG
  [    13 /     46] partly matched FAa-023,Am-005,Am-022,GCTTGGGCTTCA
  [    14 /     51] partly matched FAa-012,Am-016,Am-022
  [    15 /     71] partly matched FAa-008,Am-023,Am-062,CGTTCCCTTAGC
  [    16 /     72] partly matched FAa-019,Am-062,Am-022,ATGCAATAAAAT
  [    17 /     74] partly matched FAa-018,Am-017,Am-035,CCATTACTTCTC
  [    18 /     84] partly matched FAa-001,Am-032,Am-002,TAGATGAGTGAA
  [    19 /    114] partly matched FAa-012,Am-005,Am-008,TACCCCTCTCGC
  [    20 /    119] partly matched FAa-007,Am-062,Am-008,GATCTTAAAGTT
......
```

In the above log file, the first partly matched sequence match up to three building blocks and
the second partly matched sequence match up to three building blocks and randome sequence block but fail to match the closing sequence at the very 3' end. We can rescue partly mismatched sequences
by excluding the closing sequence and/or random sequence blocks:


Examples of modified schemes:
AAATCGATGTG{6}GAG{6}AGT{6}CGAACTGAATCTACT(12)
AAATCGATGTG{6}GAG{6}AGT{6}CGA


## Use multiple FASTQ files for comparing different datasets

A sequencing dataset of a target can be compared with a control dataset without the target protein, 
or multiple sequencing datasets of different targets can be compared.
The first FASTQ file is used to sort building blocks and the rest of the FASTQ files are compared
with this first FASTQ file.

```
python ../decode.py -e encoding.txt -b BBS.txt -r ../data/testlg.fastq.gz ../data/testsm.fastq.gz
```


## Comparison to Zhang et al. (2017)

A related C++ source code was published in the supplementary material to the Zhang et al. (2017) paper on the DNA-encoded library.

```
 Written by Yixin Zhang
 Modified by Hannes Röst, November 2009
 Modified by Fabian Buller, 14.11.2009
 ESACH VERSION 12.2.2010
 Modified by Michael Stravs, 19. Jan. 2012
 Modified by Moreno Wichert & Adrian Rabenseifner, 23. Apr. 2012
```

### structure.txt - DNA encoding scheme

It is a main input file describing DNA encoding structure.

```
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
```

- Line 1: FASTA format data file containing NGS data
- Line 2: minimum allowed length of FASTA sequence. Sequences shorter than this length are ignored.
- Line 3: number of building blocks for cycle 1, 2, 3, and 4
- Line 4: output format. Counts are reported for a combination of two building blocks. Other options are 3BB and 4BB. 
- Line 5: maximum allowed mismatches. Sequences with more mismatches are not counted.
- Line 6-: DNA encoding scheme. 

In the sixth and the following lines (DNA encoding scheme),

- Column 1: x, y, z, and $ are reserved for the building block sequences. 1, 2, 3, 4, and 5 are reserved for constant sequence block.
- Column 2: start position(number, note: the number starts from 1)
- Column 3: end position(number)
- Column 4: filename if 1st column is x, y, z or $ and constant sequence if that is 1-5. 

### bb1.txt - building block information

bb2.txt and bb3.txt are similar files for the 2nd and 3rd building blocks.

```
GCTGCC
GCTTGC
GCTGAC
GCTGAG
GCTTTC
.....
```
