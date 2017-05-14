# decode
Decode and analyze NGS(next generation sequencing) data from DNA-Encoded Library screening 

# input - building block data file
GCTGCC FAa-001 10K DEL 1
GCTTGC FAa-002 10K DEL 1
GCTGAC FAa-003 10K DEL 1
GCTGAG FAa-004 10K DEL 1
(--- omitted ---)

column 1 : sequence tag for a building block
column 2 : name or molecule code for a building block
column 3 ... : any extra information
last column : building cycle number (a same building block can be used with different sequence tags for different building cycles)

# input - fastq file
