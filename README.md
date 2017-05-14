# decode
Decode and analyze NGS(next generation sequencing) data from DNA-Encoded Library screening 

## usage
<pre>
python decode.py -b 10k-DEL-BB.txt testbig.fastq
</pre>

## building block data file format
<pre>
GCTGCC FAa-001 10K DEL 1
GCTTGC FAa-002 10K DEL 1
GCTGAC FAa-003 10K DEL 1
GCTGAG FAa-004 10K DEL 1
(--- omitted ---)
</pre>
- column 1 : sequence tag for a building block
- column 2 : name or molecule code for a building block
- column 3 ... : any extra information
- last column : building cycle number (a same building block can be used with different sequence tags for different building cycles)

## fastq file format
