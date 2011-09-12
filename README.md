##pileupExonFilter.py

###Usage
pileupExonFilter.py bedfile in_pileup > out_pileup

###Function
Remove bases which position is locating outiside of exon ranges.
Exon ranges are specified by bedfile.

###Dependency

python RangeSet.py
##pileupExonFilter.sh

###Usage
pileupExonFilter.sh bedfile in_pileup out_pileup > report

###Function
A Wrapper shellscript of pileupExonFilter.py.

###Dependency
pileupExonFilter.py

##RangeSet.py

###Usage
class file

###Function
Defines RangeSet class. 
Create normalized range list internally.

###Dependency
python

##fastqUnpairedFilter.py

###Usage
fastqUnpairedFilter.py in_seq_1 in_seq_2 out_seq_1 out_seq_2

###Function
Remove Unpaired read from pair fastq sequence files.

###Dependency
python

##fq2bam.sh

###Usage
fq2bam.sh [options] in_seq_1 in_seq_2

###Function
A lazy script to create sai,bam,pileup from paired fastq sequence data.

###Options
####-o PATH
Specify output directory.
If not specified, create out/ directory at where the location of in_seq_1.

####-r PATH
Specify the directory which contais reference files.
If not specified, uses /usr/local/share/doc/hg19/ as default.

####-p
Create pileup file additionaly.

####-a
Do NOT operate adapter clipping for the raw sequence data.

####-q
Do NOT operate qualitry trimming for the raw sequence data.

####-t INT
Chain for fastq_quality_trimmer -t option.

####-l INT
Chain for fastq_quality_trimmer -l option.

####-u
Do NOT remove unpaired reads for the raw sequence data.

####-c INT
Chain for bwa aln -t option.

####-i
Do ill2sanger convert.

###dependency
python fastqUnpairedFilter.py bwa samtools fastx_toolkit
