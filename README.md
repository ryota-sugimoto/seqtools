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

####-r PATH
Specify reference fasta.

####-p
Create pileup file additionaly.

####-a
Do not operate adapter clipping for the raw sequence data.

####-A [ATGC][ATGC]*
Specify apaptor sequence

####-M
Required minimum adaptor length of clipping.

####-q
Do not operate qualitry trimming for the raw sequence data.

####-t INT
fastq_quality_trimmer -t option.

####-l INT
fastq_quality_trimmer -l option.

####-u
Do not remove unpaired reads for the raw sequence data.

####-c INT
Chain for bwa aln -t option.

####-i
Do ill2sanger convert.

####-d
Delete internal files.

###dependency
python fastqUnpairedFilter.py bwa samtools fastx_toolkit
