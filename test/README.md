# Small test data-sets

This is a small-ish dataset for developing the pipeline.

The genomes are the first 2 chromosomes from _Aspergillus fumigatus_ isolates A1163 and Af293 and collected from <http://www.aspergillusgenome.org>.

The RNAseq data is from Watkins and Liu et al (2018; doi:10.1099/mgen.0.000154), from _Aspergillus fumigatus_ strains Af293 or CEA10.
Reads were downloaded from SRA archives SRR6163911 and SRR6163776 (BioProject PRJNA399754), which are both replicates of conidia in culture at 16 hours.
Reads were subsampled using [bbsplit](http://seqanswers.com/forums/showthread.php?t=41288) to select reads only in the first 1219878 bp of the chromosome.
