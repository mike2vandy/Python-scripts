Python-scripts
==============

Potentially useful Python scripts for file parsing or general purpose bioinformatics.

seq.file.converter.py:

Reads in fasta files and non-interleaved phylip and nexus files.
Will output fasta, interleaved and non-interleaved phylip and nexus files.

Required arguments:
-i <infile>
-inf <file_type> FASTA, NEXUS or PHYLP (all caps, not misspelled).
-outf <file_type> FASTA, NEXUS or PHYLP (all caps, not misspelled).

Optional arguments:
-prot if present will write protein to nexus file. Default: DNA
-int if present will write interleaved nexus and phylip files. Default: sequential
