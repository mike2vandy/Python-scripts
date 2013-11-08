Python-scripts
==============

Potentially useful Python scripts for file parsing or general purpose bioinformatics.

seq.file.converter.py
-------------------------

Reads in fasta files and non-interleaved phylip and nexus files.<br>
Will output fasta, interleaved and non-interleaved phylip and nexus files.<br>

Required arguments:<br>
-i < infile ><br>
-inf < file_type > FASTA, NEXUS or PHYLP (all caps, not misspelled).<br>
-outf < file_type > FASTA, NEXUS or PHYLP (all caps, not misspelled).<br>

Optional arguments:<br>
-prot if present will write protein to nexus file. Default: DNA<br>
-int if present will write interleaved nexus and phylip files. Default: sequential


pirnaClusterFinder.py
--------------------------

Will find and annotate piRNA clusters from a .SAM file<br>

Required arguments:
-i < infile.sam >
