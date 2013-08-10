#! /usr/bin/env python

import sys, argparse, re

''' 
Written in Python2.7

Code to convert among file types
fasta, nexus, and phylip without BioPython
'''

#help screen
def print_help():
   print "usage: seq.file.converter.py [-h] [-i INFILE] [-inf INFILE_FORMAT] \
  \n 			       [-outf OUTFILE_FORMAT] [-int] [-prot] \
        \n\narguments: \
        \n -h --help             help \
        \n -i INFILE             in file \
	\n -inf INFILE_FORMAT    FASTA, NEXUS, PHYLIP \
        \n -outf OUTFILE_FORMAT  FASTA, NEXUS, PHYLIP \
        \n -int                  Interleave for NEXUS and PHYLIP only [optional] \
        \n -prot                 if protein alignment, add this option "

#read and save a fasta file into memory
def build_fas_dict(f):
  
  seq_dict = {}
  for line in f:
    line = line.strip()
    #Error check. Will save to memory only if '>' is found
    #Else prints error statement
    try:
      if line.startswith('>'):
	header = line.replace('>', '')
	seq_dict[header] = ''
      else:
	seq_dict[header] += line
    except:
      sys.exit('Error: Check file, incorrect format or blank lines\n')
  
  return seq_dict

#will read white space delimited sequential nexus or phylip file
#Cannot handle interleaved input
def build_nex_phy_dict(f):
  
  seq_dict = {}
  reg = re.compile('^[-._\[\]|A-Za-z0-9]+\s+[-A-Za-z?]+$')
  for line in f:
    if reg.match(line):
      line = line.strip()
      fields = line.split()
      header, seq = fields[0], fields[1].replace(' ', '')
      seq_dict[header] = seq
  
  #Error check. if the seq_dict is empty 
  #will return error message
  if seq_dict:
    return seq_dict
  else:
    sys.exit('Error: Check file, incorrect format.\n \
      Can only read in sequential phylip or nexus formats.\n \
      A character may be unsupported.\n')

# function to print a fasta with 70 bp per line    
def print_fasta(seqs):
  
  for seq in sorted(seqs):
    start = 0
    print ">{}".format(seq)
    while start <= len(seqs[seq]):
      print seqs[seq][start:start+70]
      start += 70

#function to print nexus file with options for NEXUS header
#Supports interleaved output
def print_nex(seqs, interleaved, protein):
  
  #Condition statements for output
  if protein:
    prot = 'protein'
  else:
    prot = 'dna'
  if interleaved:
    inter = 'yes'
  else:
    inter = 'no'
  
  count = 0
  length = 0 #used for length of the alignment
  num_seqs = len(seqs) #number of sequences in file
  
  while count <= 1:
    for i in seqs:
      length = len(seqs[i])
      count += 1
  
  #Prints initial nexus block
  print '#NEXUS'
  print 'begin data;'
  print ' dimensions ntax={} nchar={};'.format(num_seqs, length)
  print ' format datatype={} missing=? gap=- interleave={};'.format(prot, inter)
  print 'matrix'
  
  #Print a tab delimited interleaved nexus file for long sequences 
  if interleaved:
    start = 0
    while start <= length:
      for seq in sorted(seqs):
	print "{}\t{}".format(seq, seqs[seq][start:start + 70])
      start += 70
      print ''
  
  #Prints the sequential version
  else:
    for seq in sorted(seqs):
      print "{}\t{}".format(seq, seqs[seq])
  
  #End of nexus file
  print ';'
  print 'end;'

# Function for phylip format Seq starts printing on 11 (12) character
def get_space(header):
  length = len(header)
  space = 11 - length
  count_1 = 0
  sp = ''
  while count_1 < space:
    sp += ' '
    count_1 += 1
  
  return sp

def print_phylip(seqs, interleaved):
  
  #Get length of alignment
  count = 0
  aln_len = 0
  num_seqs = len(seqs)
  while count <= 1:
    for seq in seqs:
      aln_len = len(seqs[seq])
      count += 1
  
  #Prints file header
  print '{}{} {}'.format(' ', num_seqs, aln_len)
  
  #Statement for interleaved file
  #Phylip is the biggest pain to write out
  if interleaved:
    start = 0
    count = 0
    
    while start <= aln_len:
      #writes first block w/names
      if count == 0:
	for seq in sorted(seqs):
	  header = ''
	  if len(seq) <= 10:
	    header = seq
	  else:
	    header = seq[0:9]
	  
	  sp = get_space(header)
	  
	  print "{}{}{}".format(header, sp, seqs[seq][start:start+70])
	  count += 1
	print ''
	start += 70
      
      #prints remainder of sequence w/out sequence name
      if count > 0:
	sp = '           '
	for seq in sorted(seqs):
	  print "{}{}".format(sp, seqs[seq][start:start+70])
	print ''
	start += 70
  
  #Prints a sequential file
  else:
    for seq in seqs:
      header = ''
      if len(seq) <= 10:
	header = seq
      else:
	header = seq[0:9]
      
      sp = get_space(header)
      print "{}{}{}".format(header, sp, seqs[seq])

##MAIN##

#read command line return help if only program is called
if len(sys.argv) == 1:
        sys.exit(print_help())

#arguement parser 
#Pain
parser = argparse.ArgumentParser()

parser.add_argument('-i', action='store', dest='infile',
                    help = 'in file')
parser.add_argument('-inf', action='store', dest='infile_format',
                    help='FASTA, NEXUS, PHYLIP', type=str)
parser.add_argument('-outf', action='store', dest='outfile_format',
                    help='FASTA, NEXUS, PHYLIP', type=str)
parser.add_argument('-int', action='store_true', dest='interleave',
                    help="Interleave for NEXUS and PHYLIP only [optional]")
parser.add_argument('-prot', action='store_true', dest='protein',
		    help='if protein alignment, add this option')
results = parser.parse_args()

#initial dictionary (hash) for sequences
seqs = {}

#Default is sequential DNA output 
protein = False
interleaved = False

#Statement to confirm infile and infile format are in place
if results.infile:
        inf=open(results.infile)
        if results.infile_format == 'NEXUS' or results.infile_format == 'PHYLIP':
	  seqs = build_nex_phy_dict(inf)
	elif results.infile_format == 'FASTA':
	  seqs = build_fas_dict(inf)
	else:
	  sys.exit(print_help())
else:
  sys.exit(print_help())

inf.close()

#Update if input is protein and/or output is interleaved
if results.protein:
  protein = True
if results.interleave:
  interleaved = True
  
#Confirm correct output format is called  
if results.outfile_format == 'NEXUS':
  print_nex(seqs, interleaved, protein)
elif results.outfile_format == 'PHYLIP':
  print_phylip(seqs, interleaved)
elif results.outfile_format == 'FASTA':
  print_fasta(seqs)
else:
  sys.exit(print_help())	
