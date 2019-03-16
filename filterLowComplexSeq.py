#! /usr/bin/env python

#Constructs an unordered dictionary for fasta files
def fasDict(i_f):
  seqDict = {}
  with open(i_f) as infile:
    for line in infile:
      line = line.strip()
      if line.startswith('>'):
        header = line.split()[0].replace('>', '')
        seqDict[header] = ''
      else:
        seqDict[header] += line

  return seDict

#MAIN 
#Read fasta file from command line
small_seqs = fasDict(sys.argv[1])

for head, seq in small_seqs.items():
  #create a compression score for each small RNA, higher scores == lower complexity 
  score = float(len(zlib.compress(seq)))/len(seq)
  if score >= 0.75:
    print ">{}\n{}".format(head, seq)
