#! /usr/bin/env python

import sys, argparse

def print_help():
	print "usage: cluster_te.py [-i INFILE] \
	\n\noptional arguments: \
	\n -h           help \
	\n -o OUTFILE   out file [optiional] \
	\n -w WINDOW    window size [3000] \
	\n -s STEP      window step size [500] \
	\n -p PIRNA     minimum numebr of piRNA per window to delcare cluster [30] \
	\n -r RATIO     5' T ratio [0] \
	\n -min MIN     minumum piRNA size [24] \
	\n -max MAX     maximum piRNA size [32] "

def make_windows(min, max, win_size, slide):
	windows = []
	win_start = min
	while win_start <= max:
		win_end = win_start + win_size
		win_frags = [win_start, win_end]
		windows.append(win_frags)
		win_start += slide
	
	return sorted(windows)
	
###MAIN###

#declare parameters	
win_size = 3000
slide = 500
min_t_rat = 0
min_pi = 24
max_pi = 32
cluster_min = 30

#read command line
if len(sys.argv) == 1:
	sys.exit(print_help())

parser = argparse.ArgumentParser()

parser.add_argument('-i', action='store', dest='infile',
                    help = 'in file')
parser.add_argument('-o', action='store', dest='outfile',
                    help='out file [optional]')					
parser.add_argument('-w', action='store', dest='window',
                    help='window size [3000]', type=int)
parser.add_argument('-s', action='store', dest='step',
                    help='window step size [500]', type=int)
parser.add_argument('-p', action='store', dest='pirna',
                    help="minimum number of piRNA per window to declare cluster [30]", type=int)
parser.add_argument('-r', action='store', dest='ratio',
                    help="5' T ratio [0]", type=float)					
parser.add_argument('-min', action='store', dest='min',
                    help="minimum piRNA size [24]", type=int)
parser.add_argument('-max', action='store', dest='max',
                    help="maximum piRNA size [32]", type=int)

results = parser.parse_args()

#update parameters
if results.infile:
	f=open(results.infile)
if results.window:
	win_size = results.window
if results.step:
	slide = results.step
if results.pirna:
	cluster_min = results.pirna
if results.ratio:
	min_t_rat = results.ratio
if results.min:
	min_pi = results.min
if results.max:
	max_pi = results.max

data = {} #contigs are hash keys and list of pirnas on the cotig are values
cluster_hash = {} #final hash of clusters

#working program starts
#read SAM into memory
for line in f:
	try:
		line = line.strip()
		fields = line.split()
		seq_id, orient, contig, start, seq_len =  fields[0], fields[1], fields[2], int(fields[3]), len(fields[9])
		first = list(fields[9])[0]
		if orient != '4' and seq_len >= min_pi and seq_len <= max_pi:
			new_list = [start, first, orient, seq_id]
			if contig not in data:
				data[contig] = []
				data[contig].append(new_list) #stores position, first base, orientation and piRNA id into data hash 
			else:
				data[contig].append(new_list)
	except:
		pass
f.close()

#parse contigs	
for chrm in data:
	data[chrm].sort() #sorts by position on contig: smallest to largest
	min = int(data[chrm][0][0]) #first piRNA on the contig
	max = int(data[chrm][-1][0]) #last piRNA on the contig
	
	#make windows
	windows = make_windows(min, max, win_size, slide)

	#intersect windows with piRNAs
	wins = []
	for window in windows:
		win_start, win_end = window[0], window[1]
		pi_count = 0
		first_nucl = []
		pis_in_window = []
		rm = []
		
		for pi in data[chrm]:
			start, first = int(pi[0]), pi[1]
			if start >= win_start and start <= win_end:
				pis_in_window.append(pi)
				pi_count += 1
				first_nucl.append(first)
			elif start > win_end:
				break
			elif start < win_start:
				rm.append(pi) #remove piRNA found in a previous window to save time
		
		for pi in rm:
			data[chrm].remove(pi)
				
		if len(first_nucl) > 0:
			t_ratio = float(first_nucl.count('T')) / len(first_nucl)
			if pi_count >= cluster_min and t_ratio >= min_t_rat:
				pis_in_window.sort()
				wins.append([win_start, win_end , pis_in_window])
	
	#merge overlapping windows together
	count = 0
	merged_clusters = {}
	index_count = 0
	while count < len(wins):
		try:
			index = chrm + '-' + str(index_count)
			if wins[count][0] >= (wins[count + 1][0] - 1000):
				if index not in merged_clusters:
					merged_clusters[index] = []
					merged_clusters[index] += wins[count][2]
				else:
					merged_clusters[index] += wins[count][2]
				count += 1
			else:
				index_count += 1
				count += 1
		except:
			break
	
	for key in sorted(merged_clusters):
		pirnas = dict((x[3], x) for x in merged_clusters[key]).values() #remove duplicate pirnas found in overlapping windows 
		pirnas.sort(key=lambda x: int(x[0]))
		plus = minus = 0
		start = pirnas[0][0]
		stop = pirnas[-1][0]
		cluster_len = stop - start
		number = len(pirnas)
		density = float(number) / cluster_len
		for pi in pirnas:
			orient = int(pi[2])
			if orient == 0:
				plus += 1
			else:
				minus += 1
			
		cluster_hash[key] = [str(start), str(stop), str(cluster_len), str(number), str(plus), str(minus), str(density)]

if results.outfile:
	out = open(results.outfile, 'w')
	out.write(str('#id\tstart\tstop\tlength\ttotal\tplus_strand\tminus_strand\tdensity\n'))
	for key in sorted(cluster_hash):
		out.write(str(key + '\t' + '\t'.join(cluster_hash[key]) + '\n'))
	out.close()
else:
	print '#id\tstart\tstop\tlength\ttotal\tplus_strand\tminus_strand\tdensity'
	for key in sorted(cluster_hash):
		print key + '\t' + '\t'.join(cluster_hash[key])
