#!/usr/bin/env python

from subprocess import call
import multiprocessing, os, sys
from optparse import OptionParser


usage = '''
getClusterML instructions
Usage: getClusterML.py [options] ref in_prefix

ref: path and filename for reference fasta file
in_prefix: path and prefix for input positions and snp files
output_prefix: path and text for output 

Options
-p, --nprocs  number of cores to use, default 1
-c, --cf   draw ClonalFrame ML trees 
-q, --iqtree	use iqtree rather than phyml
-r, --round     number of dp for rounding of scaled trees, default 1 

Example:

getClusterML.py -p 4 -c ref.fa input_prefix

'''


if __name__=='__main__':
	
	### PROCESS INPUT OPTIONS AND ARGUMENTS ###	
	parser = OptionParser()
	parser.add_option( '-p', '--nprocs', action = 'store', type='int', dest = 'nprocs', default = 1 )
	parser.add_option( "-c", "--cf", action = "store_true", dest = "draw_cf", default = False )
	parser.add_option( "-q", "--iqtree", action = "store_true", dest = "use_iqtree", default = False )
	parser.add_option( "-r", "--round", action = "store", dest = "round_dp", default = 1 )
	opts, args = parser.parse_args()
	
	#check have correct number of arguments
	if len(args) == 2:
		ref_file, input_stem = args
		nprocs = opts.nprocs
		draw_cf = opts.draw_cf
		use_iqtree = opts.use_iqtree
		round_dp = opts.round_dp
	else:
		sys.stdout.write(usage)
		sys.exit(1)
	
	# cluster_ml directory
	if not os.path.isdir('%s/cluster_ml'%input_stem):
		os.mkdir('%s/cluster_ml'%input_stem)
	
	#get list of clusters, comids
	clusterfile = '%s/clusters.txt'%input_stem
	sys.stdout.write('Opening cluster file: %s\n'%clusterfile)
	f = open(clusterfile, 'r')
	bin = f.next()
	clusterDict = dict()
	
	for l in f:
		l = l.strip().split()
		if l[0] in clusterDict.keys():
			clusterDict[l[0]].append(l[1])
		else:
			clusterDict[l[0]] = [l[1]]
	
	f.close()
	
	## use parallel loop to generate ML trees
	
	# run ML tree
	def getML(c):
		#generate ML tree
		outStem = '%s/cluster_ml/cluster_%s'%(input_stem, c)
		inStem = '%s/cluster/cluster_%s'%(input_stem, c)
		tree = 't'
		if draw_cf: tree = 'c'
		iqtree = ''
		if use_iqtree: iqtree = ' -q'
		cmd = '/home/nick/soft/runListCompare/padTree.py -%s%s -d -r %s %s %s_clean_positions.txt %s_clean_snps.fa %s'%(tree, iqtree, round_dp, ref_file, inStem, inStem, outStem)
		sys.stdout.write(cmd+'\n')
		sys.stdout.flush()
		call(cmd.split())
		
		if draw_cf:
			#generate alignment with recombination removed
			cmd = '/home/nick/soft/runListCompare/removeRecombination.py -c %s %s'%(c, input_stem)
			sys.stdout.write(cmd+'\n')
			sys.stdout.flush()
			call(cmd.split())
	
	
	# do get ML for particular sequence in list and cluster list
	def runML(i, cd, nprocs):
		for c in cd[i::nprocs]:
			getML(c)
	
	# get list of clusters
	cd = [int(c) for c in clusterDict.keys() if len(clusterDict[c])>2]
	cd.sort(reverse=True)
	#print clusterDict
	
	procs = []
	
	for i in range(nprocs):
		p = multiprocessing.Process(
				target=runML,
				args=(i, cd, nprocs))
		procs.append(p)
		p.start()
	
	# Wait for all worker processes to finish
	for p in procs:
		p.join()
