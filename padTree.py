#!/usr/bin/env python

#############################################################################
# padTree.py
# takes a reference, position file, snps file as input, and an output prefix
# generates padded file as output
# generates an ML tree
# check 13 Oct 2015
#############################################################################


from optparse import OptionParser
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
import os, sys, re
import subprocess, shlex
import gzip

usage = '''
Dependencies - Fasta2Phylip.pl, phyml, ClonalFrameML in bin

padTree instructions
Usage: padTree.py [options] reference_file position_file snps_file output_prefix

Example:

bin/padTree.py -c -d -r 1 /home/davide/ref/R00000003.fasta cluster20/cluster_200_clean_positions.txt cluster20/cluster_200_clean_snps.fa cluster20_ml/cluster_200

Options
-t, --tree		draw ML tree
-c, --clonalframe draw ClonalFrameML tree (also creates ML tree)
-d, --delete 	remove padded files
-s, --search	search strategy for PhyML [default BEST; NNI quicker, SPR more accurate, BEST best of the 2]
-r, --round     number of dp for rounding of scaled trees, default 1
-m, --mask		mask file of positions to mask, counting from 1, default no mask
-p, --nophy		don't write phylip file
-q, --iqtree	use iqtree rather than phyml

'''


if __name__=="__main__":
	
	
	### PROCESS INPUT OPTIONS AND ARGUMENTS ###
	parser = OptionParser()
	parser.add_option( "-t", "--tree", action = "store_true", dest = "draw_tree", default = False )
	parser.add_option( "-c", "--clonalframe", action = "store_true", dest = "draw_cf", default = False )
	parser.add_option( "-d", "--delete", action = "store_true", dest = "remove_pad", default = False )
	parser.add_option( "-s", "--search", action = "store", dest = "search_method", default = "BEST" )
	parser.add_option( "-r", "--round", action = "store", type='int', dest = "round_dp", default = 1 )
	parser.add_option( "-m", "--mask", action = "store", dest = "maskfile", default = "" )
	parser.add_option( "-p", "--nophy", action = "store_true", dest = "nophy", default = False )
	parser.add_option( "-q", "--iqtree", action = "store_true", dest = "use_iqtree", default = False )
	
	opts, args = parser.parse_args()
	
	#check have correct number of arguments
	if len(args) == 4:
		ref_file, pos_file, snps_file, out_prefix = args
	else:
		sys.stdout.write(usage)
		sys.exit(1)
	
	#save options	
	draw_tree = opts.draw_tree
	draw_cf = opts.draw_cf
	use_iqtree = opts.use_iqtree
	remove_pad = opts.remove_pad
	search_method = opts.search_method
	round_dp = opts.round_dp
	maskfile = opts.maskfile
	nophy = opts.nophy
	
	
	### READ IN FILES ###
	# read in reference
	if os.path.exists(ref_file):
		#refseq = SeqIO.read( ref_file, "fasta" )
		fa = SeqIO.parse( gzip.open(ref_file, 'rt'), "fasta" )
		fa_seq = Seq("".join([s.seq._data for s in fa]), generic_dna)
		fa_id = "ref"
		refseq = SeqRecord(fa_seq, id=fa_id)
		
	else:
		sys.stdout.write('Reference path incorrect, please check and try again.\n')
		sys.exit(1)
	reflen = len(refseq)
	sys.stdout.write('Read in reference\n')
	
	# read in positions file
	if os.path.exists(pos_file):
		r = open(pos_file,'r')
		pos = [int(l)-1 for l in r if l!='\n'] #convert back to being numbered from zero
		r.close()
	else:
		sys.stdout.write('Position path incorrect, please check and try again.\n')
		sys.exit(1)
	sys.stdout.write('Read in positions\n')
	
	# read in snps file
	if os.path.exists(snps_file):
		mfa = SeqIO.parse(snps_file, 'fasta')
		seqlist = [s for s in mfa]
	else:
		sys.stdout.write('SNPs path incorrect, please check and try again.\n')
		sys.exit(1)
	sys.stdout.write('Read in SNPs\n')
	
	# read in mask file if exists
	use_pos = [1]*len(pos)
	if maskfile:
		# read contents of mask file - counts from 1
		m = open(maskfile, 'r')
		mask_pos = [(int(l.strip().split()[0]), int(l.strip().split()[1])) for l in m]
		m.close()
		for i, p in enumerate(pos):
			for j in mask_pos:
				if p+1 >= j[0] and p+1 <= j[1]:
					use_pos[i] = 0
	
	### check length ###
	if len(seqlist) == 1:
		sys.stdout.write('Only 1 sequence no tree to draw.\n')
	elif len(seqlist) == 2:
		tree =  '(%s:0, %s:%s);'%(seqlist[0].id, seqlist[1].id, len(pos))
		scaled_tree = '%s_phyml_tree_scaled.tree'%out_prefix
		with open (scaled_tree,'w') as w:
			w.write(tree)
		w.close()
		sys.stdout.write('Pairwise tree written in newick format.\n')
		sys.stdout.write('Done.\n')
	else:
		### PAD SEQUENCES FILE ###
		# create a dictionary containing a list of each of the positions in the reference
		new_seq = dict()
		for s in seqlist:
			new_seq[s.id] = [i for i in refseq.seq._data]
		
		#update dictionary entries with the value of each position
		for i,b in enumerate(pos):
			for s in seqlist:
				if use_pos[i]:
					new_seq[s.id][b] = s.seq._data[i]
		
		#update the sequences to be the full padded sequences
		for s in seqlist:
			s.seq._data = ''.join(new_seq[s.id])
			s.description = ''
		
		#write the padded sequences file
		out_fa = out_prefix+'_padded.fa'
		out_phy = out_prefix+'_padded.phy'
		SeqIO.write(seqlist, out_fa, 'fasta')
		sys.stdout.write('Fasta file written\n')
		sys.stdout.flush()
                print(seqlist)
		SeqIO.write(seqlist, out_phy, 'phylip')
		sys.stdout.write('Phylip file written\n')
		sys.stdout.flush()
		
		#del the variables used for padding
		del seqlist
		del new_seq
		
		#sys.stdout.write('Deciding if to use phy\n')
		#sys.stdout.flush()
		#nophy = True
		#if not nophy:
		#	sys.stdout.write('Making phy\n')
		#	sys.stdout.flush()
		#	#convert to phy format using perl (faster than Biopython phylip writer)
		#	cmd = 'perl /users/bag/deyre/bin/Fasta2Phylip.pl %s %s'%(out_fa, out_phy)
		#	sys.stdout.write(cmd+'\n')
		#	sys.stdout.flush()
		#	bin = commands.getstatusoutput(cmd)
		#	sys.stdout.write('%s\t%s\n'%(bin[0], bin[1]))
		#	sys.stdout.flush()
		#	sys.stdout.write('Phy file written\n')
		
		
		### RUN PhyML TO DRAW TREE IF SPECIFIED ###
		
		if draw_tree or draw_cf:
			if not use_iqtree:
				sys.stdout.write('Starting PhyML tree\n')
				#run PhyML
				cmd = 'phyml -d nt -m GTR -s %s -b 0 -v 0.0 -c 4 -a e -i %s'%(search_method, out_phy)
				sys.stdout.write('%s\n'%cmd)
				sys.stdout.flush()
				subprocess.check_call(shlex.split(cmd))
				
				#scale output tree
				unscaled_tree = '%s_phyml_tree.txt'%out_phy
				scaled_tree = '%s_phyml_tree_scaled.tree'%out_prefix
			else:
				#run iqtree
				cmd = 'iqtree -s %s -m GTR+G -blmin 0.00000001 -t PARS'%out_phy
				sys.stdout.write('%s\n'%cmd)
				sys.stdout.flush()
				subprocess.check_call(shlex.split(cmd))
				
				#scale output tree
				unscaled_tree = '%s.treefile'%out_phy
				scaled_tree = '%s_iqtree_scaled.tree'%out_prefix
			
			r = open(unscaled_tree, 'r')
			string = [i.split('\n')[0] for i in r][0]
			lengths = set([i.split(')')[0] for i in [i.split(',')[0] for i in [i for i in string.split(':') if not i.startswith('(')]]])
			
			for l in lengths:
				string = string.replace(l, str(round(float(l)*reflen,round_dp)))
			
			with open(scaled_tree, 'w') as w:
				w.write(string)
			w.close()
			
		### DRAW CLONALFRAME ML TREE IF REQUIRED ###
		
		if draw_cf:
			sys.stdout.write('Starting Clonal Frame ML tree\n')
			cmd = 'ClonalFrameML %s %s %s_cf'%(unscaled_tree, out_fa, out_prefix)
			sys.stdout.write('%s\n'%cmd)
			sys.stdout.flush()
			subprocess.check_call(shlex.split(cmd))
			
			#scale output tree
			unscaled_tree = '%s_cf.labelled_tree.newick'%out_prefix
			scaled_tree = '%s_cf_scaled.tree'%out_prefix
			r = open(unscaled_tree, 'r')
			
			string = [i.split('\n')[0] for i in r][0]
			
			#remove node labels
			string = re.sub('NODE_[0-9]{0,3}','',string)
			lengths = set([i.split(')')[0] for i in [i.split(',')[0] for i in [i for i in string.split(':') if not i.startswith('(')]]])
			lengths = [l for l in lengths if l!='0'and l!='1e-07']
			
			for l in lengths:
				string = string.replace(l, str(round((float(l)*(reflen+10000))-1,round_dp)))
			
			zero = '0'
			if round_dp > 0: zero = '0.'+round_dp*'0'
			string = string.replace('1e-07', zero)
			
			with open(scaled_tree, 'w') as w:
				w.write(string)
			
			w.close()
		
		
		### DELETE PADDED FILES IF SPECIFIED ###
		
		if remove_pad:
			#clean up - remove large padded files
			sys.stdout.write('Removing padded files\n')
			os.remove(out_fa)
			os.remove(out_phy)
		
		sys.stdout.write('Done.\n')
