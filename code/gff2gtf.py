#!/usr/bin/env python
# gff2gtf.py
# convert gff gene annotation file to gtf
# Author: Tuo Zhang
# Version: 1.0
# NEW: first version
# 

import sys
from re import search

# functions
def splitAttributes(tAttrs):
	# split attributes
	tatrdic = {}
	for tx in tAttrs.split(";"):# each attribute
		tpat = search("(.*)=(.*)", tx)
		if tpat:
			tid,tval = tpat.groups()
			if tid not in tatrdic:
				tatrdic[tid] = [tval]
			else:
				tatrdic[tid].append(tval)
	return tatrdic

def getAttribute(tAttrs,tid):
	# given a attribute id, return its value if exists
	if tid in tAttrs:
		if len(tAttrs[tid]) > 1:
			print "Error: attribute id %s presents more than once in a single entry: %s"%(tid,";".join(tAttrs))
			sys.exit(2)
		return tAttrs[tid][0]
	else:
		return ""

# main
if len(sys.argv) != 4:
	print "Usage: gff2gtf.py in.gff out.gtf match.table"
	sys.exit(1)

infile = sys.argv[1]
outfile = sys.argv[2]
matchfile = sys.argv[3]

# read in gff file
fin = file(infile,'r')
gff = [[y.strip() for y in x.split("\t")] for x in fin.readlines() if not search("^#", x)]
fin.close()
print len(gff),"entries loaded from input gff file."

# build dictionary based on ID
items = {}
for k,elements in enumerate(gff):# each entry
	feature = elements[2]
	atrs = splitAttributes(elements[8])
	# get attributes:
	eid = getAttribute(atrs, "ID")
	parent = getAttribute(atrs, "Parent")
	name = getAttribute(atrs, "Name")
	gbkey = getAttribute(atrs, "gbkey")
	gene = getAttribute(atrs, "gene")
	biotype = getAttribute(atrs, "gene_biotype")
	# double check ID
	if eid == "":
		print "Error: cannot find 'ID' for %dth entry: %s"%(k+1,"\t".join(elements))
		sys.exit(1)
	# add to dictionary
	if eid not in items:
		items[eid] = (parent, name, gbkey, gene, biotype)
	else:
		print "Warning: non-unique 'ID' detected for %dth entry: %s"%(k+1, eid)
print len(items),"unique ID detected."

# process exon entries
tsps = []
mtable = {}
fout = file(outfile,'w')
for k,elements in enumerate(gff):# each entry
	feature = elements[2]
	atrs = splitAttributes(elements[8])
	eid = getAttribute(atrs, "ID")
	# only process 'exon' entry
	if feature != "exon":
		continue
	# get transcript_id and gene_id
	# get parent entry
	parent, name, gbkey, gene, biotype = items[eid]
	# get ancestry entry
	tp = ancestry = parent
	while tp != "":
		ancestry = tp
		tp = items[tp][0]
	# gene_id <-- 'ID' of ancestry entry
	gid = ancestry
	# transcript_id <-- 'ID' of parent entry
	tspid = parent
	# gene_biotype <-- 'gene_biotype' of ancestry entry
	gbiotype = items[ancestry][4]
	# gene_name <-- 'gene' of ancestry entry
	gname = items[ancestry][3]
	# transcript_biotype <-- 'gbkey' of parent entry
	tbiotype = items[parent][2]
	# aname <-- 'name' of ancestry entry
	aname = items[ancestry][1]
	# add transcript info to match table
	if tspid not in mtable:
		tsps.append(tspid)
		mtable[tspid] = [tspid, gid, aname, gbiotype, tbiotype, gname]
	##print "; ".join([parent, name, gbkey, gene, biotype])
	##print "; ".join([tp, ancestry])
	##print "; ".join([gid, tspid, gname, gbiotype, tbiotype])
	##break
	# output the first 8 columns
	fout.write("\t".join(elements[:8])+"\t")
	# add attributes
	fout.write('gene_id "%s"; transcript_id "%s"; gene_name "%s"; gene_biotype "%s"; transcript_biotype "%s";\n'%(gid,tspid,aname,gbiotype,tbiotype))
fout.close()

# output match table
fmat = file(matchfile,'w')
fmat.write("transcript_id\tgene_id\tgene_name\tgene_biotype\ttranscript_biotype\tname\n")
for tspid in tsps:# each transcript
	fmat.write("\t".join(mtable[tspid])+"\n")
fmat.close()

print "Complete!"

