#! /usr/bin/python

from Bio import SeqIO
from Bio import Entrez
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast import NCBIXML
from Bio.Alphabet import generic_dna
import numpy as np
from biom.table import Table
import random
import json

Entrez.email = "c.hahn@hull.ac.uk"
import time
import re
import sys, warnings
import argparse
from argparse import RawTextHelpFormatter
import os.path
from os import rename
from collections import defaultdict
import shlex, subprocess

#############################################################################
VERSION="0.1"
informats = {'gb': 'gb', 'genbank': 'gb', 'fasta': 'fasta', 'fa': 'fasta', 'fastq': 'fastq'}
methods = []	#this list will contain the list of methods to be applied to the queries
all_seqs = []
skipped_ref_seqs = [] #defaultdict(list)
references = {}
queries = defaultdict(dict)
seq_info = ['"seqname","accession","tax_id","species_name","is_type"']
records = {}
reference_taxa = {}
taxids = defaultdict(int)
date = time.strftime("%d-%b-%Y").upper()
files_to_barcodes = defaultdict(dict)
global_taxa = defaultdict(dict)
query_count = 0
global_taxids_hit = defaultdict(int)

parser = argparse.ArgumentParser(description='Convert pplacer\'s jplace file to biom format')

parser.add_argument("-p", "--prefix", help="prefix for output files - will also be sample ID in biom files (default=\'pplace\'", metavar="<string>", action="store", default='pplace')

phyloplace_group = parser.add_argument_group('Phylogenetic placement', 'The parameters in this group affect phylogenetic placement')
phyloplace_group.add_argument("-j", "--jplace", help='PATH to jplace file', action='store')
#phyloplace_group.add_argument("--refpkg", help="PATH to refpkg", metavar="<DIR>", action="store")

biom_group = parser.add_argument_group('BIOM OUTPUT','The arguments in this groups affect the output in BIOM format')
#biom_group.add_argument("-o","--output_prefix", help="prefix for BIOM output files (default='metaBEAT')", action="store", default="metaBEAT")
biom_group.add_argument("--mock_meta_data", help="add mock metadata to the samples in the BIOM output", action="store_true")

parser.add_argument("--version", action="version", version='%(prog)s v.'+VERSION)
args = parser.parse_args()

####START OF MAIN PROGRAM
print "\n##### INTERPRETING PPLACER RESULTS #####\n"
fh = open(args.jplace,"r")
jplace = json.load(fh)
fh.close()
#print "\n### find all placements###"
for pp_queries in jplace['placements']:
#	print pp_queries['p']
	#print pp_queries['p'][0] #just consider the first placement
	placement = pp_queries['p'][0][0] #take the first placement for now
	#print "add taxid \'%s\' to global taxid dictionary" %placement
	global_taxids_hit[placement] += len(pp_queries['nm']) #add taxid of placement to global dictionary

#print global_taxids_hit

#fh = open(args.refpkg+"/CONTENTS.json","r")
#refpkg_content = json.load(fh)
#taxonomy = args.refpkg+'/'+refpkg_content['files']['taxonomy']
#fh.close()

data_to_biom = []
observ_ids = []
for taxid in sorted(global_taxids_hit):
#	print taxid
        data_to_biom.append(global_taxids_hit[taxid])

data = np.asarray(data_to_biom).reshape(len(data_to_biom),1)
#print data

sample_id = [args.prefix]
#print sample_id

sample_metadata=[]
for q in sample_id:
	temp={}
	temp['method'] = 'pplacer'
	if args.mock_meta_data:
		r="%.1f" %random.uniform(20.0,25.0)
		temp['temperature'] = "%.1f C" %float(r)
		r="%.1f" %random.uniform(10.0,15.0)
		temp['depth'] = "%.1f m" %float(r)
		treatments = ['A','B']
		temp['treatment'] = random.choice(treatments)

	sample_metadata.append(temp)
#print sample_metadata

observation_metadata=[]

Taxonomy=defaultdict(dict)
syn = {'kingdom': 'k__', 'phylum': 'p__', 'class': 'c__', 'order': 'o__', 'family': 'f__', 'genus':'g__', 'species': 's__'}
levels = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
#print "These taxids were identified:\n%s" %global_taxids_hit
for tid in sorted(global_taxids_hit.keys()):
#	print "fetch taxonomy for taxid: %s" %tid
	handle = Entrez.efetch(db="Taxonomy", id=tid)      #search the taxonomy database for the taxon by organism name
	taxon = Entrez.read(handle)
#	print taxon[0]['Lineage']
#	print taxon[0]['ScientificName']
	ind_taxonomy = []
	observ_ids.append(taxon[0]['ScientificName'])
	
	for lev in levels:
#		print lev
		for i in range(len(taxon[0]['LineageEx'])):
			if taxon[0]['LineageEx'][i]['Rank'] == lev:
				string = taxon[0]['LineageEx'][i]['ScientificName']
				string = string.replace(' ', '_')
#				print "%s%s" %(syn[lev], string)
				ind_taxonomy.append('%s%s' %(syn[lev], string))
				break
			
	if len(ind_taxonomy) < len(levels):
#		print taxon[0]['Rank']
		if taxon[0]['Rank'] in levels:
			index = levels.index(taxon[0]['Rank'])
#			print "index: %i" %index
			ind_taxonomy.append('%s%s' %(syn[levels[index]], taxon[0]['ScientificName']))
       	

#	print ind_taxonomy			
	
	Taxonomy[taxon[0]['ScientificName']]['taxonomy'] = ind_taxonomy

#	print "Taxonomy: %s" %Taxonomy
		
for taxon in observ_ids:
#	print taxon
#	print Taxonomy[taxon]
	observation_metadata.append(Taxonomy[taxon])

#print "observation metadata:\n%s" %observation_metadata
#print len(observation_metadata)

table = Table(data, observ_ids, sample_id, observation_metadata, sample_metadata, table_id='Example Table')
print table

out=open(args.prefix+".biom","w")
table.to_json('pplacer converted by jplace_to_biom.py v.'+VERSION, direct_io=out)
out.close()

out=open(args.prefix+".tsv","w")
out.write(table.to_tsv(header_key='taxonomy', header_value='taxomomy')) #to_json('generaged by test', direct_io=out)
out.close()

print "\n##### DONE! #####\n"
