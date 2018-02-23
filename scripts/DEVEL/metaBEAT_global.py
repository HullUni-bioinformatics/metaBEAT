#! /usr/bin/python

from Bio import SeqIO
from Bio import Entrez
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
#from Bio.Alphabet import IUPAC
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast import NCBIXML
from Bio.Alphabet import generic_dna
import numpy as np
from biom.table import Table
import glob
import random
import json
import gzip
from itertools import product

Entrez.email = "" 
import time
from datetime import datetime
import re
import sys, warnings
import argparse
from argparse import RawTextHelpFormatter
import os.path
from os import rename
from collections import defaultdict
import shlex, subprocess
import shutil


#############################################################################
VERSION="0.97.10-global"
DESCRIPTION="metaBEAT - metaBarcoding and Environmental DNA Analyses tool\nversion: v."+VERSION
informats = {'gb': 'gb', 'genbank': 'gb', 'fasta': 'fasta', 'fa': 'fasta', 'fastq': 'fastq', 'uc':'uc'}
methods = []	#this list will contain the list of methods to be applied to the queries
all_seqs = []
skipped_ref_seqs = [] #defaultdict(list)
references = {}
queries = defaultdict(dict)
seq_info = ['"seqname","accession","tax_id","species_name","is_type"']
records = {}
reference_taxa = {}
#taxids = defaultdict(int)
denovo_taxa = {}
date = time.strftime("%d-%b-%Y").upper()
denovo_count = 1
files_to_barcodes = defaultdict(dict)
global_taxa = defaultdict(dict)
global_clust = defaultdict(dict)
query_count = 0
global_taxids_hit = []
metadata = defaultdict(dict)
read_stats = defaultdict(dict)
read_metrics = ['total', 'trimmed-total', 'trimmed-pe', 'trimmed-orphans', 'merged', 'cluster_thres', 'clusters_total', 'clusters_min_cov', 'cluster_above_thres', 'queries']
#primer_clip_string = ""
primer_versions = []
bl_db_extensions = ["nin", "nsq", "nhr"]
blast_dict = defaultdict(dict) #this will hold the results from parsing the BLAST output
gb_to_taxid_dict = {}
taxid_list = []
tax_dict = {}
BIOM_tables_per_method = {}
global_cluster_counts = {}
global_cluster_reads = defaultdict(list)
global_cluster_count = 0
read_counts_out = ""
queries_to_rc = []

parser = argparse.ArgumentParser(description=DESCRIPTION, prog='metaBEAT.py')
#usage = "%prog [options] REFlist"
#parser = argparse.ArgumentParser(usage='%(prog)s [options] REFlist', formatter_class=RawTextHelpFormatter)

parser.add_argument("-Q", "--querylist", help="file containing a list of query files", metavar="<FILE>", action="store")
parser.add_argument("-B", "--BIOM_input", help="OTU table in BIOM format", metavar="<FILE>", action="store")
parser.add_argument("--g_queries", help="fasta file containing query sequences (in combination with '-B' sequence headers are expected to match)", metavar="<FILE>", action="store")
parser.add_argument("-v", "--verbose", help="turn verbose output on", action="store_true")
parser.add_argument("-s", "--seqinfo", help="write out seq_info.csv file", action="store_true")
parser.add_argument("-f", "--fasta", help="write out ref.fasta file", action="store_true")
parser.add_argument("-p", "--pplace", help="perform phylogenetic placement", action="store_true")
parser.add_argument("-k", "--kraken", help="perform phylogenetic placement", action="store_true")
parser.add_argument("-t", "--taxids", help="write out taxid.txt file", action="store_true")
parser.add_argument("-b", "--blast", help="compile local blast db and blast queries", action="store_true")
parser.add_argument("-m", "--marker", help="marker ID (default: marker)", metavar="<string>", action="store", default="marker")
parser.add_argument("-n", "--n_threads", help="Number of threads (default: 1)", type=int, metavar="<INT>", action="store", default="1")
parser.add_argument("-E", "--extract_centroid_reads", help="extract centroid reads to files", action="store_true")
parser.add_argument("-e", "--extract_all_reads", help="extract reads to files", action="store_true")
parser.add_argument("--read_stats_off", help="ommit writing read stats to file", action="store_true")
query_group = parser.add_argument_group('Query preprocessing', 'The parameters in this group affect how the query sequences are processed')
query_group.add_argument("--PCR_primer", help='PCR primers (provided in fasta file) to be clipped from reads', metavar="<FILE>", action="store")
query_group.add_argument("--bc_dist", help="Number of mismatches allowed in barcode sequences", metavar="<INT>", type=int, action="store")
query_group.add_argument("--trim_adapter", help="trim adapters provided in file", metavar="<FILE>", action="store")
query_group.add_argument("--trim_qual", help="minimum phred quality score (default: 30)", metavar="<INT>", type=int, action="store", default=30)
query_group.add_argument("--phred", help="phred quality score offset - 33 or 64 (default: 33)", metavar="<INT>", type=int, action="store", default=33)
query_group.add_argument("--trim_window", help="sliding window size (default: 5) for trimming; if average quality drops below the specified minimum quality all subsequent bases are removed from the reads", metavar="<INT>", type=int, action="store", default=5)
query_group.add_argument("--read_crop", help="Crop reads to this length if they are longer than that (default: off)", metavar="<INT>", type=int, action="store", default=0)
query_group.add_argument("--trim_minlength", help="minimum length of reads to be retained after trimming (default: 50)", metavar="<INT>", type=int, action="store", default=50)
query_group.add_argument("--merge", help="attempt to merge paired-end reads", action="store_true")
query_group.add_argument("--product_length", help="estimated length of PCR product (specifying this option increases merging efficiency)", metavar="<INT>", type=int, action="store")#, default=100)
query_group.add_argument("--merged_only", help="only process successfully merged read-pairs", action="store_true")
query_group.add_argument("--forward_only", help="only process sequences that contain forward reads (i.e. unmerged forward reads and merged reads)", action="store_true")
query_group.add_argument("--length_filter", help="only process reads, which are within +/- 10 percent of this length", metavar="<INT>", type=int, action="store")
query_group.add_argument("--length_deviation", help="allowed deviation (in percent) from length specified by --length_filter (default=0.1)", metavar="<FLOAT>", type=float, action="store", default=0.1)
reference_group = parser.add_argument_group('Reference', 'The parameters in this group affect the reference to be used in the analyses')
reference_group.add_argument("-R", "--REFlist", help="file containing a list of files to be used as reference sequences", metavar="<FILE>", action="store")
reference_group.add_argument("--gb_out", help="output the corrected gb file", metavar="<FILE>", action="store", default="")
reference_group.add_argument("--rec_check", help="check records to be used as reference", action="store_true")
reference_group.add_argument("--gb_to_taxid", help="comma delimited file containing 'gb accession,taxid' for a list of taxa", metavar="<FILE>", action="store", default=os.getcwd()+"/gb_to_taxid.csv")
cluster_group = parser.add_argument_group('Query clustering options', 'The parameters in this group affect read clustering')
cluster_group.add_argument("--cluster", help="perform clustering of query sequences using vsearch", action="store_true")
cluster_group.add_argument("--clust_match", help="identity threshold for clustering in percent (default: 1)", type=float, metavar="<FLOAT>", action="store", default="1")
cluster_group.add_argument("--clust_cov", help="minimum number of records in cluster (default: 1)", type=int, metavar="<INT>", action="store", default="1")
blast_group = parser.add_argument_group('BLAST search', 'The parameters in this group affect BLAST search and BLAST based taxonomic assignment')
blast_group.add_argument("--blast_db", help="path to precompiled blast database", metavar="<PATH>", action="store", default="")
blast_group.add_argument("--blast_xml", help="path to Blast result in xml format", metavar="<PATH>", action="store", default="")
blast_group.add_argument("--update_taxonomy", help="Download/update taxonomy database. Database will be called 'taxonomy.db' and will be compiled in the same location as the metaBEAT.py script.", action="store_true", default=False)
blast_group.add_argument("--taxonomy_db", help="taxonomy database file location. In case it's not the default, which is 'taxonomy.db' in the same directory as the metaBEAT.py script.", metavar="<FILE>", action="store")
#blast_group.add_argument("--www", help="perform online BLAST search against nt database", action="store_true")
blast_group.add_argument("--min_ident", help="minimum identity threshold in percent (default: 0.80)", type=float, metavar="<FLOAT>", action="store", default="0.80")
blast_group.add_argument("--min_ali_length", help="minimum alignment length in percent of total query length (default: 0.95)", type=float, metavar="<FLOAT>", action="store", default="0.95")
blast_group.add_argument("--bitscore_skim_LCA", help="Only BLAST hits with bitscores differing by less than this factor from the top hit (bitscore skim window) will be considered for LCA (0-1; default: 0.1)", type=float, metavar="<FLOAT>", action="store", default="0.1")
blast_group.add_argument("--bitscore_skim_adjust_off", help="Per default a 100%% identity BLAST top hit across the minimum alignment length triggers an adjustment of the bitscore skim window to '0', i.e. only hits with bitscores as good as the top hit are considered for LCA. This flag switches this behaviour off.", action="store_true")
blast_group.add_argument("--min_bit", help="minimum bitscore (default: 80)", type=int, metavar="<INT>", action="store", default="80")
phyloplace_group = parser.add_argument_group('Phylogenetic placement', 'The parameters in this group affect phylogenetic placement')
phyloplace_group.add_argument("--refpkg", help="PATH to refpkg for pplacer", metavar="<DIR>", action="store")
phyloplace_group.add_argument("--jplace", help="phylogenetic placement result from prefious pplacer run in *.jplace format", metavar="<FILE>", action="store")
kraken_group = parser.add_argument_group('Kraken', 'The parameters in this group affect taxonomic assignment using Kraken')
kraken_group.add_argument("--jellyfish_hash_size", help="jellyfish hash size to control memory usage during kraken database building. A table size of '6400M' will require ~44G of RAM. '2700M' -> 20G RAM. ", type=str, default=False, metavar="<STR>", action="store")
kraken_group.add_argument("--kraken_db", help="PATH to a Kraken database", metavar="<DIR>", action="store")
kraken_group.add_argument("--kraken_score_threshold", help="minimum proportion of k-mers to support assignment (0-1; default: 0)", metavar="<FLOAT>", default=0, type=float, action="store")
kraken_group.add_argument("--rm_kraken_db", help="Remove Kraken database after successful completion", action="store_true")

biom_group = parser.add_argument_group('BIOM OUTPUT','The arguments in this groups affect the output in BIOM format')
biom_group.add_argument("-o","--output_prefix", help="prefix for BIOM output files (default='metaBEAT')", action="store", default="metaBEAT")
biom_group.add_argument("--metadata", help="comma delimited file containing metadata (optional)", action="store")
#biom_group.add_argument("--mock_meta_data", help="add mock metadata to the samples in the BIOM output", action="store_true")

Entrez_group = parser.add_argument_group('Entrez identification','metaBEAT is querying the NCBI Entrez databases, please provide an email address for identification')
Entrez_group.add_argument("-@", "--email", help='provide your email address for identification to NCBI', metavar='<email-address>', action="store", default="")

parser.add_argument("--version", action="version", version=VERSION) #'%(prog)s v.'+VERSION)
args = parser.parse_args()


if len(sys.argv) < 2:	#if the script is called without any arguments display the usage
    parser.print_usage()
    sys.exit(1)


###FUNCTIONS

## Funtions used in main script
def parse_BIOM_denovo(table):

    from biom.table import Table
    import json
	#Future additions?
	#check if there is a 'uc' file - if yes add queries[queryID]['format'] = 'uc' -> this will then check if the '*_queries.fasta' file is there.
	#if not check if there is *_queries.fasta file for the sample, if yes add queries[queryID]['format'] = 'fasta' and queries[queryID]['files'] = *_queries.fasta
    
    with open(table) as data_file:    
        data = json.load(data_file)
    t = Table.from_json(data)
    
    return t
    
def add_taxonomy_to_biom(per_tax_level_clusters, per_tax_level_trees, biom_in, method=False):

    	from biom.table import Table
	
	dictionary = {}
	levels = ['nohit', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
	biom_out = biom_in.copy()
	otu_order = []
	trees = []

	for level in reversed(levels):
		if level in per_tax_level_clusters:
#			print "level: %s" %level
			trees_ordered_by_level = []
			for tid in per_tax_level_clusters[level]:
#				print "\ttid: %s" %tid
				for i in range(len(per_tax_level_clusters[level][tid])):
					#add dummy to list
					if tid == 'nohit':
						dictionary[per_tax_level_clusters[level][tid][i]] = {'taxonomy': per_tax_level_trees[level]['unassigned']}
#						print "\t\t%s\t%s" %(per_tax_level_clusters[level][tid][i], per_tax_level_trees[level]['unassigned'])
						trees_ordered_by_level.append("|".join(per_tax_level_trees[level]['unassigned']))
					else:
						dictionary[per_tax_level_clusters[level][tid][i]] = {'taxonomy': per_tax_level_trees[level][tid]}
#						print "\t\t%s\t%s" %(per_tax_level_clusters[level][tid][i], per_tax_level_trees[level][tid])
						trees_ordered_by_level.append("|".join(per_tax_level_trees[level][tid]))

#			print "ORDER:"
			trees_ordered_by_level = sorted(list(set(trees_ordered_by_level)))
			trees.extend(trees_ordered_by_level)
			for tree in trees_ordered_by_level:
#				print "\t%s" %tree
				if level == 'nohit':
					for otu in sorted(per_tax_level_clusters[level]['nohit']):
#                                        	print "\t\t%s" %otu
						otu_order.append(otu)
				else:
					for tid in per_tax_level_clusters[level]:
						if "|".join(per_tax_level_trees[level][tid]) == tree:
							for otu in sorted(per_tax_level_clusters[level][tid]):
#								print "\t\t%s" %otu
								otu_order.append(otu)
						 

#	print "\n### OTU table with taxonomy ###\n"
	biom_out.add_metadata(md=dictionary, axis='observation')
#	print biom_out.to_tsv(header_key='taxonomy', header_value='taxomomy')
#	print "\n\n"


	#sorting by taxonomy
#	print "\n### sorted OTU table with taxonomy ###\n"
	sorted_biom_out = biom_out.sort_order(otu_order, axis='observation')
#	print sorted_biom_out.to_tsv(header_key='taxonomy', header_value='taxomomy')
#	print "\n\n"

	#rename samples - add method to sample id and metadata
	if method:
		rename = {}
		sample_metadata = {}
		for s in sorted_biom_out.ids(axis='sample'):
			rename[s] = s+'.'+method
			sample_metadata[s] = {}
			sample_metadata[s]['method'] = method
		sorted_biom_out.add_metadata(sample_metadata, axis='sample')
		sorted_biom_out.update_ids(rename, axis='sample')

	return sorted_biom_out

def collapse_biom_by_taxonomy(in_table):

    	from biom.table import Table

	bin_f = lambda id_, x: x['taxonomy'][-1].split("__")[1]
	biom_out_collapsed = in_table.collapse(bin_f, norm=False, min_group_size=1, axis='observation')
#	for taxon in biom_out_collapsed.ids(axis='observation'):
		
#	print biom_out_collapsed.to_tsv(header_key='taxonomy', header_value='taxomomy')
#	print "\n\n"

	#restore original order and taxonomy metadata
	trees = [] 
#	print in_table.metadata(axis='observation')
	for otu in in_table.ids(axis='observation'):
#		print otu
		if not "|".join(in_table.metadata(otu, axis='observation')['taxonomy']) in trees:
			trees.append("|".join(in_table.metadata(otu, axis='observation')['taxonomy']))

#	print trees
	collapsed_order = []
	collapsed_meta = {}
	for t in trees:
		collapsed_order.append(t.split("|")[-1].split("__")[-1])
		collapsed_meta[t.split("|")[-1].split("__")[-1]] = {'taxonomy': t.split("|")}
	biom_out_collapsed.add_metadata(md=collapsed_meta, axis='observation')
	biom_out_collapsed_sorted = biom_out_collapsed.sort_order(collapsed_order, axis='observation')
#	print biom_out_collapsed_sorted.to_tsv(header_key='taxonomy', header_value='taxomomy')
#	print "\n\n"

	return biom_out_collapsed_sorted

def pa_and_collapse_by_taxonomy(in_table):

    	from biom.table import Table
	
	bin_f = lambda id_, x: x['taxonomy'][-1].split("__")[1]

#	print "\n### sorted presence-absence OTU table with taxonomy ###\n"
	pa_sorted_biom = in_table.pa(inplace=False)	
	pa_sorted_biom_collapsed = pa_sorted_biom.collapse(bin_f, norm=False, min_group_size=1, axis='observation')

	#restore original order and taxonomy metadata
	trees = [] 
#	print in_table.metadata(axis='observation')
	for otu in in_table.ids(axis='observation'):
#		print otu
		if not "|".join(in_table.metadata(otu, axis='observation')['taxonomy']) in trees:
			trees.append("|".join(in_table.metadata(otu, axis='observation')['taxonomy']))

	collapsed_order = []
	collapsed_meta = {}
	for t in trees:
		collapsed_order.append(t.split("|")[-1].split("__")[-1])
		collapsed_meta[t.split("|")[-1].split("__")[-1]] = {'taxonomy': t.split("|")}

	pa_sorted_biom_collapsed.add_metadata(md=collapsed_meta, axis='observation')
	pa_sorted_biom_collapsed = pa_sorted_biom_collapsed.sort_order(collapsed_order, axis='observation')

#	print pa_sorted_biom_collapsed.to_tsv(header_key='taxonomy', header_value='taxomomy')
#	print "\n\n"

	return pa_sorted_biom_collapsed

def add_sample_metadata_to_biom(in_table, metadata, v=0):

    	from biom.table import Table
	
	sample_metadata = {}
	for s in in_table.ids(axis='sample'):
		if v:
        		print "adding metadata to: %s" %s
        	sample_metadata[s] = {}
		for meta in metadata[s].keys():
			sample_metadata[s][meta] = metadata[s][meta]

	#add metadata to individual table
	in_table.add_metadata(sample_metadata, axis='sample')


def assign_taxonomy_kraken(kraken_out, tax_dict, v=0):
    """
    finds taxonomic level of assignemt based on taxid and bins into dictionary
    """
    from collections import defaultdict

    levels = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'] #taxonomic levels of interest
    tax_count = {'kingdom':{}, 'phylum':{}, 'class':{}, 'order':{}, 'family':{}, 'genus':{}, 'species':{}}
    minimum = tax_dict["tax_id"].index("kingdom")

    print "\nInterpreting kraken results and adjust to standard taxonomic levels (%s)\n" %(", ".join(levels))

    if not tax_count.has_key('nohit'):
        tax_count['nohit'] = {'nohit':[]}

    if kraken_out.has_key('nohit'):
        tax_count['nohit']['nohit'].extend(kraken_out['nohit'])

    if kraken_out.has_key('hit'):

        for query in kraken_out['hit']:

            #find the index of the initial assignment
            index = tax_dict["tax_id"].index(tax_dict[kraken_out['hit'][query][0]][1])

            #if already the initial assignemnt is above kingdom level
            if index < minimum:
                if v:
                    print "initial assignment above level kingdom for %s -> bin as 'nohit'" %query
                tax_count['nohit']['nohit'].append(query)
                continue

            #if the initial assignment hits one of the focal taxonomic levels
            if tax_dict[kraken_out['hit'][query][0]][1] in levels and tax_dict[kraken_out['hit'][query][0]][2]:
                if v:
                    print "direct assignment for %s -> %s (%s)" %(query, kraken_out['hit'][query][0], tax_dict[kraken_out['hit'][query][0]][2])
                if not tax_count[tax_dict[kraken_out['hit'][query][0]][1]].has_key(kraken_out['hit'][query][0]):
                    tax_count[tax_dict[kraken_out['hit'][query][0]][1]][kraken_out['hit'][query][0]] = []
                tax_count[tax_dict[kraken_out['hit'][query][0]][1]][kraken_out['hit'][query][0]].append(query)

            #if the initial assignment does not hit one of the focal taxonomic levels
            else:
                if v:
                        print "invalid assignment level for %s\t%s - %s - '%s'" %(query, tax_dict[kraken_out['hit'][query][0]][1], kraken_out['hit'][query][0], tax_dict[kraken_out['hit'][query][0]][2])
                        print "\tsearching through higher taxonomic levels:"
                ok = 0

                #we'll count down, i.e. move taxonomic levels up one by one until we hit a level that is acceptable
                while index >= minimum and not ok:
                    index-=1
                    if tax_dict["tax_id"][index] in levels: #check if the current taxonomic level is a valid one
                                                    
                        if tax_dict[kraken_out['hit'][query][0]][index]: #check if there is actually a taxid availble for the lineage at this level
                            if v:
                                print "\tlevel '%s' is among the targets - taxid present -> OK!" %tax_dict["tax_id"][index]
                            ok = 1
                        else: #
                            if v:
                                print "\tlevel '%s' is among the targets, but no taxid available at this level -> moving on" %tax_dict["tax_id"][index]
                        
                            
#                    else:
#                        print "\tlevel '%s' is not acceptable" %tax_dict["tax_id"][index]
                        
                if ok:
                    taxid = tax_dict[kraken_out['hit'][query][0]][index]

                    if v:
                        print "\tadjusted assignment for %s -> %s (%s)" %(query, taxid, tax_dict[taxid][2])
                    if not tax_count[tax_dict["tax_id"][index]].has_key(taxid):
                        tax_count[tax_dict["tax_id"][index]][taxid] = []
                    tax_count[tax_dict["tax_id"][index]][taxid].append(query)
                else:
                    if v:
                        print "\tCouldn't find a valid assignment for '%s'" %query
                    tax_count['nohit']['nohit'].append(query)

    #cleanup - remove any taxonomic levels that did not get assignments from the dictionary
    for key in tax_count.keys():
        if not tax_count[key]:
            del(tax_count[key])

    return tax_count


def find_full_taxonomy(per_tax_level, taxonomy_dictionary):
        
        """
        Extracts taxonomy strings consisting of definend taxonomic levels for taxids
        """
        
        syn = {'kingdom': 'k__', 'phylum': 'p__', 'class': 'c__', 'order': 'o__', 'family': 'f__', 'genus':'g__', 'species': 's__'}
        levels = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
        level_indices = []
        tax_trees = {}

        for level in levels:
            tax_trees[level] = {}
            for i in range(len(taxonomy_dictionary['tax_id'])):
                if level == taxonomy_dictionary['tax_id'][i]:
                    level_indices.append(i)

        for level in reversed(levels):
            if level in per_tax_level:
                
                for tid in per_tax_level[level]:
#                    print "\n\ttaxid: %s" %tid
                    ind_taxonomy = []
                    for i in range(len(level_indices)):#range(len(levels)):
#                        print "seeking level %s at index: %s" %(taxonomy_dictionary['tax_id'][level_indices[i]], level_indices[i])
                        if taxonomy_dictionary[tid][level_indices[i]]:
                            ind_taxonomy.append('%s%s' %(syn[levels[i]],taxonomy_dictionary[taxonomy_dictionary[tid][level_indices[i]]][2].replace(' ','_')))
                        else:
                            ind_taxonomy.append('%sunknown' %syn[levels[i]])

                    for j in reversed(range(len(levels))):
                        if '_unknown' in ind_taxonomy[j]:
                            del ind_taxonomy[j]
                        else:
                            break

                    tax_trees[level][tid] = ind_taxonomy

        if per_tax_level['nohit']:
            tax_trees['nohit'] = {'unassigned': ['__unassigned']}

        for rank in tax_trees.keys():
            if not tax_trees[rank]:
                del tax_trees[rank]

        return tax_trees


def global_uc_to_biom(clust_dict, query_dict):
	
	import numpy as np	
	from biom.table import Table

	data_to_biom = []
	observation_ids = []
	sample_ids = []
	sample_metadata = []
	observation_metadata = []

	observation_ids = clust_dict.keys()
	sample_ids = sorted(query_dict.keys())
	
	for s_id in sample_ids:
		sample_metadata.append({})
#		sample_metadata[s_id] = {'metadata':{}}

#	print sample_ids
#	print observation_ids	

	for OTU in observation_ids:
#		print "OTU: %s" %OTU
#		print clust_dict[OTU]
		per_OTU = []
		for sample in sample_ids:
#			print "\t%s" %sample
			per_OTU.append(int(0))
			for i in range(len(clust_dict[OTU])):
				if clust_dict[OTU][i].startswith(sample+"|"):
					per_OTU[-1] += int(query_dict[sample]['cluster_counts'][clust_dict[OTU][i].split("|")[1]])
#					break

		data_to_biom.append(per_OTU)
		
	data = np.asarray(data_to_biom)

#	print "\n### FINAL CHECK: ###\n"
#	print "number of OTUs: %s" %len(observation_ids)
#	print "number of samples: %s" %len(sample_ids)
#	print "number of observations: %s" %len(data)
	
	table = Table(data, observation_ids, sample_ids, table_id='OTU table', sample_metadata=sample_metadata)
#	print table
	return table


def concatenate_for_global_clustering(queries_dict, out):
	"""
	The function concatenates the query files from all samples and adds the sample IDs to the headers in the process
	"""
	from Bio import SeqIO
	OUT = open(out, 'w')
	for sampleID in queries_dict.keys():
		
		queries = list(SeqIO.parse(open(queries_dict[sampleID]['queryfile']),'fasta'))
		
		for rec in queries:
			rec.description = sampleID+'|'+rec.id
			rec.id = rec.description
	
		SeqIO.write(queries, OUT, 'fasta')
		
	OUT.close()

def vsearch_cluster(infile, cluster_match, threads, sampleID):
	"""
	The function runs vsearch
	"""
	import shlex, subprocess

	cmd = "vsearch --cluster_fast %s --id %.2f --strand both --threads %s --centroids %s_centroids.fasta --uc %s.uc" % (infile, cluster_match, threads, sampleID, sampleID)
        print cmd
        cmdlist = shlex.split(cmd)
        stdout, stderr = subprocess.Popen(cmdlist, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate() # , stdout=subprocess.PIPE).communicate()
        if stdout:
	        print stdout


def concat_and_filter_by_length(inlist, outfile, excludefile, form='fasta', length=0, devi=0.1):
	"""
	The Function concatenates sequences from a list of files into a single file and
	 applies a length filter if specified
	"""

	import gzip
	from Bio import SeqIO


	exclude = []
	count = 0
	OUT = open(outfile,'w')
	
	for f in inlist:
		if f.endswith('.gz'):
			IN = gzip.open(f,'rb')
		else:
			IN = open(f, 'r')
    		
		for record in SeqIO.parse(IN,'fastq'):
			if length:
                		if len(record.seq) < length*(1-devi) or len(record.seq) > length*(1+devi):
					record.id += '|length_filtered'
					record.description = record.id
					exclude.append(record)
					continue
			record.id = record.description.replace(" ","_")
			record.description = record.id
			SeqIO.write(record, OUT, form)
			count += 1
	OUT.close()

	if length:
		print "\n%i sequences removed by length filter\n" %len(exclude)

	if excludefile:
		EXCLUDE = open(excludefile,'w')
		SeqIO.write(exclude, EXCLUDE, 'fasta')
		EXCLUDE.close()
	del exclude

	return count


def keep_only(inlist, pattern_list):
	"""
	The function reduced a list to only the elements that contain a pattern
	"""
	for i in reversed(range(len(inlist))):
		for j in range(len(pattern_list)):
			if pattern_list[j] in inlist[i]:
				break
			
			if j == len(pattern_list)-1:
				del inlist[i]


def check_email(mail):
        """
        The function checks that you provide an email address to Entrez
        """

        print "\nmetaBEAT may be querying NCBI's Entrez databases to fetch/verify taxonomic ids. Entrez User requirements state that you need to identify yourself by providing an email address so that NCBI can contact you in case there is a problem.\n"


        if not mail:
                print "As the mail address is not specified in the script itself (variable 'Entrez.email'), metaBEAT expects a simple text file called 'user_email.txt' that contains your email address (first line of file) in the same location as the metaBEAT.py script (in your case: %s/)\n" %os.path.dirname(sys.argv[0])
                if not os.path.isfile(os.path.dirname(sys.argv[0])+'/user_email.txt'):
                        print "Did not find the file %s/user_email.txt - you may specify your email address also via the '-@' command line option\n" %os.path.dirname(sys.argv[0])
                        sys.exit()
                now = datetime.today()
                modify_date = datetime.fromtimestamp(os.path.getmtime(os.path.dirname(sys.argv[0])+'/user_email.txt'))
                if (now-modify_date).days > 7:
                        print "%s/user_email.txt is older than 7 days - Please make sure it's up to date or specify email address via the '-@' option\n" %os.path.dirname(sys.argv[0])
                        sys.exit()
                FH = open(os.path.dirname(sys.argv[0])+'/user_email.txt','r')
                mail = FH.readline().strip()
                FH.close()
                if not '@' in mail:
                        print "\nnot sure %s is an email address - please make sure it's valid\n" %mail
                        sys.exit()
                print "found '%s' in %s/user_email.txt\n" %(mail, os.path.dirname(sys.argv[0]))

        else:
                if not '@' in mail:
                        print "\nnot sure %s is an email address - please make sure it's valid\n" %mail
                        sys.exit()
                print "You have specified: '%s'\n" %(mail)

                FH = open(os.path.dirname(sys.argv[0])+'/user_email.txt','w')
                FH.write(mail)
                FH.close()

        return mail

	
def rw_gb_to_taxid_dict(dictionary, name, mode):
    '''
    The function writes a dictionary to a file ('key,value')
    or reads a comma separated file into a dictionary
    '''
 
    if mode == 'w':
        fh = open(name, 'w')
        for key in dictionary.keys():
            fh.write("%s,%s\n" %(key, dictionary[key]))
        fh.close()
	print " ..done writing!\n"
    
    
    elif mode == 'r':
        fh = open(name, 'r')
        for line in [l.strip() for l in fh]:
            key,value = line.split(',')
            dictionary[key] = value
        fh.close()
	print " ..done parsing!\n"
        
    else:
        sys.exit("only 'r' or 'w' are allowed for mode in the rw_gb_to_taxid_dict")
        

def filter_centroid_fasta(centroid_fasta, m_cluster_size, cluster_counts, sampleID, v=False):
    "This function filters the centroid fasta file produced by vsearch"
    "and removes all centroid sequences that represent clusters of size"
    "below the specified threshold"

    from Bio import SeqIO

    badstring="" #start with an emtpy string
    good_seqs = []
    bad_seqs = []
    badcount = 0
    os.rename(centroid_fasta, centroid_fasta+"_backup")
    seqs=list(SeqIO.parse(centroid_fasta+"_backup",'fasta'))
    total_count = sum(cluster_counts.values())
    for record in seqs:
        if int(cluster_counts[record.id]) >= m_cluster_size:
	    good_seqs.append(record)
	    v_string = "%s - ok" %record.id
        else:
	    v_string = "excluding %s from further analyses - cluster_size filter\n" %record.id
	    bad_seqs.append(record)
	    badcount += 1
	if v:
	    print v_string
    
    if bad_seqs:
	bad=open(centroid_fasta+"_removed", "w")
	for record in bad_seqs:
#            badstring+=">%s|%s|%s|%.2f|minimum_cluster_filter\n%s\n" % (sampleID, record.id, cluster_counts[record.id], float(cluster_counts[record.id])/(total_count-badcount)*100, record.seq) #potetial problems with divisions by zero
            badstring+=">%s|%s|%s|minimum_cluster_filter\n%s\n" % (sampleID, record.id, cluster_counts[record.id], record.seq) #potetial problems with divisions by zero
            del cluster_counts[record.id]

	bad.write(badstring)
	bad.close()
	del bad_seqs

    total_count = sum(cluster_counts.values())
   
    f=open(centroid_fasta,"w")
    SeqIO.write(good_seqs, f, "fasta")
#    for record in good_seqs:
#        outstring=">%s|%s|%s|%.2f\n%s\n" % (sampleID, record.id, cluster_counts[record.id], float(cluster_counts[record.id])/total_count*100, record.seq)
#        f.write(outstring)
    f.close()
    del good_seqs

    return total_count



def parse_vsearch_uc(fil, cluster_counts, cluster_reads, extract_reads=0):
    "The function parses the 'uc' file from vsearch"
    f=open(fil,"r")

    read_count=0

    for line in [l.strip() for l in f]: #loop through the file one line at a time, stripping off any newline characters
        if line.startswith("C"): #process only lines that start with a "C"
            elem = line.split("\t") #split the lines at tab
            cluster_counts[elem[8]] = int(elem[2]) #write the counts to dictionary with key being the id of the centroid read

	if line.startswith('H'):
		read_count+=1
        
        if extract_reads:
            if line.startswith('S'):
                elem = line.split("\t")
                if not cluster_reads.has_key(elem[8]):
                    cluster_reads[elem[8]] = [elem[8]] #create a new key for the centroid id and add the centroid id as the first element into the list
                    
            if line.startswith('H'):
                elem = line.split("\t")
                if cluster_reads.has_key(elem[9]):
                    cluster_reads[elem[9]].append(elem[8]) #add the new read id to the centroid cluster
                else:
                    cluster_reads[elem[9]] = [elem[9]] #create a new key for the centroid id and add the centroid id as the first element into the list
                    cluster_reads[elem[9]].append(elem[8]) #add the new read id to the centroid cluster
    f.close()

    read_count += len(cluster_counts)
	
    return read_count


def assign_taxonomy_LCA(b_filtered, tax_dict, v=0):
    "The function takes a dictionary of queries and their hits"
    "and provides taxonomic assignment based on the LCA method."
    tax_count = defaultdict(dict)
    levels = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'] #taxonomic levels of interest

#    print len(b_filtered['hit'])
    if b_filtered.has_key('hit'):
            minimum = tax_dict["tax_id"].index("kingdom") #find the minimum index, i.e. the index for 'kingdom'
	    for query in b_filtered['hit'].keys():
	        if len(b_filtered['hit'][query]) == 1:
	            if v:
		            print "\ndirect assignment for %s -> %s" %(query, b_filtered['hit'][query][0])

	            if not tax_count['species'].has_key(b_filtered['hit'][query][0]):
	                tax_count['species'][b_filtered['hit'][query][0]] = []
	            tax_count['species'][b_filtered['hit'][query][0]].append(query)
	            del b_filtered['hit'][query]
	        else:
		    if v:
		            print "\nattempting LCA assignment for %s" %query
	            for index in reversed(range(len(tax_dict["tax_id"]))):
	                id_list = []
			if tax_dict["tax_id"][index] == 'below_superkingdom': #if this level is reached (that means it cannot be assigned to a kingdom) in the taxonomy then we have tested all possible levels and have not found any convergencce
			    if v:
				    print "was unable to assign LCA to %s" %query
			    if not tax_count.has_key('nohit'):
	            		tax_count['nohit'] = {'nohit':[]}
	        	    tax_count['nohit']['nohit'].append(query)
			    break
#	                print index
#	                print "\nLEVEL: %s" %tax_dict["tax_id"][index]
	                for tax in b_filtered['hit'][query]:
#			    print tax
#	                    print tax_dict[tax][index]
	                    id_list.append(tax_dict[tax][index])
	                    if not tax_dict[tax][index]:
#	                        print "nothing found at this level for %s" %tax
	                        break
	                if len(id_list) == len(b_filtered['hit'][query]):
#	                        print "ok - all have valid taxid"

	                        if len(set(id_list)) == 1:
                                    LCA_taxid=""
                                    ok=0
                                    if v:
                                        print "found LCA %s at level %s" %(id_list[0], tax_dict["tax_id"][index])
                                    while not ok and index >= minimum:
#                                       print "checking level %s (index: %i)" %(tax_dict["tax_id"][index], index)
                                        if tax_dict["tax_id"][index] in levels:
#                                           print "level ok\n"
                                            if tax_dict[id_list[0]][index]:
                                                LCA_taxid = tax_dict[id_list[0]][index]
                                                ok=1
                                                if v:
                                                        print "assigned LCA %s (taxid %s) at level %s" %(tax_dict[LCA_taxid][2], LCA_taxid, tax_dict["tax_id"][index])
                                            else:
                                                if v:
                                                        print "no valid taxid found at level %s" %tax_dict["tax_id"][index]
                                                index-=1
                                        else:
                                            if v:
                                                print "assignment to level %s (taxid %s) is not acceptable" %(tax_dict["tax_id"][index], tax_dict[id_list[0]][index])
                                            index-=1
##                                  print tax_dict[id_list[0]][1]
##                                  print tax_dict[id_list[0]][2]
                                    if LCA_taxid:
                                            if not tax_count[tax_dict["tax_id"][index]].has_key(LCA_taxid):
                                                tax_count[tax_dict["tax_id"][index]][LCA_taxid] = []
                                            tax_count[tax_dict["tax_id"][index]][LCA_taxid].append(query)
#                                           print "%s\t%s" %(tax_dict["tax_id"][index], tax_count[tax_dict["tax_id"][index]][LCA_taxid])
                                            del b_filtered['hit'][query]
                                            break
                                    else:
                                        if v:
                                                print "was unable to assign LCA to %s" %query
                                        if not tax_count.has_key('nohit'):
                                                tax_count['nohit'] = {'nohit':[]}
                                        tax_count['nohit']['nohit'].append(query)
                                        break

#	                        else:
#	                            print "not yet LCA"
                            
#	        print "\nUPDATE:\n%s" %tax_count
    	    if len(b_filtered['hit']) == 0:
	        print "\nall queries have been successfully assigned to a taxonomy"
	    else:
                print "\nLCA detection failed for %i queries\n" %(len(b_filtered['hit']))
                if v:
                        for q in b_filtered['hit']:
                                print q,b_filtered['hit'][q]


    if b_filtered.has_key('nohit'):
        if not tax_count.has_key('nohit'):
            tax_count['nohit'] = {'nohit':[]}
        tax_count['nohit']['nohit'].extend(b_filtered['nohit'])
    
        
    return tax_count

def write_taxids(tids, out='taxids.txt'):
    print "write out taxids to taxids.txt\n"
    f = open(out,"w")
    for tid in tids:
        f.write(tid + "\n")
    f.close()
	

def make_tax_dict(tids, out_tax_dict, denovo_taxa, ref_taxa):
    "The function takes a list of taxonomic ids, runs the taxtable program from the taxit suite"
    "to summarize the taxonomic information for the taxids and formats the result into a dictionary"

    import shlex, subprocess
    import re

    write_taxids(tids)

    cmd = "taxit taxtable -f taxids.txt -o taxa.csv %s" %taxonomy_db
    print "running taxit (v0.8.5) to generate reduced taxonomy table"
    if len(denovo_taxa) > 0:
        print "WARNING: any taxa without valid taxid will not be included -  NEEDS TO BE FIXED IF PHYLOGENETIC PLACEMENT IS PLANNED"
    print "\n"+cmd
    taxtable,err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    if err:
	print "something went wrong while building the reduced taxonomy:"
	print err
	sys.exit()

    taxtable = open("taxa.csv","r")
    for line in taxtable:
        line = re.sub('"','',line.strip())
        array = line.split(',')

        if len(array) > 1:
            key = array.pop(0)
#            if array[-1] == 'subspecies':
#                array.pop(-1)
            out_tax_dict[key]=array
            
            #add denovo taxa dummies to tax_dict
            if denovo_taxa.has_key(key):
                for deno in denovo_taxa[key]:
                    local = array[:]
                    local[0] = key
                    local[1] = 'species'
                    local[2] = deno.split('|')[0]
                    local[-1] = deno.split('|')[1] #reference_taxa[deno]
                    out_tax_dict[deno.split('|')[1]] = local

#	while out_tax_dict['tax_id'][-1] != 'species': #this should get rid of any 'subspecies' or 'varietas' levels
#            out_tax_dict['tax_id'].pop(-1)

def update_taxonomy(v):

    import shlex, subprocess

    filename = os.path.dirname(sys.argv[0])+'/taxonomy.db'
    print "\nUpdate/compile taxonomy database at %s - might take a minute or two\n" %filename

    if os.path.isfile(filename):
	os.remove(filename)
    
    cmd="taxit new_database --download-dir %s/ %s" %(os.path.dirname(filename), filename)
    print cmd
    cmdlist = shlex.split(cmd)
    proc = subprocess.Popen(cmdlist, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    exitcode = proc.returncode
#    print "EXITCODE: %s" %exitcode

#    stdout, stderr = subprocess.Popen(cmdlist, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate() # , stdout=subprocess.PIPE).communicate()
    if v:
        if stdout:
            print stdout
        if stderr:
            print stderr
    if exitcode:
        print stderr
        sys.exit('\nDownloading/compiling the taxonomy database failed - see error message above\n')

    path = os.path.dirname(filename)
    os.remove(os.path.dirname(filename)+'/taxdmp.zip')
    print "Done!\n"
    sys.exit()

def check_for_taxonomy_db(tax_db):
    "The function checks if the taxonomy.db is present"
    if not tax_db or not os.path.isfile(tax_db): #if no path to the taxonomy database is specifed or the path is not correct
                if os.path.isfile(os.path.dirname(sys.argv[0])+'/taxonomy.db'): #check if taxonomy.db is present in the same path as the metabeat.py script
                        tax_db = "%s/taxonomy.db" %os.path.dirname(sys.argv[0])
			print "taxonomy.db found at %s" %tax_db
			
                	now = datetime.today()
                	modify_date = datetime.fromtimestamp(os.path.getmtime(os.path.dirname(sys.argv[0])))
                	if (now-modify_date).days > 30:
                        	print "Your taxonomy database is older than 30 days - might not be a problem, but you can update it via 'metaBEAT_global.py --update_taxonomy' - good luck!\n"
			else:
				print "Taxonomy database is %s days old - good!" %(now-modify_date).days
                else:
                        print "\nIf you are attempting taxonomic assignmnt, metaBEAT.py requires a taxonomy database configured by taxtastic. Per default metaBEAT expects a file 'taxonomy.db' in the same directory that contains the metaBEAT.py script. If you have the database in a different location you can specify this via the '--taxonomy_db' option. If you don't have the database you can run 'metaBEAT_global.py --update_taxonomy', which will download and compile the database for you.\n"
                        sys.exit()
    else:
        print "taxonomy.db is present at: %s" %tax_db

    return tax_db

def gb_to_taxid(b_filtered, all_taxids, processed, gb_to_taxid_file=args.gb_to_taxid, v=0):
    "This function takes the resulting dictionary from the blast_filter function"
    "and fetches the corresponding taxonomic ids to the gb accessions if necessary"
    "(i.e. if the format is specified as 'gb')."
    "It returns a list of unique taxonomic ids."
    if b_filtered['format'] == 'gb':
	total = 0
        for query in b_filtered['hit'].keys():
	    new = 0
	    old = 0
            taxids=[]
            if v:
                print "processing query: %s" %query
            for hit in b_filtered['hit'][query]:
#                print "gi: %s" %hit
                if processed.has_key(hit):
                    taxid = processed[hit]
                    old += 1
                    if v:
                        print "have seen gb accession '%s' before" %(hit)
                else:
                    
                    if v:
                        print "querying Genbank to fetch taxid for gb accession '%s'" %(hit),
		    done = False
	            while not done:
			try:
				handle = Entrez.efetch(db="nucleotide", id=hit, rettype="fasta", retmode="xml")
				done = True
			except:
				print "connection closed - retrying entrez for query '%s'.." %hit

                    taxon = Entrez.read(handle)
                    taxid = taxon[0]['TSeq_taxid']
                    processed[hit] = taxon[0]['TSeq_taxid']
                    time.sleep(0.25)
                    new += 1
		    if v:
			print " .. success!"
                        
#                print "taxid: %s" %taxid                
                taxids.append(taxid)
		total += 1
                if not v:
                    if total%100 == 0:
                        print "accessions processed:\t%i" %total
                
#            print b_filtered['hit'][query]
            b_filtered['hit'][query] = list(set(taxids))
#            print len(taxids)
            all_taxids.extend(taxids)
#            print len(all_taxids)

	    if new > 0:
		if v:
		    print "detected %i new accessions -> updating local gb_to_taxid record" %new
		rw_gb_to_taxid_dict(dictionary=processed, name=gb_to_taxid_file, mode='w')
    
        if v:
            print "\nDone! - gathered taxids for %i accessions" %total
        
    elif b_filtered['format'] == 'taxid':
        for query in b_filtered['hit'].keys():
            all_taxids.extend(b_filtered['hit'][query])

#    print all_taxids
#    print "make to set"
    all_taxids = list(set(all_taxids))
    return all_taxids


def pre_pplacer_filter(blast_handle, m_ident=0.8, m_ali_length=0.95, v=0):
	"""
	Parse blast xml and bin queries with <80% similarity. 
	"""
	
	result = {'format':'taxid', 'nohit':[], 'to_pplacer':[], 'strand_reversed':[]}
	count = 0
	for res in blast_handle:
		count+=1
		if not res.alignments:
			result['nohit'].append(res.query)
		elif (len(res.alignments[0].hsps[0].query) < int(res.query_length*m_ali_length)): #check minimum alignment length
			result['nohit'].append(res.query)
		elif (float(res.alignments[0].hsps[0].identities)/len(res.alignments[0].hsps[0].query) < m_ident): #check if hit is under the identity threshold
			result['nohit'].append(res.query)
		else: #this is a hit worth sending to pplacer
			result['to_pplacer'].append(res.query)
			if res.alignments[0].hsps[0].sbjct_start > res.alignments[0].hsps[0].sbjct_end:
                		result['strand_reversed'].append(res.query)

	print " number of queries processed: %i" %count
	print "\nNumber of queries passed to pplacer:\t%i\n(%i will be reverse complemented first)" %(len(result['to_pplacer']), len(result['strand_reversed']))
	print "Number of queries directly binned to 'unassigned':\t%i" %len(result['nohit'])

	return result		

def makeblastdb(in_fasta, dbtype, out_prefix, v=0):

	import os
	import shlex, subprocess

	cmd="makeblastdb -in %s -dbtype %s -out %s_blast_db" %(in_fasta, dbtype, out_prefix)
        print cmd
        cmdlist = shlex.split(cmd)
        cmd = subprocess.Popen(cmdlist, stdout=subprocess.PIPE)
        stdout = cmd.communicate()[0]
        if v:
		print stdout
        
	return os.path.abspath('.')+"/"+"%s_blast_db" %out_prefix

		
def blast_filter(b_result, v=0, m_bitscore=80, m_ident=0.8, m_ali_length=0.95, bit_score_cutoff=0.1, bit_score_cutoff_adjust_off=False, strand_reversed=None):
    """
    The function interprets a BLAST results handle and filters results for subsequent taxonomic assignment using LCA
    """
   
    #Summarize filter parameters:
    print "\tFiltering parameters:"
    print "\tMinimum bitscore: %s" %m_bitscore
    print "\tMinimum alignment length: %s" %m_ali_length
    print "\tMinimum identity across alignment length: %s" %m_ident
    print "\tBitscore skim window for LCA: %s" %bit_score_cutoff
    print "\tBitscore skim window adjustment (for 100% identity matches): ",
    if bit_score_cutoff_adjust_off:
	print "off\n"
    else:
	print "on\n"
    ###############  
 
    result = {'format':''}
    count=0
    for res in b_result:
        bit_score_cutoff_set = 1-bit_score_cutoff #specify the bit score cutoff (n
        count += 1
        if v:
            print "\nquery: %s" % res.query #the current query
            print "query_length: %s" %str(res.query_length)
        if not res.alignments:  #if no alignment was found for the query
            if v:
                print "no hit - done"
            if not result.has_key('nohit'):
                result['nohit'] = []
            result['nohit'].append(res.query)

        elif (len(res.alignments[0].hsps[0].query) < int(res.query_length*m_ali_length)):
                if v:
                        print "alignment length (%i) below threshold (%i * %s -> %i)" %(len(res.alignments[0].hsps[0].query), res.query_length, m_ali_length, (res.query_length*m_ali_length))      
                if not result.has_key('nohit'):
                        result['nohit'] = []
                result['nohit'].append(res.query)
        elif (float(res.alignments[0].hsps[0].identities)/len(res.alignments[0].hsps[0].query) < m_ident) or res.alignments[0].hsps[0].bits < m_bitscore:
            if v:
                print "no significant hit - done"
            if not result.has_key('nohit'):
                result['nohit'] = []
            result['nohit'].append(res.query)

        else: #if a reasonable hit was found
	    
#	    print res.query,res.alignments[0].hsps[0].query_start,res.alignments[0].hsps[0].query_end,res.alignments[0].hsps[0].sbjct_start,res.alignments[0].hsps[0].sbjct_end
	    if res.alignments[0].hsps[0].sbjct_start > res.alignments[0].hsps[0].sbjct_end:
		if v:
			print "query inverted"
		strand_reversed.append(res.query)

            if not result.has_key('hit'):
                result['hit'] = {}

            if (float(res.alignments[0].hsps[0].identities)/len(res.alignments[0].hsps[0].query) == 1): #if we have a 100 % match of at least the minimum length, adjust the bitscore cutoff so that only hits with the same bitscore as the top hit are considered for the LCA
		message = "100%% match - Alignment length: %s" %str(len(res.alignments[0].hsps[0].query))
		if not bit_score_cutoff_adjust_off:
                        bit_score_cutoff_set = 1
                        if v:
                                message+=" -> adjusting bitscore cutoff"
		if v:
                        print message
                
            max_bit_score = res.alignments[0].hsps[0].bits #record the maximum bitscore

            result['hit'][res.query]=[] #create empty list for query
            for alignment in res.alignments: #for each hit of the current query
#                print alignment.hsps[0].bits
                if not result['format']: #determining the format of the blast database (only once per blast file)
                    if alignment.title.startswith('gi') or alignment.title.startswith('gb'):
                        result['format']='gb'

			found_index = False
			for tag in ['gb','dbj']:
				if tag in alignment.title.split("|"):
					gb_index = alignment.title.split("|").index(tag)+1
					found_index = True

			if not found_index:
				sys.exit("\nCan't recognize the header format in the BLAST database: %s\n" %alignment.title.split(" ")[0])
                    else:
                        result['format']='taxid'

                if alignment.hsps[0].bits >= (max_bit_score*bit_score_cutoff_set): #if a hit has a bitscore that falls within the specified window it will be collected for later LCA
#                    print alignment.title.split("|")[1]
                    if result['format'] == 'gb':
                        result['hit'][res.query].append(alignment.title.split("|")[gb_index])
                    elif result['format'] == 'taxid':
                        result['hit'][res.query].append(alignment.title.split("|")[-2])

                    if v:
                        print "%s\t%s" %(alignment.hsps[0].bits, alignment.title)
                else:
#                    print result['hit'][res.query]
                    break
#            print "%s: %s" %(res.query, result['hit'][res.query])

    print "%i queries processed" %count

    if not result['format']:    #if no database format could be determined at this point that means there was no significant hit, so we set the format to unknown
        result['format'] = 'unknown'

    return result

def file_check(file_to_test, optional_message=None):
	"tests if a file exists"
	file_to_test = os.path.abspath(file_to_test)
	if not os.path.isfile(file_to_test):
		exit_message = "%s is not a valid file" %file_to_test
		if optional_message:
			exit_message = optional_message
#		parser.print_help()
		sys.exit(exit_message)
	return True

def write_out_refs_to_fasta(ref_seqs, ref_taxids = {}):
	print "write out reference sequences to refs.fasta\n"
	OUT = open("refs.fasta","w")
#	OUT_temp = open ("temp_refs.fasta","w")
	for record in ref_seqs:
#		outstring = ">%s\n%s" % (record.features[0].qualifiers['db_xref'][0].split(":")[1], record.seq) #note that the header of the sequence will be the taxid
#		if record.features[0].qualifiers.has_key('db_xref'):
#			outstring = ">%s|%s|%s\n%s" % (record.id, record.features[0].qualifiers['db_xref'][0].split(":")[1], record.features[0].qualifiers['organism'][0], record.seq) #note that the header of the sequence will be the taxid
#		else:
#			outstring = ">%s|%s|%s\n%s" % (record.name, ref_taxids[record.features[0].qualifiers['organism'][0]], record.features[0].qualifiers['organism'][0], record.seq)
		
		for db_xref in record.features[0].qualifiers['db_xref']:
			if 'taxon:' in db_xref:
				loc_taxid = db_xref.split(":")[1]
#		outstring = ">%s|%s|%s\n%s" % (record.name, ref_taxids[record.features[0].qualifiers['organism'][0]], record.features[0].qualifiers['organism'][0], record.seq)
		outstring = ">%s|%s|%s\n%s" % (record.name, loc_taxid, record.features[0].qualifiers['organism'][0], record.seq)

#		outstring2 = ">%s\n%s" % (record.id, record.seq)
		OUT.write(outstring + "\n")
#		OUT_temp.write(outstring2 + "\n")

	OUT.close()
#	OUT_temp.close()

def annotation_BIOM_table_with_taxonmy(BIOM_table, hits_per_tax_level, taxonomy_dictionary, method, version=VERSION, prefix=args.output_prefix):
	"""
	Annoates a denovo BIOM table with a taxonomy and writes the new tables out
	"""
	from biom.table import Table

	print "\n\n##### FORMATTING AND WRITING BIOM OUTPUT #####\n"

	BIOM_tables_loc = {}
	taxonomic_trees = {}
	taxonomic_trees = find_full_taxonomy(per_tax_level=hits_per_tax_level, taxonomy_dictionary=taxonomy_dictionary)

	BIOM_tables_loc['OTU_taxonomy'] = add_taxonomy_to_biom(per_tax_level_clusters=hits_per_tax_level, per_tax_level_trees=taxonomic_trees, biom_in=BIOM_table, method=method)
	BIOM_tables_loc['taxon_taxonomy'] = collapse_biom_by_taxonomy(in_table=BIOM_tables_loc['OTU_taxonomy'])
	BIOM_tables_loc['cluster_taxonomy'] = pa_and_collapse_by_taxonomy(in_table=BIOM_tables_loc['OTU_taxonomy'])

#	print tables
	out=open(prefix+"-OTU-taxonomy."+method+".biom","w")
	BIOM_tables_loc['OTU_taxonomy'].to_json('metaBEAT v.'+version, direct_io=out)
	out.close()

	out=open(prefix+"-OTU-taxonomy."+method+".tsv","w")
	out.write(BIOM_tables_loc['OTU_taxonomy'].to_tsv(header_key='taxonomy', header_value='taxomomy')) #to_json('generaged by test', direct_io=out)
	out.close()

	out=open(prefix+"-by-taxonomy-readcounts."+method+".biom","w")
	BIOM_tables_loc['taxon_taxonomy'].to_json('metaBEAT v.'+version, direct_io=out)
	out.close()

	out=open(prefix+"-by-taxonomy-readcounts."+method+".tsv","w")
	out.write(BIOM_tables_loc['taxon_taxonomy'].to_tsv(header_key='taxonomy', header_value='taxomomy')) #to_json('generaged by test', direct_io=out)
	out.close()

	out=open(prefix+"-by-taxonomy-clustercounts."+method+".biom","w")
	BIOM_tables_loc['cluster_taxonomy'].to_json('metaBEAT v.'+version, direct_io=out)
	out.close()

	out=open(prefix+"-by-taxonomy-clustercounts."+method+".tsv","w")
	out.write(BIOM_tables_loc['cluster_taxonomy'].to_tsv(header_key='taxonomy', header_value='taxomomy')) #to_json('generaged by test', direct_io=out)
	out.close()

	return BIOM_tables_loc

def gb_to_kraken_db(Sequences, fasta_for_kraken):
    	"""
    	convert Genbank data to fasta file ready to be formatted as kraken database
    	"""
    
	from Bio import SeqIO
	Seqs_new = []

	for r in Sequences:
		source=r.features[0]
		for t in source.qualifiers['db_xref']:
			if 'taxon' in t:
				taxid = t.split(":")[1]
		r.id = "%s|kraken:taxid|%s" %(r.id,taxid)
    		r.description = r.id
    		Seqs_new.append(r)
    
	out = open(fasta_for_kraken+'.kraken.fasta','w')
	SeqIO.write(Seqs_new, out, 'fasta')
	out.close()

def initiate_custom_kraken_db(db_path='./', db_name='custom', v=0):
    """
    Initiates a custom kraken database 
    """
    
    import shlex, subprocess
    
    print "\n## Initiate kraken database and download taxonomy ##\n"
    
    cmd="kraken-build --download-taxonomy --db %s/%s\n" %(db_path,db_name)
    
    print cmd
    cmdlist = shlex.split(cmd)
    proc = subprocess.Popen(cmdlist, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    exitcode = proc.returncode #exitcode is '0' if process finished without error

#    stdout, stderr = subprocess.Popen(cmdlist, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate() # , stdout=subprocess.PIPE).communicate()
#    exitcode = 
    if v:
	if stdout:
            print stdout
#	if stderr:
#	    print stderr
    if exitcode:
        print stderr
	sys.exit('\nkraken-build failed - see error message above\n')

def add_fasta_to_kraken_db(db_in_fasta, db_path='./', db_name='custom', v=0):
    """
    adds a fasta file (formatted according to kraken guidelines) to kraken db
    """
    
    import shlex, subprocess
    
    print "## Add custom data to kraken database ##\n"
    cmd="kraken-build --add-to-library %s --db %s/%s\n" %(db_in_fasta, db_path, db_name)
    
    print cmd
    cmdlist = shlex.split(cmd)
    proc = subprocess.Popen(cmdlist, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    exitcode = proc.returncode
#    print "EXITCODE: %s" %exitcode

#    stdout, stderr = subprocess.Popen(cmdlist, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate() # , stdout=subprocess.PIPE).communicate()
    if v:
	if stdout:
            print stdout
	if stderr:
	    print stderr
    if exitcode:
        print stderr
	sys.exit('\nkraken-build failed - see error message above\n')
        
def build_custom_kraken_db(db_path='./', db_name='custom', threads=1, jellyfish_hash_size=False, v=0):
    """
    adds a fasta file (formatted according to kraken guidelines) to kraken db
    """
    
    import shlex, subprocess
    
    print "## Build kraken database ##\n"
    cmd="kraken-build --build --threads %s --db %s/%s" %(threads, db_path, db_name)
    if jellyfish_hash_size:
	cmd+=' --jellyfish-hash-size %s' %jellyfish_hash_size
    
    print cmd
    cmdlist = shlex.split(cmd)
    proc = subprocess.Popen(cmdlist, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    exitcode = proc.returncode
#    print "EXITCODE: %s" %exitcode

#    stdout, stderr = subprocess.Popen(cmdlist, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate() # , stdout=subprocess.PIPE).communicate()
    if v:
	if stdout:
            print stdout
	if stderr:
	    print stderr
    if exitcode:
        print stderr
	sys.exit('\nkraken-build failed - see error message above\n')
    
def full_build_custom_kraken_db(db_in_fasta, db_path='.', db_name='KRAKEN_DB', threads=1, v=args.verbose, jellyfish_hash_size=False):
    
    print "\n### BUILD CUSTOM KRAKEN DATABASE ###\n"
    
    initiate_custom_kraken_db(db_path, db_name, v)
    add_fasta_to_kraken_db(db_in_fasta, db_path, db_name, v)
    build_custom_kraken_db(db_path, db_name, threads, jellyfish_hash_size, v)
    
def run_kraken(kraken_db, queries_fasta, threads=args.n_threads):
    """
    Run kraken
    """
    
    import shlex, subprocess
    import os.path

    print "## Running Kraken ##\n"
    cmd="kraken --threads %s --db %s %s\n" %(str(threads), kraken_db, queries_fasta)
    
    print cmd
    cmdlist = shlex.split(cmd)
    proc = subprocess.Popen(cmdlist, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    exitcode = proc.returncode

#    stdout, stderr = subprocess.Popen(cmdlist, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate() # , stdout=subprocess.PIPE).communicate()
    if exitcode:
        print stderr
	sys.exit('\nkraken failed - see error message above\n')
    
    if stderr:
        print stderr

    out = open('kraken.out','w')
    out.write(stdout)
    out.close()
    print "Written results to %s\n" %os.path.abspath('kraken.out')
    
    return os.path.abspath('kraken.out')

def kraken_filter(kraken_db, kraken_out, score_threshold=0):
    """
    Run kraken-filter. Infer confidence score for kraken results and or filter by score.
    """

    import shlex, subprocess
    import os.path

    print "## Inferring kraken confidence scores / filter by score ##\n"
    cmd="kraken-filter --db %s --threshold %s %s\n" %(kraken_db, score_threshold, kraken_out)

    print cmd
    cmdlist = shlex.split(cmd)
    proc = subprocess.Popen(cmdlist, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    exitcode = proc.returncode

#    stdout, stderr = subprocess.Popen(cmdlist, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate() # , stdout=subprocess.PIPE).communicate()
    if exitcode:
        print stderr
        sys.exit('\nkraken failed - see error message above\n')

    if stderr:
        print stderr

    out = open('kraken.filtered.score-'+str(score_threshold)+'.out','w')
    out.write(stdout)
    out.close()
    print "Written results to %s\n" %os.path.abspath('kraken.filtered.score-'+str(score_threshold)+'.out')

    return os.path.abspath('kraken.filtered.score-'+str(score_threshold)+'.out')

def parse_kraken_out(kraken_out):
    """
    The function parses the raw kraken output table into a dictionary
    """
    
    result_dict = {'hit':{}, 'nohit':[], 'format':'taxid'}
    
    in_table = open(kraken_out,'r')
    for l in in_table:
#        print l.strip()
        if l.startswith('U'):
            result_dict['nohit'].append(l.split("\t")[1])
        else:
#            print "HIT: %s" %l.split("\t")[1]
            result_dict['hit'][l.split("\t")[1]] = []
            result_dict['hit'][l.split("\t")[1]] = [l.split("\t")[2]]
        
    in_table.close()
    
    if len(result_dict['nohit']) == 0:
        del(result_dict['nohit'])
    
    return result_dict

def extract_taxid_list_from_result(dictionary):
    """
    The function extracts the values of a dictionary, which are themselves lists to a list
    """

    out_list = []

    for val in dictionary.values():
	out_list.extend(val)

    return out_list

def extract_queries_plus_rc(infile, good_list, bad_list, rc_list, out_prefix='pplacer'):
	"""
	Bins reads into two files according to lists containing the seqeuence ids.
	"""
	from Bio import SeqIO

	good_seqs = []
	bad_seqs = []
	handle = open(infile,'r')

	for r in SeqIO.parse(handle,'fasta'):
		if r.id in good_list:
			good_seqs.append(r)
		elif r.id in bad_list:
			bad_seqs.append(r)

	if rc_list:
		for r in good_seqs:
			if r.id in rc_list:
				r.seq = r.seq.reverse_complement()

	good_out = open(out_prefix+'.queries.fasta','w')
	bad_out = open(out_prefix+'.excluded.fasta','w')
	SeqIO.write(good_seqs, good_out, 'fasta')
	SeqIO.write(bad_seqs, bad_out, 'fasta')
	good_out.close()
	bad_out.close()

def parse_refpkg_json(refpkg):
	"""
	parse pplacer refpkg and return full path to alignment and hmm profile
	"""
	import json

	print "\nparsing pplacer refpkg: %s\n" %refpkg
	fh = open(refpkg+"/CONTENTS.json","r")
	refpkg_content = json.load(fh)
	profile = refpkg+'/'+refpkg_content['files']['profile']
	print "found hmm profile: %s" %profile
	aln_fasta = refpkg+'/'+refpkg_content['files']['aln_fasta']
	print "found alignment: %s" %aln_fasta

	return profile,aln_fasta

def run_hmmalign(queries, aln_fasta, profile, out='queries_to_profile.sto'):
	"""
	Run basic hmmalign to align sequences to alignment and output in *.sto format
	"""
	import shlex, subprocess
	
	cmd = "hmmalign -o %s --mapali %s %s %s" %(out, aln_fasta, profile, queries)
	print cmd
	cmdlist = shlex.split(cmd)
	stdout, stderr = subprocess.Popen(cmdlist, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate() # , stdout=subprocess.PIPE).communicate()
	print stdout	#this currently doesnt print anything, i.e. there is not stdout
	if stderr:
		print stderr
		sys.exit('some error from hmmalign')

def run_pplacer(refpkg, queries_to_profile, prefix='pplacer', v=0):
	"""
	Run pplacer
	"""

	cmd = "pplacer -c %s %s -p --keep-at-most 3 -o %s.jplace" %(refpkg, queries_to_profile, prefix)
	print cmd+'\n'
	cmdlist = shlex.split(cmd)
	stdout, stderr = subprocess.Popen(cmdlist, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate() # , stdout=subprocess.PIPE).communicate()
	if v:
		print stdout
	if stderr:
		print stderr
		sys.exit('some error from placer')

def parse_pplacer(jplace_file, result_dict): 
        """ 
        parse pplacer result (*.jplace) and write results to dictionary 
        """
        import json

        result_dict['hit'] = {}
        
        fh = open(jplace_file,"r")
        jplace = json.load(fh)
        for pp_placement in jplace['placements']:
#                print "PLACEMENT\n%s" %pp_placement
#                print pp_placement['p'][0][0]
#               print pp_queries['p'][0] #just consider the first placement
                placement = pp_placement['p'][0][0] #take the first placement for now
                for query in pp_placement['nm']:
                        query = str(query[0].replace('\\',''))
#                        print "%s\t%s" %(query,placement)
                        result_dict['hit'][query] = [placement]

        return result_dict
    

def assign_taxonomy_pplacer(pplacer_out, tax_dict, v=0):
    """
    finds taxonomic level of assignemt based on taxid and bins into dictionary
    """
    from collections import defaultdict
  
      
    levels = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'] #taxonomic levels of interest
    tax_count = {'kingdom':{}, 'phylum':{}, 'class':{}, 'order':{}, 'family':{}, 'genus':{}, 'species':{}}
    minimum = tax_dict["tax_id"].index("kingdom")
  
    print "\nInterpreting pplacer results and adjust to standard taxonomic levels (%s)\n" %", ".join(levels)
 
    if pplacer_out.has_key('nohit'):
        if not tax_count.has_key('nohit'):
            tax_count['nohit'] = {'nohit':[]}
        tax_count['nohit']['nohit'].extend(pplacer_out['nohit'])

    if pplacer_out.has_key('hit'):
        
        for query in pplacer_out['hit']:

            #find the index of the initial assignment
            index = tax_dict["tax_id"].index(tax_dict[pplacer_out['hit'][query][0]][1])

            #if already the initial assignemnt is above kingdom level
            if index < minimum:
                if v:
                    print "initial assignment above level kingdom for %s -> bin as 'nohit'" %query
                tax_count['nohit']['nohit'].append(query)
                continue

            #if the initial assignment hits one of the focal taxonomic level
            if tax_dict[pplacer_out['hit'][query][0]][1] in levels and tax_dict[pplacer_out['hit'][query][0]][2]:
                if v:
                    print "\ndirect assignment for %s -> %s" %(query, pplacer_out['hit'][query][0], tax_dict[pplacer_out['hit'][query][0]][2])
                if not tax_count[tax_dict[pplacer_out['hit'][query][0]][1]].has_key(pplacer_out['hit'][query][0]):
                    tax_count[tax_dict[pplacer_out['hit'][query][0]][1]][pplacer_out['hit'][query][0]] = []
                tax_count[tax_dict[pplacer_out['hit'][query][0]][1]][pplacer_out['hit'][query][0]].append(query)

            #if the initial assignment does not hit one of the focal taxonomic levels
            else:
                if v:
                        print "invalid assignment level for %s\t%s - %s - '%s'" %(query, tax_dict[pplacer_out['hit'][query][0]][1], pplacer_out['hit'][query][0], tax_dict[pplacer_out['hit'][query][0]][2])
                        print "\tsearching through higher taxonomic levels:"
                ok = 0
                
                #we'll count down, i.e. move taxonomic levels up one by one until we hit a level that is acceptable
                while index >= minimum and not ok:
                    index-=1
                    if tax_dict["tax_id"][index] in levels: #check if the current taxonomic level is a valid one

                        if tax_dict[pplacer_out['hit'][query][0]][index]: #check if there is actually a taxid availble for the lineage at this level
                            if v:
                                print "\tlevel '%s' is among the targets - taxid present -> OK!" %tax_dict["tax_id"][index]
                            ok = 1
                        else: #
                            if v:
                                print "\tlevel '%s' is among the targets, but no taxid available at this level -> moving on" %tax_dict["tax_id"][index]

                if ok:
                    taxid = tax_dict[pplacer_out['hit'][query][0]][index]

                    if v:
                        print "\tadjusted assignment for %s -> %s (%s)" %(query, taxid, tax_dict[taxid][2])
                    if not tax_count[tax_dict["tax_id"][index]].has_key(taxid):
                        tax_count[tax_dict["tax_id"][index]][taxid] = []
                    tax_count[tax_dict["tax_id"][index]][taxid].append(query)
                else:
                    if v:
                        print "\tCouldn't find a valid assignment for '%s'" %query
                    tax_count['nohit']['nohit'].append(query)

    ## cleanup - remove any taxonomic levels that were not hit
    for key in tax_count.keys():
        if not tax_count[key]:
            del(tax_count[key])
    
    return tax_count


###########


####START OF MAIN PROGRAM

print "\n%s\n" %DESCRIPTION

print '\n'+time.strftime("%c")+'\n'
print "%s\n" % (' '.join(sys.argv))

if args.update_taxonomy:
	update_taxonomy(v=args.verbose)

if args.email:
	Entrez.email = args.email
Entrez.email = check_email(mail = Entrez.email)#mail=Entrez.email)	

if args.read_crop:
	args.read_crop = " CROP:%s" %args.read_crop
else:
	args.read_crop = ""

if args.kraken_db:
	args.kraken_db = os.path.abspath(args.kraken_db)

#if args.blast or args.blast_xml or args.pplace:
if args.blast or args.kraken or args.pplace:
	taxonomy_db = check_for_taxonomy_db(tax_db=args.taxonomy_db)


if args.pplace:
	if not (args.refpkg or args.jplace):
		print "\nTo perform phylogenetic placement with pplacer metaBEAT currently expects a reference package to be specified via the --refpkg flag\n"
		sys.exit()
	if args.refpkg:
		if not os.path.isdir(args.refpkg):
			print "\nThe specified reference package does not seem to be valid\n"
			sys.exit()
		else:
			args.refpkg = os.path.abspath(args.refpkg) 
	if args.jplace:
		args.jplace = os.path.abspath(args.jplace)

	if not (args.blast or args.blast_xml):
		print "\nPhylogenetic placement currently requires a blast search first to determine the set of queries for which phylogenetic placement is attempted - add the --blast flag to your command\n"
		sys.exit()	

if args.product_length:
	args.product_length = '-M %s' %args.product_length
else:
	args.product_length = ''

if (args.min_ali_length > 1 or args.min_ali_length < 0) or (args.min_ident > 1 or args.min_ident < 0):
        print "\nThe options --min_ali_length and --min_ident expect values between 0 and 1.\n"
        sys.exit()

if args.REFlist and args.blast_db:
	print "\nPlease provide either a set of custom reference sequences OR a precompiled BLAST database\n"
	sys.exit()

if args.metadata:
	args.metadata = os.path.abspath(args.metadata)	

if args.blast_db:
	args.blast_db = os.path.abspath(args.blast_db)
	for f in glob.glob(args.blast_db+"*"):
		if not bl_db_extensions:
			print "ok - seems the precompiled BLAST database contains all necessary files\n"
			break
		for i in range(len(bl_db_extensions)):
			if f.endswith(bl_db_extensions[i]):
				del bl_db_extensions[i]
				break

	if bl_db_extensions:	#if I get to this point and there is a BLAST database extension that has not yet encountert then there is no file with this extension
		print "precompiled BLAST database is not intact\n"
		sys.exit()

elif args.REFlist:
	file_check(file_to_test=args.REFlist, optional_message="\nprovided reference file is not a valid file")
	
	refdata_files = [line.strip() for line in open(args.REFlist)]
	if not refdata_files:
		print "reference file is empty\n"
		sys.exit(0)
	for line in refdata_files:
		refdata = line.split("\t")
		if not len(refdata)==2:
			print "\nreference file is in incorrect format - we expect a tab delimited file:\n<filepath>\t<format>\n"
			sys.exit(0)
		else:
			if not refdata[1] in informats:
				print "\nthe format provided for %s is not valid\naccepted are: %s" % (refdata[0], informats.keys())
				sys.exit(0)
			elif not os.path.isfile(refdata[0]):
				print "\n%s is not a valid file\n" % refdata[0]
				sys.exit(0)
			else:
				references[refdata[0]] = refdata[1]

if args.querylist:
	print "\nParsing querylist file\n"
	data_formats = defaultdict(int)
	barcode_count = 0
	crop_count = 0
	if not os.path.isfile(args.querylist):
	        print "no valid query file supplied\n"
	        sys.exit(0)
	else:
	        query_files = [line.strip() for line in open(args.querylist)]
		if not query_files:
			print "query file is empty\n"
			sys.exit(0)
		for line in query_files:
			querydata = line.split("\t")
			if not len(querydata)>=3:
				print "\nquery file is in incorrect format - we expect a tab delimited file:\n<sample_ID>\t<format>\t<file>\t<optional file>"
				sys.exit(0)
			else:
				per_sample_query_format=""
				per_sample_barcodes=[]
				per_sample_query_files=[]
				crops = []
#				print "processing sample: %s" %querydata[0]
				if not informats.has_key(querydata[1]):
					print "query file format for sample %s is invalid" % querydata[0]
					per_sample_query_format=querydata[1]
					sys.exit(0) 
				for i in range(2,len(querydata)):
					if re.match("^[ACTG]*$", querydata[i]):
#						print "column %i looks like a barcode" %i
						per_sample_barcodes.append(querydata[i])
					elif re.match("^[0-9]*$", querydata[i]):
						crops.append(int(querydata[i]))
					else:
						if not os.path.isfile(querydata[i]):
							print "%s is not a valid file" %querydata[i]
							sys.exit(0)
						per_sample_query_files.append(os.path.abspath(querydata[i]))
				
				queries[querydata[0]]['format'] = querydata[1]
				queries[querydata[0]]['files'] = per_sample_query_files
				if crops:
					queries[querydata[0]]['crop_bases'] = crops
					crop_count+=1
				if per_sample_barcodes:
#					print "barcodes %s found for sample %s" %("-".join(per_sample_barcodes), querydata[0])
					queries[querydata[0]]['barcodes'] = per_sample_barcodes
					files_to_barcodes["|".join(per_sample_query_files)]["-".join(per_sample_barcodes)] = querydata[0]
					barcode_count+=1
				data_formats[queries[querydata[0]]['format']]+=1
#			print queries[querydata[0]]
#			print files_to_barcodes

	print "Number of samples to process: %s" %len(queries)
	print "Sequence input format: %s" %data_formats
	print "Barcodes for demultiplexing provided for %s samples" %barcode_count
	print "Cropping instructions provided for %s samples\n" %crop_count

	if args.PCR_primer and crops:
		print "\ncurrently you are requesting both PCR primer removal AND sequence clipping at the same time - note that sequences will be clipped first so primer removal may be redundant\n"
#		sys.exit()


	if args.PCR_primer:
		if not os.path.isfile(args.PCR_primer):
			print "PCR primer file is not valid"
			parser.print_help()
			sys.exit(0)
		args.PCR_primer = os.path.abspath(args.PCR_primer)
		#check for IUPAC and expand if necessary
		IUPAC = {'R': ['A','G'], 'W': ['A','T'], 'M': ['A','C'], 'S': ['C','G'], 'Y': ['C','T'], 'K': ['G','T']}
		t_seqs = list(SeqIO.parse(open(args.PCR_primer,'r'),'fasta'))
		for r in t_seqs:
			pos = []
			base = []
			states = []
#			print r
			for i in range(len(r.seq)):
				if not r.seq[i] in 'AGCTagct':
#					print r.seq[i]
					pos.append(i)
					if r.seq[i] == 'I':
						base.append('K')
					else:
						base.append(r.seq[i])

			if len(pos) > 0:
				states = list(product('01',repeat=len(pos)))
			
#				print len(states)
				state_count=0
				for s in states:
#					print s
#					print r.seq 
					new_seq = list(str(r.seq))
					for i in range(len(pos)):
#						print "This is ambiguous site %i: %i" %(i, pos[i])
#						print "This is the original base: %s" %base[i]
#						print "These are the options for this base: %s" %IUPAC[base[i]]
#						print "Will be set to %s (state: %s)" %(IUPAC[base[i]][int(s[i])], s[i])
						new_seq[pos[i]] = IUPAC[base[i]][int(s[i])]
					rec = SeqRecord(Seq("".join(new_seq)), id=r.id+'_v'+str(state_count), description=r.id+'_v'+str(state_count))
#					print rec
					primer_versions.append(rec)
#					primer_versions[r.id+'_v'+str(state_count)] = "".join(new_seq)
#					print ">%s_v%i" %(r.id, state_count)
#					print "".join(new_seq)
					state_count += 1

			else:
				primer_versions.append(r)

		##add reverse complements for all primers
		rc_versions = []
		for rec in primer_versions:
#			print rec
			new_rec = SeqRecord(rec.seq.reverse_complement(), id=rec.id+'_rc', description=rec.id+'_rc')
			rc_versions.append(new_rec)
		primer_versions.extend(rc_versions)
		
	if args.trim_adapter:
		args.trim_adapter = os.path.abspath(args.trim_adapter)
		if not os.path.isfile(args.trim_adapter):
			print "adapter file is not a valid file"
			parser.print_help()
			sys.exit(0)

if args.BIOM_input:
	print "\nparsing BIOM table\n"
	args.BIOM_input = os.path.abspath(args.BIOM_input)
	BIOM_tables_per_method['OTU_denovo'] = parse_BIOM_denovo(table=args.BIOM_input)

#check if gb_to_taxid file is there and read if yes
if os.path.isfile(args.gb_to_taxid):
	print "\nfound gb_to_taxid file at %s" %args.gb_to_taxid,
	rw_gb_to_taxid_dict(dictionary=gb_to_taxid_dict, name=args.gb_to_taxid, mode='r')
	

#print references
if args.REFlist:
	print "\n######## PROCESSING REFERENCE DATA ########\n"

	for reffile, refformat in references.items():
		seqs=list(SeqIO.parse(reffile,refformat, alphabet=generic_dna))
		print "\nprocessing %s (containing %i records)" % (reffile,len(seqs))
		for i in reversed(range(len(seqs))):
#			print seqs[i].description
			taxid = ''
			if not seqs[i].name in records:	#records list keeps track of all records.
				records[seqs[i].name] = 1
			else:				#if a record id occurs again, say in the case of custom fasta sequences the record id is the first name in the fasta header, e.g. the Genus name
				records[seqs[i].name] += 1
				seqs[i].name += "_%s" % records[seqs[i].name]
	
			if not seqs[i].features:		#if the record does not have a source feature yet. Will happen for records originating from custom fasta files
#				print "the record has no features\n"	
				seqs[i].features.append(SeqFeature(FeatureLocation(0,len(seqs[i].seq),strand=1), type="source"))	#create source feature
				seqs[i].features[0].qualifiers['original_desc'] = seqs[i].description
				
			if not seqs[i].features[0].qualifiers.has_key('organism'):	#if the record does not have an organism key. Again will happen for custom fasta sequences
#				print seqs[i].description
				seqs[i].features[0].qualifiers['organism'] = ["%s %s" % (seqs[i].description.split(" ")[0], seqs[i].description.split(" ")[1])] #add an organism key and make it the first 2 words in the record description will are per default the first 2 words in the fasta header
#				seqs[i].features[0].qualifiers['organism'] = ["%s" % seqs[i].description.split(" ")[0]]	#add an organism key and make it the first word in the record description will are per default the first 2 words in the fasta header
#				sp_search = re.compile('^sp$|^sp[.]$|^sp[0-9]+$')
#				if not sp_search.match(seqs[i].description.split(" ")[1]):
#					seqs[i].features[0].qualifiers['organism'] = ["%s %s" % (seqs[i].features[0].qualifiers['organism'][0], seqs[i].description.split(" ")[1])]	#if the second word in the description is not 'sp' or 'sp.', i.e. it is unknown, then add also the second word to the organism field
	
			if not seqs[i].annotations.has_key('date'):	#set the date in the record if it does not exist
				seqs[i].annotations['date'] = date
				
			if not seqs[i].features[0].qualifiers['organism'][0] in reference_taxa:	#if the organism name has not yet been encountert the following lines will find the corresponding taxid, add it to the record and store it in the reference_taxa dictionary
#				print "seeing %s for the first time" %seqs[i].features[0].qualifiers['organism'][0]
				if seqs[i].features[0].qualifiers.has_key('db_xref'):	#there may be more than one 'db_xref' qualifiers, e.g. sometimes there are BOLD ids also in 'db_xref'
					for j in range(len(seqs[i].features[0].qualifiers['db_xref'])):
						if seqs[i].features[0].qualifiers['db_xref'][j].startswith('taxon'):
							taxid = seqs[i].features[0].qualifiers['db_xref'][j].split(":")[-1]
#							print "found valid taxid in record: %s" %taxid
							break
	
				if not taxid or args.rec_check:
					done = False
       		                    	while not done:
	       	                 		try:
							handle = Entrez.esearch(db="Taxonomy", term=seqs[i].features[0].qualifiers['organism'][0])	#search the taxonomy database for the taxon by organism name
	                                		done = True
	                        		except:
	                                		print "connection closed - retrying entrez for query '%s'.." %seqs[i].features[0].qualifiers['organism'][0]
	
					taxon = Entrez.read(handle)
					if not taxon["IdList"]:	#if the search term has not yielded any result
						print "\nWARNING: \"%s\" in record '%s' was not recognized as a valid taxon name on Genbank." % (seqs[i].features[0].qualifiers['organism'][0], seqs[i].name)
						print "I'll try if just \"%s\" is valid" %seqs[i].features[0].qualifiers['organism'][0].split(" ")[0]
						done = False
			                        while not done:
			                        	try:
								handle = Entrez.esearch(db="Taxonomy", term=seqs[i].features[0].qualifiers['organism'][0].split(" ")[0])
	                                			done = True
	                        			except:
	                                			print "connection closed - retrying entrez for query '%s'.." %seqs[i].features[0].qualifiers['organism'][0].split(" ")[0]
	
						taxon = Entrez.read(handle)
						if taxon['IdList']: #if the search with just the first term of the taxon name yielded a result
							tax_rank = ''
							done = False
				                        while not done:
					                        try:
									handle = Entrez.efetch(db="Taxonomy", id=taxon['IdList'][0]) #fetch the taxonomy to this taxid
					                                done = True
					                        except:
	                                				print "connection closed - retrying entrez for query '%s'.." %taxon['IdList'][0]
	
							recs = Entrez.read(handle)
							tax_rank = recs[0]["Rank"]
							print "ok - found it with taxid \"%s\" at taxonomic rank \"%s\"" %(taxon['IdList'][0],tax_rank)
							print "I am interpreting \"%s\" as a valid species name and will assign a dummy taxid to it" %seqs[i].features[0].qualifiers['organism'][0]
#							taxids[taxon['IdList'][0]]+= 1
							if denovo_taxa.has_key(taxon['IdList'][0]):
								denovo_taxa[taxon['IdList'][0]].append(seqs[i].features[0].qualifiers['organism'][0]+'|denovo'+str(denovo_count))  #seqs[i].features[0].qualifiers['organism'][0])
							else:
								denovo_taxa[taxon['IdList'][0]] = [seqs[i].features[0].qualifiers['organism'][0]+'|denovo'+str(denovo_count)] #seqs[i].features[0].qualifiers['organism'][0]]  #add species that has not been found in taxdump. The value is the taxid of the genus
							reference_taxa[seqs[i].features[0].qualifiers['organism'][0]] = taxon['IdList'][0] #"denovo"+str(denovo_count)
							seqs[i].features[0].qualifiers['db_xref'] = ['taxon:denovo'+str(denovo_count)]
							denovo_count += 1
							print "WARNING: WILL NOT ADD ANY INFORMATION TO seq_info LIST AT THIS STAGE - NEEDS TO BE FIXED IF PHYLOGENETIC PLACEMENT IS PLANNED"
						else:
							print "WARNING: Nothing found. Something's wrong with this record - typo? I'll put the sequence record in the file 'skipped_seqs.fasta' in case you want to fix it."
							skipped_ref_seqs.append(seqs.pop(i))
						
		
					else:
						seqs[i].id = seqs[i].name
						taxid = taxon["IdList"][0]	#the taxid is the first element in the dictionary that is returned from Entrez
						if seqs[i].features[0].qualifiers.has_key('db_xref'):
							for j in range(len(seqs[i].features[0].qualifiers['db_xref'])):
								if seqs[i].features[0].qualifiers['db_xref'][j].startswith('taxon'):
									seqs[i].features[0].qualifiers['db_xref'][j] = "taxon:" + taxid	#update the db_rxref qualifier with the correct taxid
						else:
							seqs[i].features[0].qualifiers['db_xref'] = ["taxon:" + taxid]
#						taxids[taxid] += 1
#						print "adding %s with taxid: %s to the dictionary" %(seqs[i].features[0].qualifiers['organism'][0], taxid)
						reference_taxa[seqs[i].features[0].qualifiers['organism'][0]] = taxid
				else:
					seqs[i].id = seqs[i].name
					reference_taxa[seqs[i].features[0].qualifiers['organism'][0]] = taxid
#					taxids[taxid] += 1
					
			else:
				taxid = reference_taxa[seqs[i].features[0].qualifiers['organism'][0]]
				if denovo_taxa.has_key(taxid):
					for den in denovo_taxa[taxid]:
						if den.split("|")[0] == seqs[i].features[0].qualifiers['organism'][0]:
							taxid = den.split("|")[1]
							break
	
				seqs[i].id = seqs[i].name
				if seqs[i].features[0].qualifiers.has_key('db_xref'):
					for j in range(len(seqs[i].features[0].qualifiers['db_xref'])):
						if seqs[i].features[0].qualifiers['db_xref'][j].startswith('taxon'):
							seqs[i].features[0].qualifiers['db_xref'][j] = "taxon:" + taxid	#update the db_rxref qualifier with the correct taxid
				else:
					seqs[i].features[0].qualifiers['db_xref'] = ["taxon:" + taxid]
	
#				print "Have seen %s before. The taxid is: %s" %(seqs[i].features[0].qualifiers['organism'][0], reference_taxa[seqs[i].features[0].qualifiers['organism'][0]])
	
			
			if len(taxid) > 0:
#				taxids[taxid] += 1
				current = '"%s","%s","%s","%s",0' % (seqs[i].id, seqs[i].id, taxid, seqs[i].features[0].qualifiers['organism'][0]) #producing a string for the current record for the seq_info file that is needed for the reference package required by pplacer
				seq_info.extend([current])	#add the string to the seq_info array
	
		if skipped_ref_seqs:
			SKIPPED = open("skipped_seqs.fasta","w")
			for rec in skipped_ref_seqs:
				outstring = ">%s\n%s" %(rec.name, rec.seq)
				SKIPPED.write(outstring + "\n")
				
			SKIPPED.close()
				
		all_seqs.extend(seqs)
	
		if args.gb_out:
			gb_out=open(args.gb_out,"w")
#			gb_out.write(seqs.format("genbank"))
			SeqIO.write(all_seqs, gb_out, "genbank")
#fas_out.write(unknown_seqs_dict[ID].format("fasta"))
			gb_out.close()
	
		print "\ntotal number of valid records: %i\n" % len(all_seqs)
	
if args.seqinfo:
	print "write out seq_info.csv\n"
	f = open("seq_info.csv","w")
	for elem in seq_info:
#		print elem
		f.write(elem + "\n")
	f.close()

if args.taxids:
	write_taxids(reference_taxa.values())

if args.fasta:	#this bit writes out the sequences that will become the reference dataset
	write_out_refs_to_fasta(ref_seqs=all_seqs, ref_taxids=reference_taxa)

if args.blast_xml:
	args.blast_xml = os.path.abspath(args.blast_xml)

print '\n'+time.strftime("%c")+'\n'
querycount = defaultdict(int)

if files_to_barcodes:
	print "\n### DEMULTIPLEXING ###\n"
#	print "Barcodes detected - initialize demultiplexing"
	if not os.path.exists('demultiplexed'):
		os.makedirs('demultiplexed')
		
	os.chdir('demultiplexed')
	to_append=''
	for lib in files_to_barcodes.keys():
		data = lib.split("|")
		compr_mode = ""
		print "assessing basic characteristics"
#		print files_to_barcodes[lib]
		barcodes_num = len(files_to_barcodes[lib].keys()[0].split("-"))
		if barcodes_num == 2:
#			print "interpreting barcodes as forward and reverse inline"
			barcode_mode = "--inline_inline"
		elif barcodes_num == 1:
#			print "interpreting barcodes as forward inline only"
			barcode_mode = "--inline_null"
				
		if lib.endswith(".gz"):
			compr_mode = "gzfastq"
#			print "set mode to %s" %compr_mode
		else:
			compr_mode = "fastq"
#			print "set mode to %s" %compr_mode
#		print "create temporal barcode file"
		if args.bc_dist:
			to_append=' -r --barcode_dist %s' %args.bc_dist


		BC = open("barcodes","w")
		for sample in files_to_barcodes[lib].keys():
#			print "This is sample: %s" %sample
			bc = sample.split("-")
#			print bc
			if len(bc) != barcodes_num:
				print "expecting %i barcode(s) for all samples - something's wrong with sample %s" %(barcodes_num, sample)
				sys.exit()
			BC.write("\t".join(bc)+"\n")

		BC.close()
		if len(data)>1:
			print "data comes as paired end - ok"
			demultiplex_cmd="process_shortreads -1 %s -2 %s -i %s -o . -b barcodes %s -y fastq" %(data[0], data[1], compr_mode, barcode_mode)
			
			if to_append:
				demultiplex_cmd+=to_append

			print demultiplex_cmd
			cmdlist = shlex.split(demultiplex_cmd)
			cmd = subprocess.Popen(cmdlist, stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
			stdout,stderr = cmd.communicate() #[0]
			if args.verbose:
				print stdout
			if stderr:
				print stderr

#			print "renaming results"
			for sample in files_to_barcodes[lib].keys():
#				print files_to_barcodes[lib][sample]
				new_file_list=[]
				for a in range(1,3):
					old = "sample_%s.%i.fq" %(sample,a)
#					print old
					new = "%s_%i.fastq" %(files_to_barcodes[lib][sample],a)
#					print new
					os.rename(old,new)
					os.remove("sample_%s.rem.%i.fq" %(sample,a))
					new_file_list.append(os.path.abspath('.')+"/"+new)
					
				queries[files_to_barcodes[lib][sample]]['files'] = new_file_list
		elif len(data)==1:
			print "data comes in single end - ok"
			print "currently not supported"
			sys.exit()
	os.chdir('../')
#print queries

for queryID in sorted(queries):
#	print queries
#	query_count += 1
	print "\n##### processing query ID: %s #####\n" % (queryID)
	if not os.path.exists(queryID):
		os.makedirs(queryID)
		
#	print queries[queryID]['files']
		
	os.chdir(queryID)

	if (queries[queryID]['format']=="fastq"):
		#determine total read number per sample
		read_stats[queryID]['total'] = int(0)
		for f in queries[queryID]['files']:
			if f.endswith('gz'):
				t_seqs = list(SeqIO.parse(gzip.open(f, 'rb'),queries[queryID]['format']))
			else:
				t_seqs = list(SeqIO.parse(open(f, 'r'),queries[queryID]['format']))
			read_stats[queryID]['total'] += len(t_seqs)
		
#		print "Total number of reads for sample %s: %i" %(queryID, read_stats[queryID][0])

		print "\n### READ QUALITY TRIMMING ###\n"
		if len(queries[queryID]['files'])==2:
			trimmomatic_path="java -jar /usr/bin/trimmomatic-0.32.jar"
			print "\ntrimming PE reads with trimmomatic"
			trimmomatic_exec=trimmomatic_path+" PE -threads %i -phred%i -trimlog trimmomatic.log %s %s %s_forward.paired.fastq.gz %s_forward.singletons.fastq.gz %s_reverse.paired.fastq.gz %s_reverse.singletons.fastq.gz ILLUMINACLIP:%s:5:5:5 TRAILING:%i LEADING:%i SLIDINGWINDOW:%i:%i%s MINLEN:%i" % (args.n_threads, args.phred, queries[queryID]['files'][0], queries[queryID]['files'][1], queryID, queryID, queryID, queryID, args.trim_adapter, args.trim_qual, args.trim_qual, args.trim_window, args.trim_qual, args.read_crop, args.trim_minlength)
			trimmed_files=[queryID+'_forward.paired.fastq.gz', queryID+'_forward.singletons.fastq.gz', queryID+'_reverse.paired.fastq.gz', queryID+'_reverse.singletons.fastq.gz']
		elif len(queries[queryID]['files'])==1:
			print "\ntrimming SE reads with trimmomatic"
			trimmomatic_exec= trimmomatic_path+" SE -threads %i -phred%i %s %s_trimmed.fastq.gz ILLUMINACLIP:%s:5:5:5 TRAILING:%i LEADING:%i SLIDINGWINDOW:%i:%i%s MINLEN:%i" % (args.n_threads, args.phred, queries[queryID]['files'][0], queryID, args.trim_adapter, args.trim_qual, args.trim_qual, args.trim_window, args.trim_qual, args.read_crop, args.trim_minlength)
			trimmed_files=[queryID+'_trimmed.fastq.gz']

		trimmomatic="%s" % trimmomatic_exec
		cmdlist = shlex.split(trimmomatic)
		cmd = subprocess.Popen(cmdlist, stdout=subprocess.PIPE, stderr=subprocess.PIPE)	
		stdout,stderr = cmd.communicate() #[0]
		print trimmomatic
		if args.verbose:
			print stdout
		if stderr:
			print stderr

		if 'crop_bases' in queries[queryID]:
			to_clip = trimmed_files[:]
			print "\n### READ CLIPPING ###\n"
			for i in range(len(to_clip)):
				f = to_clip[i]
				if i > len(queries[queryID]['crop_bases'])/2:
					c = queries[queryID]['crop_bases'][1]
				else:
					c = queries[queryID]['crop_bases'][0]
					
                       		if f.endswith('gz'):
                               		t_seqs = list(SeqIO.parse(gzip.open(f, 'rb'),queries[queryID]['format']))
                       		else:
                               		t_seqs = list(SeqIO.parse(open(f, 'r'),queries[queryID]['format']))
				print "Cropping the first %s bases off reads in file: %s" %(c, f)
				for j in range(len(t_seqs)):
					t_seqs[j] = t_seqs[j][c:]
				
				OUT=gzip.open(f.replace(queries[queryID]['format'], "headclipped."+queries[queryID]['format']), 'wb')
				SeqIO.write(t_seqs, OUT, "fastq")
				OUT.close()
				trimmed_files[i] = f.replace(queries[queryID]['format'], "headclipped."+queries[queryID]['format']) 

			
		if args.PCR_primer:
			temp_primer_out = open("temp_primers.fasta","w")
			SeqIO.write(primer_versions, temp_primer_out, "fasta")
			temp_primer_out.close()

			for f in trimmed_files:
				cmd = "zcat %s | fastx_reverse_complement -Q %i | gzip > %s" %(f, args.phred, f+'.rc.gz')
				print cmd
				cmdlist = shlex.split(cmd)
				cmd = subprocess.call(cmd, shell=True)

			print "\nCLIP PRIMER SEQUENCES\n"
			trimmomatic_exec=trimmomatic_path+" PE -threads %i -phred%i %s %s temp_forward.paired.fastq.gz temp_forward.singletons.fastq.gz temp_reverse.paired.fastq.gz temp_reverse.singletons.fastq.gz ILLUMINACLIP:temp_primers.fasta:2:5:5 MINLEN:%i" % (args.n_threads, args.phred, trimmed_files[0]+'.rc.gz', trimmed_files[2]+'.rc.gz', args.trim_minlength)

			trimmomatic="%s" % trimmomatic_exec
			print trimmomatic
			cmdlist = shlex.split(trimmomatic)
			cmd = subprocess.Popen(cmdlist, stdout=subprocess.PIPE, stderr=subprocess.PIPE)	
			stdout,stderr = cmd.communicate() #[0]
			print stdout
			print stderr

			cmd = 'zcat *singletons.fastq.gz.rc.gz | gzip > singletons_rc.fastq.gz'
			print cmd
			cmdlist = shlex.split(cmd)
			cmd = subprocess.call(cmd, shell=True)
			
			trimmomatic_exec=trimmomatic_path+" SE -threads %i -phred%i singletons_rc.fastq.gz singletons_rc_clipped.fastq.gz ILLUMINACLIP:temp_primers.fasta:2:5:10 MINLEN:%i" % (args.n_threads, args.phred, args.trim_minlength)

			trimmomatic="%s" % trimmomatic_exec
			print trimmomatic
			cmdlist = shlex.split(trimmomatic)
			cmd = subprocess.Popen(cmdlist, stdout=subprocess.PIPE, stderr=subprocess.PIPE)	
			stdout,stderr = cmd.communicate() #[0]
			print stdout
			print stderr


			cmd = "zcat temp_forward.paired.fastq.gz | fastx_reverse_complement -Q %i | gzip > %s" %(args.phred, trimmed_files[0])
			print cmd
			cmdlist = shlex.split(cmd)
			cmd = subprocess.call(cmd, shell=True)

			cmd = "zcat temp_reverse.paired.fastq.gz | fastx_reverse_complement -Q %i | gzip > %s" %(args.phred, trimmed_files[2])
			print cmd
			cmdlist = shlex.split(cmd)
			cmd = subprocess.call(cmd, shell=True)

			cmd = "zcat singletons_rc_clipped.fastq.gz temp_forward.singletons.fastq.gz | fastx_reverse_complement -Q %i | gzip > %s" %(args.phred, trimmed_files[1])
			print cmd
			cmdlist = shlex.split(cmd)
			cmd = subprocess.call(cmd, shell=True)

			cmd = "zcat temp_reverse.singletons.fastq.gz | fastx_reverse_complement -Q %i | gzip > %s" %(args.phred, trimmed_files[3])
			print cmd
			cmdlist = shlex.split(cmd)
			cmd = subprocess.call(cmd, shell=True)

#			for f in queries[queryID]['files']:
#				cmd = "cat %s | fastx_reverse_complement -Q %i |fastx_clipper -Q %i -n %s | fastx_reverse_complement -Q %i > %s" %(f, args.phred, args.phred, primer_clip_string,  args.phred, f+'.no_primer') 
#				print cmd
#				cmdlist = shlex.split(cmd)
#				cmd = subprocess.call(cmd, shell=True)

		if len(trimmed_files)==4:
			#determining read counts
			read_stats[queryID]['trimmed-total'] = int(0)
			read_stats[queryID]['trimmed-pe'] = int(0)
			read_stats[queryID]['trimmed-orphans'] = int(0)
			for s in trimmed_files:
				t_seqs = list(SeqIO.parse(gzip.open(s,'rb'),'fastq'))
				read_stats[queryID]['trimmed-total'] += len(t_seqs)
				if 'paired' in s:
					read_stats[queryID]['trimmed-pe'] += len(t_seqs)
					
				elif 'singleton' in s:				
					read_stats[queryID]['trimmed-orphans'] += len(t_seqs)



			if args.merge:
				print "\n### MERGING READ PAIRS ###\n"
				print "merging paired-end reads with flash\n"
				cmd="flash %s %s %s -t %i -p %i -o %s -z" % ( trimmed_files[0], trimmed_files[2], args.product_length, args.n_threads, args.phred, queryID)
				print cmd
				cmdlist = shlex.split(cmd)
				cmd = subprocess.Popen(cmdlist, stdout=subprocess.PIPE, stderr=subprocess.PIPE)	
				stdout,stderr = cmd.communicate() #[0]
				print stdout
				print stderr
				trimmed_files[0]= queryID+'.notCombined_1.fastq.gz'
				trimmed_files[2]= queryID+'.notCombined_2.fastq.gz'
				trimmed_files.insert(0, queryID+'.extendedFrags.fastq.gz')

				#change names of extended reads
				ext_seqs=[]
				t_seqs = list(SeqIO.parse(gzip.open(queryID+'.extendedFrags.fastq.gz', 'rb'),'fastq'))
				read_stats[queryID]['merged'] = 2*len(t_seqs)
				for record in t_seqs:
					record.id = record.id+"_ex"
					record.description = record.id
					ext_seqs.append(record)
					
				ex_out = gzip.open("ex_temp.fastq.gz","wb")
				SeqIO.write(ext_seqs, ex_out, "fastq")
				ex_out.close()
				os.rename("ex_temp.fastq.gz", queryID+'.extendedFrags.fastq.gz')						
				
#			print read_stats
			print "\n### Fastq to fasta (plus filtering for length, merged-only, etc. if applicable) ###\n"

			if args.merged_only:
				print "\nkeeping only merged reads for subsequent analyses"
				keep_only(inlist=trimmed_files, pattern_list=['.extendedFrags.fastq.gz'])

			if args.forward_only:
				print "\nkeeping only forward sequences (non merged forward and merged reads) for subsequent analyses"
				keep_only(inlist=trimmed_files, pattern_list=['.extendedFrags.fastq.gz', '_forward', 'notCombined_1'])
			
			survived = concat_and_filter_by_length(inlist=trimmed_files, outfile=queryID+'_trimmed.fasta', excludefile=queryID+'_excluded.fasta', form='fasta', length=args.length_filter, devi=args.length_deviation)
			print "\n%i sequences survived filtering\n" %survived

		else:
			survived = concat_and_filter_by_length(inlist=trimmed_files, outfile=queryID+'_trimmed.fasta', excludefile=False, form='fasta', length=0)

		#cleanup
		cont = os.listdir(".")
		to_remove = []
		for f in cont:
			if f.startswith('temp'):
				to_remove.append(f)
			
			if f.endswith('.gz.rc.gz') or f.endswith('_rc_clipped.fastq.gz') or f.endswith('_rc.fastq.gz'):
				to_remove.append(f)
			
		for f in to_remove:
			os.remove(f) #trimmed.fasta")
			

		queryfile=os.path.abspath('.')+"/"+queryID+'_trimmed.fasta'

	elif queries[queryID]['format']=="fasta":
		queryfile = queries[queryID]['files'][0]
		print "\nWARNING - EXPECTING input data to come in SINGLE FASTA FILE\n"


		#getting number of sequences to be queried		

	elif queries[queryID]['format'] == "uc":
		to_check=queries[queryID]['files'][0][:-3]+"_centroids.fasta"
		print("\nprovided uc file as input - expecting to find the file: "+to_check+" - checking .. "),
#		print("\nprovided uc file as input - expecting to find the file: "+queryID+"_queries.fasta - checking .. "),
		
		if os.path.exists(to_check):
			print "ok!"
			queryfile = os.path.abspath('.')+"/"+queryID+"_centroids.fasta"
			queries[queryID]['uc'] = os.path.abspath('.')+"/"+queryID+".uc"

	total_queries = len(list(SeqIO.parse(queryfile,'fasta')))

	if args.cluster or 'uc' in queries[queryID]['format']:

		print "\n### CLUSTERING ###\n"

		if args.cluster:
			##running clustering
			print "\nclustering using vsearch"
			vsearch_cluster(infile=queryfile, cluster_match=args.clust_match, threads=args.n_threads, sampleID=queryID)
		
			queries[queryID]['uc'] = os.path.abspath('.')+"/"+queryID+'.uc'
			queryfile = os.path.abspath('.')+"/"+"%s_centroids.fasta" %queryID #"../%s_centroids.fasta" % queryID
	
			#read in the results from the clustering, i.e. number of reads per cluster, later I could read in the actual read ids to allow for retrievel of reads assigned to a given taxon
#			all_clust_count = int(0)
		else:
			print "using existing clustering result: %s" %queries[queryID]['uc']


		if queries[queryID]['uc']:
	
			print "\nparse vsearch uc file\n"
			cluster_counts = {}
			cluster_reads = defaultdict(list)
			total_queries = parse_vsearch_uc(fil=queryID+".uc", cluster_counts=cluster_counts, extract_reads=args.extract_all_reads, cluster_reads=cluster_reads)

			
#			print cluster_reads.keys()
	
#			total_queries = querycount[queryID]
			total_clusters = len(cluster_counts)
		
			if args.clust_cov>1:
				print "\nreduce cluster files - minimum coverage: %s\n" %args.clust_cov
				querycount[queryID] = filter_centroid_fasta(centroid_fasta=queryfile, m_cluster_size=args.clust_cov, cluster_counts=cluster_counts, sampleID=queryID, v=args.verbose)		
#		for ID in unknown_seqs_dict.keys():
#			if cluster_counts.has_key(ID):
#				unknown_seqs_dict[ID].description = "%s|%s|%s|%.2f" %(queryID, unknown_seqs_dict[ID].id, cluster_counts[unknown_seqs_dict[ID].id], float(cluster_counts[unknown_seqs_dict[ID].id])/querycount[queryID]*100)
##			print unknown_seqs_dict[ID].description

			unknown_seqs_dict = SeqIO.to_dict(SeqIO.parse(queryfile,'fasta'))
			if not querycount[queryID]:
				querycount[queryID] = len(unknown_seqs_dict)
			for ID in unknown_seqs_dict.keys():
	                	if cluster_counts.has_key(ID):
	                        	unknown_seqs_dict[ID].description = "%s|%s|%s|%.2f" %(queryID, unknown_seqs_dict[ID].id, cluster_counts[unknown_seqs_dict[ID].id], float(cluster_counts[unknown_seqs_dict[ID].id])/querycount[queryID]*100)
#	                       print unknown_seqs_dict[ID].description

			print "vsearch processed %i sequences and identified %i clusters (clustering threshold %.2f) - %i clusters (representing %i sequences) passed the filter criteria (minimum of %i sequences per cluster) and will be used in subsequent analyses\n" % (total_queries, total_clusters, float(args.clust_match), len(cluster_counts), querycount[queryID], args.clust_cov)
		
			read_stats[queryID]['clusters_total'] = total_clusters
			read_stats[queryID]['cluster_thres'] = args.clust_match
			read_stats[queryID]['clusters_min_cov'] = args.clust_cov
			read_stats[queryID]['cluster_above_thres'] = len(cluster_counts)
			read_stats[queryID]['queries'] = querycount[queryID]



	else:
		unknown_seqs_dict = SeqIO.to_dict(SeqIO.parse(queryfile,'fasta')) #read in query sequences, atm only fasta format is supported. later I will check at this stage if there are already sequences in memory from prior quality filtering
		querycount[queryID] += len(unknown_seqs_dict)
		print "The queryfile contains %i query sequences.\n" %querycount[queryID]
		
		cluster_counts = {}
		cluster_reads = defaultdict(list)

		for ID in unknown_seqs_dict.keys():
#			print sequence
			cluster_counts[unknown_seqs_dict[ID].description] = 1
			cluster_reads[unknown_seqs_dict[ID].description] = [unknown_seqs_dict[ID].description]
#			unknown_seqs_dict[ID].description = "%s|%s|%s|%.2f" %(queryID, unknown_seqs_dict[ID].id, cluster_counts[unknown_seqs_dict[ID].id], float(cluster_counts[unknown_seqs_dict[ID].id])/querycount[queryID]*100)
			unknown_seqs_dict[ID].description = "%s|%s|%s|%.2f" %(queryID, unknown_seqs_dict[ID].id, cluster_counts[unknown_seqs_dict[ID].id], float(cluster_counts[unknown_seqs_dict[ID].id])/querycount[queryID]*100)
		read_stats[queryID]['queries'] = querycount[queryID]


	if not queryfile.split("/")[-1] == queryID+'_queries.fasta':  
		shutil.copy(queryfile, queryID+'_queries.fasta')
	queryfile = os.path.abspath('.')+"/"+queryID+'_queries.fasta' 


	queries[queryID]['queryfile'] = queryfile
	queries[queryID]['cluster_counts'] = dict(cluster_counts)
	queries[queryID]['cluster_reads'] = dict(cluster_reads)
#	print "\n%s\n" %queries[queryID]
	
	os.chdir('../') #leave directory for query

	if not args.read_stats_off:
		print "WRITING BASIC READ STATS TO FILE\n"

		#start read stats output
		if not read_counts_out: #create the file (overwrite existing) and add header line
			read_counts_out = open("./"+args.output_prefix+"_read_stats.csv","w")
			outstring = "sample,"+",".join(read_metrics)
			read_counts_out.write(outstring+"\n")
			read_counts_out.close()
			
		read_counts_out = open("./"+args.output_prefix+"_read_stats.csv","a")
		outstring = queryID
		for m in read_metrics:
			if read_stats[queryID].has_key(m):
				outstring += ','+str(read_stats[queryID][m])
			else:
				outstring += ',NA'
		read_counts_out.write(outstring+"\n")
		read_counts_out.close()
	


############



if not 'OTU_denovo' in BIOM_tables_per_method:

	cluster_check = 0
	for q in queries.keys():
		if 'cluster_counts' in queries[queryID]:
			cluster_check+=1

	if not cluster_check == len(queries):
		print "\nAborted before global clustering\n"
		sys.exit()

	if not os.path.exists('GLOBAL'):
		os.makedirs('GLOBAL')
	os.chdir('GLOBAL')

	print "\n### GLOBAL CLUSTERING ###\n"

	global_queries = os.getcwd()+'/global_queries.fasta'
	global_centroids = os.getcwd()+'/global_centroids.fasta'
	global_uc = os.getcwd()+'/global.uc'

	#print global_queries

	print "\nConcatenating for global clustering\n"
	concatenate_for_global_clustering(queries_dict=queries, out=global_queries)

	print "\nPerforming global clustering\n"
	vsearch_cluster(infile=global_queries, cluster_match=args.clust_match, threads=args.n_threads, sampleID='global')

	global_cluster_count=0
	for i in SeqIO.parse(open(global_centroids), 'fasta'):
		global_cluster_count += 1

	print "\nparse vsearch uc file\n"
	global_cluster_counts = {}
	global_cluster_reads = defaultdict(list)
	parse_vsearch_uc(fil=global_uc, cluster_counts=global_cluster_counts, extract_reads=True, cluster_reads=global_cluster_reads)

	BIOM_tables_per_method['OTU_denovo'] = global_uc_to_biom(clust_dict=global_cluster_reads, query_dict=queries)

	print "\nwriting denovo OTU table\n"

	out=open(args.output_prefix+"-OTU-denovo.biom","w")
	BIOM_tables_per_method['OTU_denovo'].to_json('metaBEAT v.'+VERSION, direct_io=out)
	out.close()

	out=open(args.output_prefix+"-OTU-denovo.tsv","w")
	out.write(BIOM_tables_per_method['OTU_denovo'].to_tsv()) #to_json('generaged by test', direct_io=out)
	out.close()
else: #I would get there if I gave the denovo BIOM table via the command line
	print "\ndenovo OTU table provided - checking if all other files that are required for taxonomic assignment are present\n" 
	print "expecting to find query sequences .."
	#check if global clustering result is there
	if args.g_queries:
		global_centroids = os.path.abspath(args.g_queries)
		print "manually specified query file via the '--g_queries' flag: %s" %global_centroids 
	elif os.path.isfile("/".join(args.BIOM_input.split("/")[:-1])+"/global_centroids.fasta"):
		global_centroids = "/".join(args.BIOM_input.split("/")[:-1])+"/global_centroids.fasta"
		print "query sequences detected automatically: %s" %global_centroids
	else:
		print "Can't find query sequences.\nYou may use the '--g_queries' flag to specify the fasta file.\nBye for now!\n"
		sys.exit()
	
        for i in SeqIO.parse(open(global_centroids), 'fasta'):
                global_cluster_count += 1
	
	if not os.path.exists('GLOBAL'):
		os.makedirs('GLOBAL')
	os.chdir('GLOBAL')
	print "creating local copies of relevant files if necessary"
	if not os.path.isfile(args.BIOM_input.split("/")[-1]):
		shutil.copy(args.BIOM_input, './')
	if not os.path.isfile(global_centroids.split("/")[-1]):
		shutil.copy(global_centroids, './')

	
if args.metadata:
	print "\nprovided metadata in file: %s .. " %args.metadata,
	fh = open(args.metadata, "r")
	header_line = fh.readline().strip()
	headers = header_line.split(",")
	for line in fh:
		line = line.strip()
		cols = line.split(",")
		if not len(cols) == len(headers):
			print "\n\tsample %s in metadata file has an invalid number of columns - should have %i / has %i\n" %(cols[0], len(headers), len(cols))
			sys.exit()
		for i in range(1,len(cols)):
			metadata[cols[0]][headers[i]] = cols[i]
	for sID in BIOM_tables_per_method['OTU_denovo'].ids(axis='sample'):
		if not metadata.has_key(sID):
			print "\n\tThe sample %s has no metadata available\n" %sID
			sys.exit()
	fh.close()
	print "ok!\n"
	add_sample_metadata_to_biom(in_table=BIOM_tables_per_method['OTU_denovo'], metadata=metadata, v=args.verbose)

if args.blast or args.blast_xml or args.pplace or args.kraken:
	if args.blast: # or args.blast_xml:
                methods.append('blast')
        if args.pplace:
                methods.append('pplacer')
	if args.kraken:
		methods.append('kraken')

	for approach in methods:
		if approach == 'blast' and ( global_cluster_count > 0): # or args.blast_xml ):
	
			print "\n##### BLAST #####\n"
	
			if not os.path.exists('BLAST_'+str(args.min_ident)):
				os.makedirs('BLAST_'+str(args.min_ident))
			os.chdir('BLAST_'+str(args.min_ident))

			if args.REFlist:
				print "\n### BUILDING BLAST DATABASE ###\n"
				write_out_refs_to_fasta(ref_seqs=all_seqs, ref_taxids=reference_taxa)
				args.blast_db = makeblastdb(in_fasta='refs.fasta', dbtype='nucl', out_prefix=args.marker)
			elif args.blast_db:
				print "\n### USING PRECOMPILED BLAST DATABASE (%s) ###\n" %args.blast_db

			if not args.blast_xml:
				print "\n### RUNNING BLAST ###\n"

				print "running blast search against database %s" % args.blast_db
				blast_out = "global_blastn.out.xml"# % (args.marker, queryID)
				blast_cmd = "blastn -query %s -db %s -evalue 1e-20 -outfmt 5 -out %s -num_threads %i -max_target_seqs 50" % (global_centroids, args.blast_db, blast_out, args.n_threads) 
				print blast_cmd #this is just for the output, the actual blast command is run using the NCBI module below 

				blast_handle = NcbiblastxCommandline(cmd='blastn', query=global_centroids, db=args.blast_db, evalue=1e-20, outfmt=5, out=blast_out, num_threads=args.n_threads, max_target_seqs=50)		
				stdout, stderr = blast_handle()

				args.blast_xml = os.path.abspath('global_blastn.out.xml')
				print '\n'+time.strftime("%c")+'\n'
			else:
				print "\nusing existing BLAST result from '%s'\n" %args.blast_xml
				blast_out = args.blast_xml

			print "\n### INTERPRETING BLAST RESULTS ###\n"
			blast_result_handle = open(blast_out) #read in blast result (in xml format)

			print "\nparse blast xml file\n"
			blast_results = NCBIXML.parse(blast_result_handle) #parse blast xml file

			print "\nfiltering BLAST output\n"
			res_dict = {}
			res_dict = blast_filter(b_result=blast_results, v=args.verbose, m_ident=args.min_ident, m_bitscore=args.min_bit, m_ali_length=args.min_ali_length, bit_score_cutoff=args.bitscore_skim_LCA, bit_score_cutoff_adjust_off=args.bitscore_skim_adjust_off, strand_reversed=queries_to_rc)

			if res_dict['format'] == 'gb':	#This is the case if BLAST was performed against Genbank
				print '\n'+time.strftime("%c")+'\n'
				print "\nfetch taxids for gb accessions\n"
				taxid_list = gb_to_taxid(b_filtered=res_dict, all_taxids=taxid_list, processed=gb_to_taxid_dict, v=args.verbose)

			elif res_dict['format'] == 'unknown': #blast_filter function did not assign a format because no good hit was found to actually assess the format
				if not tax_dict.has_key('tax_id'): #in this rare case I check if the tax_dict is in the proper format and if not
                                       	tax_dict['tax_id'] = [] #just add the tax_id key and an empty list 

			print "\nestablishing taxonomy for reference sequences from custom database\n"
#			if not taxid_list:
#		                if reference_taxa.values():
#					taxid_list = reference_taxa.values()
#				else:
#					taxid_list = extract_taxid_list_from_result(res_dict['hit'])
			taxid_list = extract_taxid_list_from_result(res_dict['hit'])

                	make_tax_dict(tids=taxid_list, out_tax_dict=tax_dict, denovo_taxa=denovo_taxa, ref_taxa=reference_taxa)

			print "\nassign taxonomy\n"
			taxonomy_count = assign_taxonomy_LCA(b_filtered=res_dict, tax_dict=tax_dict, v=args.verbose)

			if not tax_dict.has_key('tax_id'): #this is only relevant in the rare case when no valid sequence made it throuh the trimming/clustering and no reference taxids have been produced earlier from a custom reference
                		tax_dict['tax_id'] = ['nohit']
			
			BIOM_tables_per_method[approach] = annotation_BIOM_table_with_taxonmy(BIOM_table=BIOM_tables_per_method['OTU_denovo'], hits_per_tax_level=taxonomy_count, taxonomy_dictionary=tax_dict, prefix=args.output_prefix, method=approach)	

			print "\n##### BLAST ANALYSIS DONE #####\n"
			print "Find result tables in '%s'\n\n" %os.path.abspath('.')
			os.chdir('../')
                	print '\n'+time.strftime("%c")+'\n'


		elif approach == 'pplacer' and global_cluster_count > 0:
			taxonomy_count = defaultdict(dict)
			pplacer_out_dict = {}
                        print "\n##### PHYLOGENETIC PLACEMENT WITH PPLACER #####\n"
			if not os.path.exists('PPLACER'):
				os.makedirs('PPLACER')
			os.chdir('PPLACER')

			print "\nestablishing taxonomy for reference sequences from custom database\n"
                        taxid_list = reference_taxa.values()
                        make_tax_dict(tids=taxid_list, out_tax_dict=tax_dict, denovo_taxa=denovo_taxa, ref_taxa=reference_taxa)

			if not args.jplace:
				print "\n### PREPARING PPLACER INPUT ###\n" 
			
				print "pplacer analysis will be limited to queries with at least %s %% similarity to the reference database\n" %(args.min_ident*100)
				print "This will be determined based on a blast analysis .."
				if args.blast_xml:
					print "Filtering queries for pplacer based on results from previous BLAST result : %s" %args.blast_xml
				else:
					print "Please provide a blast result in xml format '--blast_xml' or simply do a full blast analysis (because you can) via '--blast'"
					sys.exit()

				blast_result_handle = open(args.blast_xml) #read in blast result (in xml format)

				print "parse blast xml file"
				blast_results = NCBIXML.parse(blast_result_handle) #parse blast xml file
	
				print "pre-bin based on blast output .. ",
				pplacer_out_dict = pre_pplacer_filter(blast_handle=blast_results, m_ident=0.8, m_ali_length=args.min_ali_length, v=args.verbose)
				
				extract_queries_plus_rc(infile=global_centroids, good_list=pplacer_out_dict['to_pplacer'], bad_list=pplacer_out_dict['nohit'], rc_list=pplacer_out_dict['strand_reversed'], out_prefix='pplacer')
				profile,aln_fasta = parse_refpkg_json(args.refpkg)
				print "\n### ALIGN QUERIES TO HMM PROFILE ###\n"
				run_hmmalign(os.path.abspath('pplacer.queries.fasta'), aln_fasta, profile, out='queries_to_profile.sto')
	
				print "\n### RUNNING PPLACER ###\n"
				run_pplacer(args.refpkg, queries_to_profile='queries_to_profile.sto', prefix='pplacer', v=args.verbose)
				args.jplace = os.path.abspath('pplacer.jplace')

			else:
				print "using existing pplacer result: %s" %args.jplace
	
			print "\n### INTERPRETING PPLACER RESULTS ###\n"
			pplacer_out_dict = parse_pplacer(jplace_file=args.jplace, result_dict=pplacer_out_dict)
			taxonomy_count = assign_taxonomy_pplacer(pplacer_out=pplacer_out_dict, tax_dict=tax_dict, v=args.verbose)

			if not tax_dict.has_key('tax_id'): #this is only relevant in the rare case when no valid sequence made it throuh the trimming/clustering and no reference taxids have been produced earlier from a custom reference
                		tax_dict['tax_id'] = ['nohit']
			
			BIOM_tables_per_method[approach] = annotation_BIOM_table_with_taxonmy(BIOM_table=BIOM_tables_per_method['OTU_denovo'], hits_per_tax_level=taxonomy_count, taxonomy_dictionary=tax_dict, prefix=args.output_prefix, method=approach)	

			print BIOM_tables_per_method

			print "\n##### PPLACER ANALYSIS DONE #####\n"
			print "Find result tables in '%s'\n\n" %os.path.abspath('.')
			os.chdir('../')
                	print '\n'+time.strftime("%c")+'\n'

		elif approach == 'kraken' and global_cluster_count > 0:
			taxonomy_count = defaultdict(dict)
			kraken_out_dict = {}
                        print "\n##### RUNNING KRAKEN #####\n"
		
			if not os.path.exists('KRAKEN'):
				os.makedirs('KRAKEN')
			os.chdir('KRAKEN')
			
			if not args.kraken_db:
				# Format Sequences for input to kraken database build
				gb_to_kraken_db(all_seqs, args.marker)
				#Build database
				full_build_custom_kraken_db(args.marker+'.kraken.fasta', db_path='.', db_name='KRAKEN_DB', threads=args.n_threads, jellyfish_hash_size=args.jellyfish_hash_size, v=args.verbose)
				args.kraken_db = os.path.abspath('./KRAKEN_DB')
				print "\nSuccessuflly generated at: %s\n" %args.kraken_db
			else:
				print "\nUsing precompiled kraken database: %s\n" %args.kraken_db

			kraken_out_file = run_kraken(kraken_db=args.kraken_db, queries_fasta=global_centroids, threads=args.n_threads)
			kraken_out_file = kraken_filter(kraken_db=args.kraken_db, kraken_out=kraken_out_file, score_threshold=args.kraken_score_threshold)

			kraken_out_dict = parse_kraken_out(kraken_out=kraken_out_file)
			
			print "\nestablishing taxonomy for reference sequences from custom database\n"
#	                if reference_taxa.values():
#				taxid_list = reference_taxa.values()
#			else:
#				taxid_list = extract_taxid_list_from_result(kraken_out_dict['hit'])
#			
			taxid_list = extract_taxid_list_from_result(kraken_out_dict['hit'])

                	make_tax_dict(tids=taxid_list, out_tax_dict=tax_dict, denovo_taxa=denovo_taxa, ref_taxa=reference_taxa)

			taxonomy_count = assign_taxonomy_kraken(kraken_out=kraken_out_dict, tax_dict=tax_dict, v=args.verbose)

			if not tax_dict.has_key('tax_id'): #this is only relevant in the rare case when no valid sequence made it throuh the trimming/clustering and no reference taxids have been produced earlier from a custom reference
                		tax_dict['tax_id'] = ['nohit']
			
			BIOM_tables_per_method[approach] = annotation_BIOM_table_with_taxonmy(BIOM_table=BIOM_tables_per_method['OTU_denovo'], hits_per_tax_level=taxonomy_count, taxonomy_dictionary=tax_dict, prefix=args.output_prefix, method=approach)	

			if args.rm_kraken_db:
				print "\nRemoving Kraken database\n"
				shutil.rmtree(args.kraken_db)

			print "\n##### KRAKEN ANALYSIS DONE #####\n"
			print "Find result tables in '%s'\n\n" %os.path.abspath('.')
			os.chdir('../')
                	print '\n'+time.strftime("%c")+'\n'


else:
	print "\nPlease select a method for taxonomic assignment:\n--blast\n--pplace\n--kraken\n"

os.chdir('../')
os.utime('GLOBAL', None)


print "\n##### DONE! #####\n"
print '\n'+time.strftime("%c")+'\n'

sys.exit()

#
#add functionality to iterate over parameters, e.g. cluster_coverage 10,20,30,40,50 and treat each as different sample in the biom file

