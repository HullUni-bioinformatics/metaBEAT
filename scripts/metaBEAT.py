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



##############set this, or put a file called taxonomy.db in the same directory as the metaBEAT.py script############
taxonomy_db = '/home/chrishah/src/taxtastic/taxonomy_db/taxonomy.db'
	
#############################################################################
VERSION="0.8"
informats = {'gb': 'gb', 'genbank': 'gb', 'fasta': 'fasta', 'fa': 'fasta', 'fastq': 'fastq'}
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
gi_to_taxid_dict = {}
taxid_list = []
tax_dict = {}

parser = argparse.ArgumentParser(description='metaBEAT - metaBarcoding and Environmental DNA Analyses tool', prog='metaBEAT.py')
#usage = "%prog [options] REFlist"
#parser = argparse.ArgumentParser(usage='%(prog)s [options] REFlist', formatter_class=RawTextHelpFormatter)

parser.add_argument("-Q", "--querylist", help="file containing a list of query files", metavar="<FILE>", action="store")
parser.add_argument("-v", "--verbose", help="turn verbose output on", action="store_true")
parser.add_argument("-s", "--seqinfo", help="write out seq_info.csv file", action="store_true")
parser.add_argument("-f", "--fasta", help="write out ref.fasta file", action="store_true")
parser.add_argument("-p", "--phyloplace", help="perform phylogenetic placement", action="store_true")
#parser.add_argument("-t", "--taxids", help="write out taxid.txt file", action="store_true")
parser.add_argument("-b", "--blast", help="compile local blast db and blast queries", action="store_true")
parser.add_argument("-m", "--marker", help="marker ID (default: marker)", metavar="<string>", action="store", default="marker")
parser.add_argument("-n", "--n_threads", help="Number of threads (default: 1)", type=int, metavar="<INT>", action="store", default="1")
parser.add_argument("-E", "--extract_centroid_reads", help="extract centroid reads to files", action="store_true")
parser.add_argument("-e", "--extract_all_reads", help="extract reads to files", action="store_true")
query_group = parser.add_argument_group('Query preprocessing', 'The parameters in this group affect how the query sequences are processed')
query_group.add_argument("--PCR_primer", help='PCR primers (provided in fasta file) to be clipped from reads', metavar="<FILE>", action="store")
query_group.add_argument("--trim_adapter", help="trim adapters provided in file", metavar="<FILE>", action="store")
query_group.add_argument("--trim_qual", help="minimum phred quality score (default: 30)", metavar="<INT>", type=int, action="store", default=30)
query_group.add_argument("--phred", help="phred quality score offset - 33 or 64 (default: 33)", metavar="<INT>", type=int, action="store", default=33)
query_group.add_argument("--trim_window", help="sliding window size (default: 5) for trimming; if average quality drops below the specified minimum quality all subsequent bases are removed from the reads", metavar="<INT>", type=int, action="store", default=5)
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
reference_group.add_argument("--blast_db", help="path to precompiled blast database", metavar="<PATH>", action="store", default="")
reference_group.add_argument("--gi_to_taxid", help="comma delimited file containing 'gi accession,taxid' for a list of taxa", metavar="<FILE>", action="store", default=os.getcwd()+"/gi_to_taxid.csv")
cluster_group = parser.add_argument_group('Query clustering options', 'The parameters in this group affect read clustering')
cluster_group.add_argument("--cluster", help="perform clustering of query sequences using vsearch", action="store_true")
cluster_group.add_argument("--clust_match", help="identity threshold for clustering in percent (default: 1)", type=float, metavar="<FLOAT>", action="store", default="1")
cluster_group.add_argument("--clust_cov", help="minimum number of records in cluster (default: 1)", type=int, metavar="<INT>", action="store", default="1")
blast_group = parser.add_argument_group('BLAST search', 'The parameters in this group affect BLAST search and BLAST based taxonomic assignment')
blast_group.add_argument("--www", help="perform online BLAST search against nt database", action="store_true")
blast_group.add_argument("--min_ident", help="minimum identity threshold in percent (default: 0.80)", type=float, metavar="<FLOAT>", action="store", default="0.80")
blast_group.add_argument("--min_ali_length", help="minimum alignment length in percent of total query length (default: 0.95)", type=float, metavar="<FLOAT>", action="store", default="0.95")
blast_group.add_argument("--min_bit", help="minimum bitscore (default: 80)", type=int, metavar="<INT>", action="store", default="80")
phyloplace_group = parser.add_argument_group('Phylogenetic placement', 'The parameters in this group affect phylogenetic placement')
phyloplace_group.add_argument("--refpkg", help="PATH to refpkg", metavar="<DIR>", action="store")

biom_group = parser.add_argument_group('BIOM OUTPUT','The arguments in this groups affect the output in BIOM format')
biom_group.add_argument("-o","--output_prefix", help="prefix for BIOM output files (default='metaBEAT')", action="store", default="metaBEAT")
biom_group.add_argument("--metadata", help="comma delimited file containing metadata (optional)", action="store")
biom_group.add_argument("--mock_meta_data", help="add mock metadata to the samples in the BIOM output", action="store_true")

Entrez_group = parser.add_argument_group('Entrez identification','metaBEAT is querying the NCBI Entrez databases, please provide an email address for identification')
Entrez_group.add_argument("-@", "--email", help='provide your email address for identification to NCBI', metavar='<email-address>', action="store", default="")

parser.add_argument("--version", action="version", version='%(prog)s v.'+VERSION)
args = parser.parse_args()


if len(sys.argv) < 2:	#if the script is called without any arguments display the usage
    parser.print_usage()
    sys.exit(1)


###FUNCTIONS

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
			print "\nnot sure %s is an email address\n" %mail
			sys.exit()
		print "You have specified: '%s'\n" %(mail)

		FH = open(os.path.dirname(sys.argv[0])+'/user_email.txt','w')
		FH.write(mail)
		FH.close()

	return mail
	

def rw_gi_to_taxid_dict(dictionary, name, mode):
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
        sys.exit("only 'r' or 'w' are allowed for mode in the rw_gi_to_taxid_dict")
        


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
	    v_string = "excluding %s from further analyses - cluster_size filter" %record.id
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



def parse_vsearch_uc(fil, cluster_counts, cluster_reads, extract_reads=1):
    "The function parses the 'uc' file from vsearch"
    f=open(fil,"r")
    for line in [l.strip() for l in f]: #loop through the file one line at a time, stripping off any newline characters
        if line.startswith("C"): #process only lines that start with a "C"
            elem = line.split("\t") #split the lines at tab
            cluster_counts[elem[8]] = int(elem[2]) #write the counts to dictionary with key being the id of the centroid read
        
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


def assign_taxonomy_LCA(b_filtered, tax_dict, v=0):
    "The function takes a dictionary of queries and their hits"
    "and provides taxonomic assignment based on the LCA method."
    tax_count = defaultdict(dict)

#    print len(b_filtered['hit'])
    if b_filtered.has_key('hit'):
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
#	                    print tax_dict[tax][index]
	                    id_list.append(tax_dict[tax][index])
	                    if not tax_dict[tax][index]:
#	                        print "nothing found at this level for %s" %tax
	                        break
	                if len(id_list) == len(b_filtered['hit'][query]):
#	                        print "ok - all have valid taxid"

	                        if len(set(id_list)) == 1:
				    if v:
		                            print "found LCA %s at level %s" %(id_list[0], tax_dict["tax_id"][index])
#	                            print tax_dict[id_list[0]][1]
#	                            print tax_dict[id_list[0]][2]
	                            if not tax_count[tax_dict["tax_id"][index]].has_key(id_list[0]):
	                                tax_count[tax_dict["tax_id"][index]][id_list[0]] = []
 	                            tax_count[tax_dict["tax_id"][index]][id_list[0]].append(query)
	                            del b_filtered['hit'][query]
	                            break
#	                        else:
#	                            print "not yet LCA"
                            
#	        print "\nUPDATE:\n%s" %tax_count
    	    if len(b_filtered['hit']) == 0:
	        print "all queries have been successfully assigned to a taxonomy"
	    else:
	        print "%i queries failed:\n%s" %(len(b_filtered['hit']), b_filtered['hit'])

    if b_filtered.has_key('nohit'):
        if not tax_count.has_key('nohit'):
            tax_count['nohit'] = {'nohit':[]}
        tax_count['nohit']['nohit'].extend(b_filtered['nohit'])
    
        
    return tax_count



def make_tax_dict(tids, out_tax_dict, denovo_taxa, ref_taxa):
    "The function takes a list of taxonomic ids, runs the taxtable program from the taxit suite"
    "to summarize the taxonomic information for the taxids and formats the result into a dictionary"

    print "write out taxids to taxids.txt\n"
    f = open("taxids.txt","w")
    for tid in tids:
        f.write(tid + "\n")
    f.close()

    cmd = "taxit taxtable -d %s -t taxids.txt -o taxa.csv" %taxonomy_db
    print "running taxit to generate reduced taxonomy table"
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

#	while array[-1] != 'species': #this should get rid of any 'subspecies' or 'varietas' levels
#	    array.pop(-1)
		
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

	while out_tax_dict['tax_id'][-1] != 'species': #this should get rid of any 'subspecies' or 'varietas' levels
            out_tax_dict['tax_id'].pop(-1)



def check_for_taxonomy_db(tax_db):
    "The function checks if the taxonomy.db is present"
    if not tax_db or not os.path.isfile(tax_db): #if no path to the taxonomy database is specifed or the path is not correct
                if os.path.isfile(os.path.dirname(sys.argv[0])+'/taxonomy.db'): #check if taxonomy.db is present in the same path as the metabeat.py script
                        tax_db = "%s/taxonomy.db" %os.path.dirname(sys.argv[0])
			print "taxonomy.db found at %s" %tax_db
                else:
                        print "\nmetaBEAT.py requires a taxonomy database. Please configure it using taxtastic. Per default metaBEAT expects a file 'taxonomy.db' in the same directory that contains the metaBEAT.py script. If you are not happy with that please change the variable 'taxonomy_db' at the top of this script to the correct path to your taxonomy.db file.\n"
                        sys.exit()
    else:
        print "taxonomy.db is present at: %s" %tax_db

    return tax_db

def gi_to_taxid(b_filtered, all_taxids, processed, v=0):
    "This function takes the resulting dictionary from the blast_filter function"
    "and fetches the corresponding taxonomic ids to the gi accessions if necessary"
    "(i.e. if the format is specified as 'gi')."
    "It returns a list of unique taxonomic ids."
    done = 0
    if b_filtered['format'] == 'gi':
        for query in b_filtered['hit'].keys():
            if v:
                print "processing query: %s" %query
            taxids=[]
            for hit in b_filtered['hit'][query]:
#                print "gi: %s" %hit
                if processed.has_key(hit):
                    taxid = processed[hit]
                    done += 1
                    if v:
                        print "have seen gi %s (%i) before" %(hit, done)
                else:
                    
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
                    done += 1
                    if v:
                        print "fetched taxid for gi: %s (%i) using Entrez" %(hit, done)
                        
#                print "taxid: %s" %taxid                
                taxids.append(taxid)
                if not v:
                    if done%100 == 0:
                        print ".", #the comma at the end ommits the newline when printing
                
#            print b_filtered['hit'][query]
            b_filtered['hit'][query] = list(set(taxids))
#            print len(taxids)
            all_taxids.extend(taxids)
#            print len(all_taxids)
            
        if v:
            print "\nfetched taxids for %i gis" %done
        
    elif b_filtered['format'] == 'taxid':
        for query in b_filtered['hit'].keys():
            all_taxids.extend(b_filtered['hit'][query])

#    print all_taxids
#    print "make to set"
    all_taxids = list(set(all_taxids))
    return all_taxids


def blast_filter(b_result, v=0, m_bitscore=80, m_ident=0.8, m_ali_length=0.95):
    "The function interprets a BLAST results handle and filters results subsequent taxonomic assignment using LCA"
    result = {'format':''}
    count=0
    for res in b_result:
	bit_score_cutoff = 0.9 #that's the default setting, i.e. top 10% of top bit score
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

	elif (len(res.alignments[0].hsps[0].query) < res.query_length*m_ali_length):
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
            
        else: #if a hit has been found
	    if not result.has_key('hit'):
		result['hit'] = {}

            if (float(res.alignments[0].hsps[0].identities)/len(res.alignments[0].hsps[0].query) == 1): #if we have a full length 100 % match adjust the bitscore so window so that only this hit is considered
		bit_score_cutoff = 1
		if v:
			print "\n100% match:\n"
			print "Query length: %i\nAlignment length: %s\nnumber of identitites: %s" %(res.query_length, str(len(res.alignments[0].hsps[0].query)), str(res.alignments[0].hsps[0].identities))
            max_bit_score = res.alignments[0].hsps[0].bits #record the maximum bitscore

            result['hit'][res.query]=[] #create empty list for query
            for alignment in res.alignments: #for each hit of the current query
#                print alignment.hsps[0].bits
                if not result['format']: #determining the format of the blast database (only once per blast file)
                    if alignment.title.startswith('gi'):
                        result['format']='gi'
                    else:
                        result['format']='taxid'
                        
                if alignment.hsps[0].bits >= (max_bit_score*bit_score_cutoff): #if a hit has a bitscore that falls within the top 90 % of the bitscores recorded
#                    print alignment.title.split("|")[1]
		    if result['format'] == 'gi':
                        result['hit'][res.query].append(alignment.title.split("|")[1])
		    elif result['format'] == 'taxid':
                        result['hit'][res.query].append(alignment.title.split("|")[-2])
			
                    if v:
                        print "%s\t%s" %(alignment.hsps[0].bits, alignment.title)
                else:
#                    print result['hit'][res.query]
                    break
#            print "%s: %s" %(res.query, result['hit'][res.query])

    print "%i queries processed" %count
    
    if not result['format']:	#if no database format could be determined at this point that means there was no significant hit, so we set the format to unknown
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


###########


####START OF MAIN PROGRAM

print '\n'+time.strftime("%c")+'\n'
print "%s\n" % (' '.join(sys.argv))

if args.email:
	Entrez.email = args.email
Entrez.email = check_email(mail = Entrez.email)#mail=Entrez.email)	

if args.blast or args.phyloplace:
	taxonomy_db = check_for_taxonomy_db(tax_db=taxonomy_db)

if args.phyloplace:
	if not args.refpkg:
		print "\nTo perform phylogenetic placement with pplacer metaBEAT currently expects a reference package to be specified via the --refpkg flag\n"
		parser.print_help()
		sys.exit()
	
	if not os.path.isdir(args.refpkg):
		print "\nThe specified reference package does not seem to be valid\n"
		sys.exit()
	if not args.blast:
		print "\nPhylogenetic placement currently requires a blast search first to determine the set of queries for which phylogenetic placement is attempted - add the --blast flag to your command\n"
		sys.exit()	
	args.refpkg = os.path.abspath(args.refpkg) 

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
elif args.blast_db:
	for f in glob.glob(args.blast_db+"*"):
		if not bl_db_extensions:
			args.blast_db = os.path.abspath(args.blast_db)
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
#				print "processing sample: %s" %querydata[0]
				if not informats.has_key(querydata[1]):
					print "query file format for sample %s is invalid" % querydata[0]
					per_sample_query_format=querydata[1]
					sys.exit(0) 
				for i in range(2,len(querydata)):
					if re.match("^[ACTG]*$", querydata[i]):
#						print "column %i looks like a barcode" %i
						per_sample_barcodes.append(querydata[i])
					else:
						if not os.path.isfile(querydata[i]):
							print "%s is not a valid file" %querydata[i]
							sys.exit(0)
						per_sample_query_files.append(os.path.abspath(querydata[i]))
				
				queries[querydata[0]]['format'] = querydata[1]
				queries[querydata[0]]['files'] = per_sample_query_files
				if per_sample_barcodes:
#					print "barcodes %s found for sample %s" %("-".join(per_sample_barcodes), querydata[0])
					queries[querydata[0]]['barcodes'] = per_sample_barcodes
					files_to_barcodes["|".join(per_sample_query_files)]["-".join(per_sample_barcodes)] = querydata[0]
#			print queries[querydata[0]]
#			print files_to_barcodes
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

	if args.metadata:
		fh = open(args.metadata, "r")
		header_line = fh.readline().strip()
		headers = header_line.split(",")
		for line in fh:
			line = line.strip()
			cols = line.split(",")
			if not len(cols) == len(headers):
				print "sample %s in metadata file has an invalid number of columns - should have %i / has %i\n" %(cols[0], len(headers), len(cols))
				sys.exit()
			for i in range(1,len(cols)):
#				print cols[i]
				metadata[cols[0]][headers[i]] = cols[i]
#		print metadata
		for sID in queries.keys():
			if not metadata.has_key(sID):
				print "The sample %s has no metadata available\n" %sID
#			else:
#				print metadata[sID]
		fh.close()

#print len(queries)
#print len(files_to_barcodes[files_to_barcodes.keys()[0]])
#print files_to_barcodes

#check if gi_to_taxid dictonary is there and read if yes
#print "looking for gi_to_taxid file at %s" %args.gi_to_taxid
if os.path.isfile(args.gi_to_taxid):
	print "\nfound gi_to_taxid file at %s" %args.gi_to_taxid,
	rw_gi_to_taxid_dict(dictionary=gi_to_taxid_dict, name=args.gi_to_taxid, mode='r')
	

#print references
print "\n######## PROCESSING REFERENCE DATA ########\n"

for reffile, refformat in references.items():
	seqs=list(SeqIO.parse(reffile,refformat, alphabet=generic_dna))
	print "\nprocessing %s (containing %i records)" % (reffile,len(seqs))
	for i in reversed(range(len(seqs))):
#		print seqs[i].description
		taxid = ''
		if not seqs[i].name in records:	#records list keeps track of all records.
			records[seqs[i].name] = 1
		else:				#if a record id occurs again, say in the case of custom fasta sequences the record id is the first name in the fasta header, e.g. the Genus name
			records[seqs[i].name] += 1
			seqs[i].name += "_%s" % records[seqs[i].name]

		if not seqs[i].features:		#if the record does not have a source feature yet. Will happen for records originating from custom fasta files
#			print "the record has no features\n"	
			seqs[i].features.append(SeqFeature(FeatureLocation(0,len(seqs[i].seq),strand=1), type="source"))	#create source feature
			seqs[i].features[0].qualifiers['original_desc'] = seqs[i].description
			
		if not seqs[i].features[0].qualifiers.has_key('organism'):	#if the record does not have an organism key. Again will happen for custom fasta sequences
#			print seqs[i].description
			seqs[i].features[0].qualifiers['organism'] = ["%s %s" % (seqs[i].description.split(" ")[0], seqs[i].description.split(" ")[1])] #add an organism key and make it the first 2 words in the record description will are per default the first 2 words in the fasta header
#			seqs[i].features[0].qualifiers['organism'] = ["%s" % seqs[i].description.split(" ")[0]]	#add an organism key and make it the first word in the record description will are per default the first 2 words in the fasta header
#			sp_search = re.compile('^sp$|^sp[.]$|^sp[0-9]+$')
#			if not sp_search.match(seqs[i].description.split(" ")[1]):
#				seqs[i].features[0].qualifiers['organism'] = ["%s %s" % (seqs[i].features[0].qualifiers['organism'][0], seqs[i].description.split(" ")[1])]	#if the second word in the description is not 'sp' or 'sp.', i.e. it is unknown, then add also the second word to the organism field

		if not seqs[i].annotations.has_key('date'):	#set the date in the record if it does not exist
			seqs[i].annotations['date'] = date
			
		if not seqs[i].features[0].qualifiers['organism'][0] in reference_taxa:	#if the organism name has not yet been encountert the following lines will find the corresponding taxid, add it to the record and store it in the reference_taxa dictionary
#			print "seeing %s for the first time" %seqs[i].features[0].qualifiers['organism'][0]
			if seqs[i].features[0].qualifiers.has_key('db_xref'):	#there may be more than one 'db_xref' qualifiers, e.g. sometimes there are BOLD ids also in 'db_xref'
				for j in range(len(seqs[i].features[0].qualifiers['db_xref'])):
					if seqs[i].features[0].qualifiers['db_xref'][j].startswith('taxon'):
						taxid = seqs[i].features[0].qualifiers['db_xref'][j].split(":")[-1]
#						print "found valid taxid in record: %s" %taxid
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
#						taxids[taxon['IdList'][0]]+= 1
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
#					taxids[taxid] += 1
#					print "adding %s with taxid: %s to the dictionary" %(seqs[i].features[0].qualifiers['organism'][0], taxid)
					reference_taxa[seqs[i].features[0].qualifiers['organism'][0]] = taxid
			else:
				seqs[i].id = seqs[i].name
				reference_taxa[seqs[i].features[0].qualifiers['organism'][0]] = taxid
#				taxids[taxid] += 1
				
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

#			print "Have seen %s before. The taxid is: %s" %(seqs[i].features[0].qualifiers['organism'][0], reference_taxa[seqs[i].features[0].qualifiers['organism'][0]])

		
		if len(taxid) > 0:
#			taxids[taxid] += 1
			current = '"%s","%s","%s","%s",0' % (seqs[i].name, seqs[i].name, taxid, seqs[i].features[0].qualifiers['organism'][0]) #producing a string for the current record for the seq_info file that is needed for the reference package required by pplacer
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
#		gb_out.write(seqs.format("genbank"))
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

if args.fasta:	#this bit writes out the sequences that will become the reference dataset
	write_out_refs_to_fasta(ref_seqs=all_seqs, ref_taxids=reference_taxa)

if args.blast:
	if args.REFlist:
		print "\nestablishing taxonomy for reference sequences from custom database\n"
		taxid_list = reference_taxa.values()
#		taxid_list = [i for i in reference_taxa.values() if not 'denovo' in i] # reference_taxa.values()
		make_tax_dict(tids=taxid_list, out_tax_dict=tax_dict, denovo_taxa=denovo_taxa, ref_taxa=reference_taxa)

		print "\n### BUILDING BLAST DATABASE ###\n"
		if not args.fasta:	#if the reference sequences have not yet been written out as fasta it is done here
			write_out_refs_to_fasta(ref_seqs=all_seqs, ref_taxids=reference_taxa)
#		print "building blast db for marker %s\n" % args.marker
		makeblastdb="makeblastdb -in refs.fasta -dbtype nucl -out %s_blast_db" % args.marker
		print makeblastdb
		cmdlist = shlex.split(makeblastdb)
		cmd = subprocess.Popen(cmdlist, stdout=subprocess.PIPE)
		stdout = cmd.communicate()[0]
		if args.verbose:
			print stdout
	
		blast_db = os.path.abspath('.')+"/"+"%s_blast_db" % args.marker
	elif args.blast_db:
		print "\n### USING PRECOMPILED BLAST DATABASE (%s) ###\n" %args.blast_db
		blast_db = args.blast_db


#start read stats output
read_counts_out = open(args.output_prefix+"_read_stats.csv","w")
outstring = "sample,"+",".join(read_metrics)
read_counts_out.write(outstring+"\n")
read_counts_out.close()


print '\n'+time.strftime("%c")+'\n'
querycount = defaultdict(int)
if args.blast or args.phyloplace or args.merge or args.cluster:
	##determine which methods will be applied
	if args.blast:
		methods.append('blast')
	if args.phyloplace:
		methods.append('pplacer')

	print "\n######## PROCESSING QUERIES ########\n"
	if files_to_barcodes:
		print "\n### DEMULTIPLEXING ###\n"
#		print "Barcodes detected - initialize demultiplexing"
		if not os.path.exists('demultiplexed'):
			os.makedirs('demultiplexed')
			
		os.chdir('demultiplexed')
		for lib in files_to_barcodes.keys():
			data = lib.split("|")
			compr_mode = ""
			print "assessing basic characteristics"
#			print files_to_barcodes[lib]
			barcodes_num = len(files_to_barcodes[lib].keys()[0].split("-"))
			if barcodes_num == 2:
#				print "interpreting barcodes as forward and reverse inline"
				barcode_mode = "--inline_inline"
			elif barcodes_num == 1:
#				print "interpreting barcodes as forward inline only"
				barcode_mode = "--inline_null"
				
			if lib.endswith(".gz"):
				compr_mode = "gzfastq"
#				print "set mode to %s" %compr_mode
			else:
				compr_mode = "fastq"
#				print "set mode to %s" %compr_mode

#			print "create temporal barcode file"
			BC = open("barcodes","w")
			for sample in files_to_barcodes[lib].keys():
#				print "This is sample: %s" %sample
				bc = sample.split("-")
#				print bc
				if len(bc) != barcodes_num:
					print "expecting %i barcode(s) for all samples - something's wrong with sample %s" %(barcodes_num, sample)
					sys.exit()
				BC.write("\t".join(bc)+"\n")

			BC.close()
			if len(data)>1:
				print "data comes as paired end - ok"
				demultiplex_cmd="process_shortreads -1 %s -2 %s -i %s -o . -b barcodes %s -y fastq" %(data[0], data[1], compr_mode, barcode_mode)
				print demultiplex_cmd
				cmdlist = shlex.split(demultiplex_cmd)
				cmd = subprocess.Popen(cmdlist, stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
				stdout,stderr = cmd.communicate() #[0]
				print stdout
				print stderr

#				print "renaming results"
				for sample in files_to_barcodes[lib].keys():
#					print files_to_barcodes[lib][sample]
					new_file_list=[]
					for a in range(1,3):
						old = "sample_%s.%i.fq" %(sample,a)
#						print old
						new = "%s_%i.fastq" %(files_to_barcodes[lib][sample],a)
#						print new
						os.rename(old,new)
						os.remove("sample_%s.rem.%i.fq" %(sample,a))
						new_file_list.append(os.path.abspath('.')+"/"+new)
						
					queries[files_to_barcodes[lib][sample]]['files'] = new_file_list
			elif len(data)==1:
				print "data comes in single end - ok"
				print "currently not supported"
				sys.exit()

		os.chdir('../')
#	print queries
		

	for queryID in sorted(queries):
#		print queries
#		query_count += 1
#	for queryID, querydata in sorted(queries.items()): #loop through the query files and read in the current query id and the path to the files
		print "\n##### processing query ID: %s #####\n" % (queryID)
		species_count = defaultdict(list)
		taxonomy_count = defaultdict(dict)
		if not os.path.exists(queryID):
			os.makedirs(queryID)
		
#		print queries[queryID]['files']
		
		#determine total read number per sample
		read_stats[queryID]['total'] = int(0)
		for f in queries[queryID]['files']:
			if f.endswith('gz'):
				t_seqs = list(SeqIO.parse(gzip.open(f, 'rb'),queries[queryID]['format']))
			else:
				t_seqs = list(SeqIO.parse(open(f, 'r'),queries[queryID]['format']))
			read_stats[queryID]['total'] += len(t_seqs)
		
#		print "Total number of reads for sample %s: %i" %(queryID, read_stats[queryID][0])

		os.chdir(queryID)

		if (queries[queryID]['format']=="fastq"):
			print "\n### READ QUALITY TRIMMING ###\n"
			if len(queries[queryID]['files'])==2:
				trimmomatic_path="java -jar /usr/bin/trimmomatic-0.32.jar"
				print "\ntrimming PE reads with trimmomatic"
				trimmomatic_exec=trimmomatic_path+" PE -threads %i -phred%i -trimlog trimmomatic.log %s %s %s_forward.paired.fastq.gz %s_forward.singletons.fastq.gz %s_reverse.paired.fastq.gz %s_reverse.singletons.fastq.gz ILLUMINACLIP:%s:5:5:5 TRAILING:%i LEADING:%i SLIDINGWINDOW:%i:%i MINLEN:%i" % (args.n_threads, args.phred, queries[queryID]['files'][0], queries[queryID]['files'][1], queryID, queryID, queryID, queryID, args.trim_adapter, args.trim_qual, args.trim_qual, args.trim_window, args.trim_qual, args.trim_minlength)
				trimmed_files=[queryID+'_forward.paired.fastq.gz', queryID+'_forward.singletons.fastq.gz', queryID+'_reverse.paired.fastq.gz', queryID+'_reverse.singletons.fastq.gz']
			elif len(queries[queryID]['files'])==1:
				print "\ntrimming SE reads with trimmomatic"
				trimmomatic_exec= trimmomatic_path+" SE -threads %i -phred%i %s %s_trimmed.fastq.gz ILLUMINACLIP:%s:5:5:5 TRAILING:%i LEADING:%i SLIDINGWINDOW:%i:%i MINLEN:%i" % (args.n_threads, args.phred, queries[queryID]['files'][0], queryID, args.trim_adapter, args.trim_qual, args.trim_qual, args.trim_window, args.trim_qual, args.trim_minlength)
				trimmed_files=[queryID+'_trimmed.fastq.gz']

			trimmomatic="%s" % trimmomatic_exec
			cmdlist = shlex.split(trimmomatic)
			cmd = subprocess.Popen(cmdlist, stdout=subprocess.PIPE, stderr=subprocess.PIPE)	
			stdout,stderr = cmd.communicate() #[0]
			print trimmomatic
			print stdout
			print stderr
			
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

#				for f in queries[queryID]['files']:
#					cmd = "cat %s | fastx_reverse_complement -Q %i |fastx_clipper -Q %i -n %s | fastx_reverse_complement -Q %i > %s" %(f, args.phred, args.phred, primer_clip_string,  args.phred, f+'.no_primer') 
#					print cmd
#					cmdlist = shlex.split(cmd)
#					cmd = subprocess.call(cmd, shell=True)

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
				
#				print read_stats

				print "\n### Fastq to fasta (plus filtering for length, merged-only, etc. if applicable) ###\n"

				if args.merged_only:
					print "\nkeeping only merged reads for subsequent analyses"
					keep_only(inlist=trimmed_files, pattern_list=['.extendedFrags.fastq.gz'])

				if args.forward_only:
					print "\nkeeping only forward sequences (non merged forward and merged reads) for subsequent analyses"
					keep_only(inlist=trimmed_files, pattern_list=['.extendedFrags.fastq.gz', '_forward', 'notCombined_1'])
			
				survived = concat_and_filter_by_length(inlist=trimmed_files, outfile=queryID+'_trimmed.fasta', excludefile=queryID+'_excluded.fasta', form='fasta', length=args.length_filter, devi=args.length_deviation)
				print "\n%i sequences survived filtering\n" %survived
#				sys.exit()


#				files = " ".join(trimmed_files[-2:])
#				cmd="zcat %s | fastx_reverse_complement -Q %i| fastq_to_fasta -Q %i > temp2.fasta" % (files, args.phred, args.phred)
##				print cmd
#				cmdlist = shlex.split(cmd)
#				cmd = subprocess.call(cmd, shell=True)
#				
#				files = " ".join(trimmed_files[:-2])
#				cmd="zcat %s | fastq_to_fasta -Q %i > temp1.fasta" % (files, args.phred)
##				print cmd
#				cmdlist = shlex.split(cmd)
#				cmd = subprocess.call(cmd, shell=True)
#
#				cmd="cat temp1.fasta temp2.fasta > temp_trimmed.fasta"
##				print cmd
#				cmdlist = shlex.split(cmd)
#				cmd = subprocess.call(cmd, shell=True)
#				
#				os.remove("temp1.fasta")
#				os.remove("temp2.fasta")
			else:
				survived = concat_and_filter_by_length(inlist=trimmed_files, outfile=queryID+'_trimmed.fasta', excludefile=False, form='fasta', length=0)
#				cmd="zcat %s | fastq_to_fasta > temp_trimmed.fasta" % trimmed_files[0]
#				print cmd
#				cmdlist = shlex.split(cmd)
#				cmd = subprocess.call(cmd, shell=True)
				
				#determine the number of reads surviving trimming
#				t_seqs = list(SeqIO.parse('temp_trimmed.fasta','fasta'))
#				read_stats['description'].append('trimmed-total')
#				read_stats[queryID]['trimmed-total'] = survived

#			fin=open("temp_trimmed.fasta","r")
#			f_out=open(queryID+'_trimmed.fasta',"w")
#			for line in fin:
#				f_out.write(line.replace(" ","_"))
				
#			fin.close()
#			f_out.close()

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

		else:
			queryfile = queries[queryID]['files'][0]
			print "\nWARNING - EXPECTING input data to come in SINGLE FASTA FILE\n"


		#getting number of sequences to be queried		
		total_queries = len(list(SeqIO.parse(queryfile,'fasta')))

		if args.cluster:
			##running clustering
			print "\n### CLUSTERING ###\n"
			print "\nclustering using vsearch"
			cmd = "vsearch --cluster_fast %s --id %.2f --strand both --threads %s --centroids %s_centroids.fasta --uc %s.uc" % (queryfile, args.clust_match, args.n_threads, queryID, queryID )
			print cmd
			cmdlist = shlex.split(cmd)
		        stdout, stderr = subprocess.Popen(cmdlist, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate() # , stdout=subprocess.PIPE).communicate()
			if stdout:
				print stdout
	
			#read in the results from the clustering, i.e. number of reads per cluster, later I could read in the actual read ids to allow for retrievel of reads assigned to a given taxon
#			all_clust_count = int(0)
			
			print "\nparse vsearch uc file\n"
			cluster_counts = {}
			cluster_reads = defaultdict(list)
			parse_vsearch_uc(fil=queryID+".uc", cluster_counts=cluster_counts, extract_reads=args.extract_all_reads, cluster_reads=cluster_reads)
			
#			total_queries = querycount[queryID]
			total_clusters = len(cluster_counts)
			
			if args.clust_cov>1:
				print "\nreduce cluster files\n"
			querycount[queryID] = filter_centroid_fasta(centroid_fasta=queryID+'_centroids.fasta', m_cluster_size=args.clust_cov, cluster_counts=cluster_counts, sampleID=queryID, v=args.verbose)		
#			for ID in unknown_seqs_dict.keys():
#				if cluster_counts.has_key(ID):
#					unknown_seqs_dict[ID].description = "%s|%s|%s|%.2f" %(queryID, unknown_seqs_dict[ID].id, cluster_counts[unknown_seqs_dict[ID].id], float(cluster_counts[unknown_seqs_dict[ID].id])/querycount[queryID]*100)
##				print unknown_seqs_dict[ID].description

			print "vsearch processed %i sequences and identified %i clusters (clustering threshold %.2f) - %i clusters (representing %i sequences) passed the filter criteria (minimum of %i sequences per cluster) and will be used in subsequent analyses\n" % (total_queries, total_clusters, float(args.clust_match), len(cluster_counts), querycount[queryID], args.clust_cov)
		
			read_stats[queryID]['clusters_total'] = total_clusters
			read_stats[queryID]['cluster_thres'] = args.clust_match
			read_stats[queryID]['clusters_min_cov'] = args.clust_cov
			read_stats[queryID]['cluster_above_thres'] = len(cluster_counts)
			read_stats[queryID]['queries'] = querycount[queryID]

			queryfile = os.path.abspath('.')+"/"+"%s_centroids.fasta" %queryID #"../%s_centroids.fasta" % queryID

			unknown_seqs_dict = SeqIO.to_dict(SeqIO.parse(queryfile,'fasta'))
			for ID in unknown_seqs_dict.keys():
                                if cluster_counts.has_key(ID):
                                        unknown_seqs_dict[ID].description = "%s|%s|%s|%.2f" %(queryID, unknown_seqs_dict[ID].id, cluster_counts[unknown_seqs_dict[ID].id], float(cluster_counts[unknown_seqs_dict[ID].id])/querycount[queryID]*100)
#                               print unknown_seqs_dict[ID].description

			queries[queryID]['uc'] = os.path.abspath('.')+"/"+queryID+'.uc'

		else:
			unknown_seqs_dict = SeqIO.to_dict(SeqIO.parse(queryfile,'fasta')) #read in query sequences, atm only fasta format is supported. later I will check at this stage if there are already sequences in memory from prior quality filtering
			querycount[queryID] += len(unknown_seqs_dict)
			print "The queryfile contains %i query sequences.\n" %querycount[queryID]
		
			cluster_counts = {}
			cluster_reads = defaultdict(list)

			for ID in unknown_seqs_dict.keys():
#				print sequence
				cluster_counts[unknown_seqs_dict[ID].description] = 1
				cluster_reads[unknown_seqs_dict[ID].description] = [unknown_seqs_dict[ID].description]
#				unknown_seqs_dict[ID].description = "%s|%s|%s|%.2f" %(queryID, unknown_seqs_dict[ID].id, cluster_counts[unknown_seqs_dict[ID].id], float(cluster_counts[unknown_seqs_dict[ID].id])/querycount[queryID]*100)
				unknown_seqs_dict[ID].description = "%s|%s|%s|%.2f" %(queryID, unknown_seqs_dict[ID].id, cluster_counts[unknown_seqs_dict[ID].id], float(cluster_counts[unknown_seqs_dict[ID].id])/querycount[queryID]*100)
#				print unknown_seqs_dict[ID]
			read_stats[queryID]['queries'] = querycount[queryID]


		print "WRITING BASIC READ STATS TO FILE\n"
		read_counts_out = open("../"+args.output_prefix+"_read_stats.csv","a")
		outstring = queryID
		for m in read_metrics:
			if read_stats[queryID].has_key(m):
				outstring += ','+str(read_stats[queryID][m])
			else:
				outstring += ',NA'
		read_counts_out.write(outstring+"\n")
		read_counts_out.close()

		queries[queryID]['queryfile'] = queryfile
		print "\n%s\n" %queries[queryID]
	
		
#		print "BEFORE BLAST"
#		print cluster_counts
#		print cluster_reads
		#running blast search against previously build database
		blast_out = "%s_%s_blastn.out.xml" % (args.marker, queryID)

		some_hit = []	#This list will contain the ids of all queries that get a blast hit, i.e. even if it is non-signficant given the user specfied thresholds. For this set of reads we will attempt phylogenetic placement. No point of even trying if they dont get any blast hit.
		for approach in methods:
			if approach == 'pplacer' and len(cluster_counts) > 0:
				query_count += 1
				taxonomy_count = defaultdict(dict)
				species_count = defaultdict(list) #this is actually not needed in the pplacer parsing, but I just need to empty it here in case it still contains stuff from the blast assignemtns, NEEDS TO BE CLEANED UP LATER
				nohit_count = defaultdict(list) #again this is not needed here and I just empty it
				print "\n##### PHYLOGENETIC PLACEMENT #####\n"
				if not os.path.exists('PHYLOPLACE'):
					os.makedirs('PHYLOPLACE')
				os.chdir('PHYLOPLACE')
				print "number of queries to be processed: %i" %len(some_hit)
				print "write queries to file: \'%s_pplacer_queries.fasta\'\n" %queryID
#				print some_hit

				fas_out = open(queryID+'_pplacer_queries.fasta','w')
				for read in some_hit:
##					backup = unknown_seqs_dict[read].id
##					unknown_seqs_dict[read].id = unknown_seqs_dict[read].description
					backup = unknown_seqs_dict[read].description
					unknown_seqs_dict[read].description = unknown_seqs_dict[read].id
#					print unknown_seqs_dict[read].id
					fas_out.write(unknown_seqs_dict[read].format("fasta"))
##					unknown_seqs_dict[read].id = backup
					unknown_seqs_dict[read].description = backup

				fas_out.close()
	
				print "\n##### ALIGN QUERIES TO HMM PROFILE #####\n"
				print "detecting reference hmm profile in reference package\n"
				fh = open(args.refpkg+"/CONTENTS.json","r")
				refpkg_content = json.load(fh)
				profile = args.refpkg+'/'+refpkg_content['files']['profile']
#				print profile
				aln_fasta = args.refpkg+'/'+refpkg_content['files']['aln_fasta']
#				print aln_fasta				
				cmd = "hmmalign -o %s_queries_to_profile.sto --mapali %s %s %s_pplacer_queries.fasta" %(queryID, aln_fasta, profile, queryID)
				print cmd
				cmdlist = shlex.split(cmd)
			        stdout, stderr = subprocess.Popen(cmdlist, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate() # , stdout=subprocess.PIPE).communicate()
				print stdout	#this currently doesnt print anything, i.e. there is not stdout
				if stderr:
					print stderr
					sys.exit('some error from hmmalign')

				print "\n##### RUNNING PHYLOGENETIC PLACEMENT #####\n"
				cmd = "pplacer -c %s %s_queries_to_profile.sto -p --keep-at-most 3 -o %s.jplace" %(args.refpkg, queryID, queryID)
				print cmd+'\n'
				cmdlist = shlex.split(cmd)
			        stdout, stderr = subprocess.Popen(cmdlist, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate() # , stdout=subprocess.PIPE).communicate()
				print stdout	#this currently doesnt print anything, i.e. there is not stdout
				if stderr:
					print stderr
					sys.exit('some error from placer')
				
				print "\n##### INTERPRETING PPLACER RESULTS #####\n"
				fh = open(queryID+".jplace","r")
				jplace = json.load(fh)
#				print "\n### find all placements###"
				for pp_queries in jplace['placements']:
#					print pp_queries['p'][0] #just consider the first placement
					placement = pp_queries['p'][0][0] #take the first placement for now
#					print "add taxid \'%s\' to global taxid dictionary" %placement
					global_taxids_hit.append(placement) #add taxid of placement to global list
					
					if tax_dict[placement][1] == 'subspecies':
						rank = 'species'
					else:
						rank = tax_dict[placement][1]

					if not taxonomy_count.has_key(rank):
						taxonomy_count[rank] = defaultdict(list)
						for q in pp_queries['nm']:
							taxonomy_count[rank][tax_dict[placement][2]].append(q[0])
					else:
						for q in pp_queries['nm']:
							taxonomy_count[rank][tax_dict[placement][2]].append(q[0])
#						taxonomy_count[rank][tax_dict[placement][2]].append('dummy')
#					print tax_dict[placement][1]	#this is the taxonomy rank
#					print tax_dict[placement][2]	#this is the scientific name

#				print taxonomy_count

#				sys.exit()


			elif approach == 'blast' and len(cluster_counts) > 0:
				query_count += 1
				if not os.path.exists('BLAST_'+str(args.min_ident)):
					os.makedirs('BLAST_'+str(args.min_ident))
				os.chdir('BLAST_'+str(args.min_ident))
				print "\n### RUNNING LOCAL BLAST ###\n"

				print "running blast search against local database %s" % blast_db
#				print "QUERIES: %s" %queryfile

				blast_cmd = "blastn -query %s -db %s -evalue 1e-20 -outfmt 5 -out %s -num_threads %i -max_target_seqs 50" % (queryfile, blast_db, blast_out, args.n_threads) 
				print blast_cmd #this is just for the output, the actual blast command is run using the NCBI module below 


#				blast_cmd = "blastn -query %s -db %s -evalue 1e-20 -out %s -num_threads %i -max_target_seqs 50" % (queryfile, blast_db, "blastn.standard.out", args.n_threads) 
#				cmdlist = shlex.split(blast_cmd)
#				stdout, stderr = subprocess.Popen(cmdlist, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate() 


				blast_handle = NcbiblastxCommandline(cmd='blastn', query=queryfile, db=blast_db, evalue=1e-20, outfmt=5, out=blast_out, num_threads=args.n_threads, max_target_seqs=50)		
				stdout, stderr = blast_handle()

				print '\n'+time.strftime("%c")+'\n'
				print "\n### INTERPRETING BLAST RESULTS ###\n"
#				if 
				blast_result_handle = open(blast_out) #read in blast result (in xml format)

				print "\nparse blast xml file\n"
				blast_results = NCBIXML.parse(blast_result_handle) #parse blast xml file

				print "get top 10% hits from blast output .. ",
				res_dict = {}
				res_dict = blast_filter(b_result=blast_results, v=args.verbose, m_ident=args.min_ident, m_bitscore=args.min_bit, m_ali_length=args.min_ali_length)

				if res_dict['format'] == 'gi':
					print '\n'+time.strftime("%c")+'\n'
					print "\nfetch taxids for gis\n"
#					print "current length of taxid list: %i" %len(taxid_list)
					taxid_list = gi_to_taxid(b_filtered=res_dict, all_taxids=taxid_list, processed=gi_to_taxid_dict, v=args.verbose)
#					print "current length of taxid list: %i" %len(taxid_list)

					if len(gi_to_taxid_dict) > 0:
						print "\nwriting 'gi_to_taxid' dictionary to file %s" %args.gi_to_taxid,
						rw_gi_to_taxid_dict(dictionary=gi_to_taxid_dict, name=args.gi_to_taxid, mode='w')

					print "\ngenerate taxonomy dictionary\n"
					#tax_dict = {}
					make_tax_dict(tids=taxid_list, out_tax_dict=tax_dict, denovo_taxa=denovo_taxa, ref_taxa=reference_taxa)
				elif res_dict['format'] == 'unknown': #blast_filter function did not assign a format because no good hit was found to actually assess the format
					if not tax_dict.has_key('tax_id'): #in this rare case I check if the tax_dict is in the proper format and if not
                                        	tax_dict['tax_id'] = [] #just add the tax_id key and an empty list 

				print "\nassign taxonomy\n"
				taxonomy_count = assign_taxonomy_LCA(b_filtered=res_dict, tax_dict=tax_dict, v=args.verbose)
				###PREVIOUS VERSION HAS collected all query ids that had a hit with the database (even if it was not significant) in a list 'some_hit' limited phylogenetic placement to only these queries

			if len(cluster_counts) == 0: #this is only relevant in the rare case when no valid sequence made it throuh the trimming/clustering
				query_count += 1

			if taxonomy_count.has_key('nohit'):
				tax_dict["tax_id"].insert(0,'nohit')


			print '\n'+time.strftime("%c")+'\n'
			print "\n######## RESULT SUMMARY ########\n"
			out=open(queryID+"-results.txt","w")
			outstring="\nThe sample %s contained %i valid query sequences:\n" % (queryID, querycount[queryID])
			print outstring
			out.write(outstring+"\n")


#			print tax_dict['tax_id']
			if not tax_dict.has_key('tax_id'): #this is only relevant in the rare case when no valid sequence made it throuh the trimming/clustering and no reference taxids have been produced earlier from a custom reference
				tax_dict['tax_id'] = ['nohit']

			for tax_rank in reversed(tax_dict["tax_id"]):
#				print "The current rank to be checked is: %s" %tax_rank
				if taxonomy_count.has_key(tax_rank):
#					print "The current rank: %s" %tax_rank
#					print taxonomy_count[tax_rank]
					total_per_rank_count=int(0)
					output=[]
					for hit in sorted(taxonomy_count[tax_rank].keys()):
						if not tax_rank == 'nohit':
							global_taxids_hit.append(hit)

#						print hit
#						print "\t%s: %i" % (hit, len(taxonomy_count[tax_rank][hit])) #This is the count of unique clusters that were assigned
						current_count=int(0)
						current_reads=[]
						for read in taxonomy_count[tax_rank][hit]:
#							print "centroid: %s" % read
#							print "read count in cluster: %i" % cluster_counts[read]
##							print read.split('|')[1]
#							current_count+=cluster_counts[read.split('|')[1]]	#sum up the number of actual reads in the clusters assigned to this taxon
							current_count+=cluster_counts[read]	#sum up the number of actual reads in the clusters assigned to this taxon
							if args.extract_centroid_reads:
#								current_reads.append(read.split('|')[1])
								current_reads.append(read)

							elif args.extract_all_reads:
#								current_reads.extend(cluster_reads[read.split('|')[1]])
								current_reads.extend(cluster_reads[read])

#							total_per_rank_count+=cluster_counts[read.split('|')[1]]
							total_per_rank_count+=cluster_counts[read]
						if not tax_rank == 'nohit':
							output.append("\t%s: %i (%.2f %%)" % (tax_dict[hit][2], current_count, 100*float(current_count)/querycount[queryID]))
#						print "\t%s: %i (%.2f %%)" % (hit, current_count, 100*float(current_count)/querycount[queryID])

						#### add data to global dictionary for biom table
#						if tax_rank != 'nohit':	#dont include the nohits in the biom output
						if not global_taxa[tax_rank].has_key(hit):
							global_taxa[tax_rank][hit] = []
							global_clust[tax_rank][hit] = []
							if query_count >= 2: #if this is already the 2+ query and the taxon has not been seen so far, I need to fill up the previous samples with count 0 for this taxon
								for i in range(query_count-1):
#										print "fill: %s" %hit
									global_taxa[tax_rank][hit].append(int(0))
									global_clust[tax_rank][hit].append(int(0))

								
							global_taxa[tax_rank][hit].append(int(current_count))
							global_clust[tax_rank][hit].append(int(len(taxonomy_count[tax_rank][hit])))
						else:
							global_taxa[tax_rank][hit].append((current_count))
							global_clust[tax_rank][hit].append((len(taxonomy_count[tax_rank][hit])))
						
						### print out reads
						if current_reads: #This list is only non empty if either -e or -E was specified
#							
							if args.extract_centroid_reads:	#if the user has specified to extract centroid reads
								#the following will sort the read ids by the number of reads in the respective cluster
								temp_dict = defaultdict(list) #temporal dictionary 
								for ID in current_reads:	#populate temporal dictionary with read counts as keys and centroid ids as values in a list (in case there are clusters that happen to have the same number of reads
									temp_dict[cluster_counts[unknown_seqs_dict[ID].id]].append(ID)
#									print "%s: %s" %(cluster_counts[unknown_seqs_dict[ID].id], ID)
									
								current_reads=[]	#empty the original list
								for count in sorted(temp_dict.keys(), key=int, reverse=True):	#loop through the temporal dictionary reverse sorted by keys
#									print temp_dict[count]
									for read in temp_dict[count]:	#because the read ids are in a list I ll loop through here
#										print "%s: %s" %(count, read)
										current_reads.append(read)	#put the read ids back into the list sorted in descending order of the number reads in the respective clusters
									

							if hit == 'nohit':
								fas_out = open('nohit.fasta',"w")
							else:
								fas_out = open(tax_dict[hit][2].replace(" ", "_")+'.fasta',"w")


							for ID in current_reads:
								backup = unknown_seqs_dict[ID].id
								unknown_seqs_dict[ID].id = unknown_seqs_dict[ID].description
##								unknown_seqs_dict[ID].description = "%s|%s" %(queryID, unknown_seqs_dict[ID].id)
##								if args.extract_centroid_reads:
##									unknown_seqs_dict[ID].id = "%s|%s|%i|%.2f" % (queryID, unknown_seqs_dict[ID].id, cluster_counts[unknown_seqs_dict[ID].id], float(cluster_counts[unknown_seqs_dict[ID].id])/querycount[queryID]*100)
##									unknown_seqs_dict[ID].description = unknown_seqs_dict[ID].id
#								print unknown_seqs_dict[ID]
								fas_out.write(unknown_seqs_dict[ID].format("fasta"))
								unknown_seqs_dict[ID].id = backup
#								fas_out.write(">%s|%i\n%s\n" % ( unknown_seqs_dict[ID].id, cluster_counts[unknown_seqs_dict[ID].id], unknown_seqs_dict[ID].seq ))
							fas_out.close()

					outstring="%s: %i (%.2f %%)" % (tax_rank, total_per_rank_count, 100*float(total_per_rank_count)/querycount[queryID])
					print outstring
					out.write(outstring+"\n")
					if tax_rank!='nohit':
						for line in output:
							print line
							out.write(line+"\n")

			out.close()
			print "\n\n"
			if len(cluster_counts) > 0:
				os.chdir("../")#leave directory of current analysis
				if tax_dict["tax_id"][0] == 'nohit':
					del tax_dict["tax_id"][0] #remove the first element, i.e. 'nohit'
#			del tax_dict["tax_id"][-1] #remove the last element, i.e. 'species'
		
			print '\n'+time.strftime("%c")+'\n'

			####this is the end
#			print "\nTHIS IS THE END"
#			print '%s: %s' %(queryID,queries[queryID])
#			print "query_count: %i" %query_count
#			print 'global_taxa: %s' %global_taxa

			for rank in global_taxa: #this is just to make sure that all taxa have the same number of counts
				for taxon in global_taxa[rank]:
					if len(global_taxa[rank][taxon]) < query_count:	#if the number of entries for the current taxon is smaller than the current index
						global_taxa[rank][taxon].append(int(0))
						global_clust[rank][taxon].append(int(0))
						

#			print "TAXIDS hit by query %s: %s" %(queryID, global_taxids_hit)

		os.chdir('../')	#leave directory for query



if args.blast or args.phyloplace:

	print "\n##### DONE PROCESSING ALL SAMPLES #####"		
	print "##### FORMATTING AND WRITING BIOM OUTPUT #####\n"
	clust_to_biom = []
	data_to_biom = []
	observ_ids=[]
	for tax_rank in reversed(tax_dict["tax_id"]):
		if global_taxa.has_key(tax_rank):
#			print tax_rank
#			print global_taxa[tax_rank]
			for out in sorted(global_taxa[tax_rank]):
#				print "%s: %s" %(out, global_taxa[tax_rank][out])
				observ_ids.append(tax_dict[out][2])
				data_to_biom.append(global_taxa[tax_rank][out])
				clust_to_biom.append(global_clust[tax_rank][out])

	#take care of nohits if any exist
	if global_taxa.has_key('nohit'):
		observ_ids.append('unassigned')
		data_to_biom.append(global_taxa['nohit']['nohit'])
		clust_to_biom.append(global_clust['nohit']['nohit'])

	#print data_to_biom
	data = np.asarray(data_to_biom)
	clust_data = np.asarray(clust_to_biom)
#	print "data:\n%s" %data
#	print len(data)

	sample_ids = []	
	for key in sorted(queries.keys()):
#		print key
		for met in methods:
#			print met
			sample_ids.append(key+"."+met)
#	print "sample_ids:\n%s" %sample_ids
	#print len(sample_ids)

	sample_metadata=[]
	for q in sample_ids:
		temp={}
		temp['method'] = q.split(".")[-1]
#		print queries[q]
#		sample_metadata.append(queries[q])
		if args.mock_meta_data:
			r="%.1f" %random.uniform(20.0,25.0)
			temp['temperature'] = "%.1f C" %float(r)
			r="%.1f" %random.uniform(10.0,15.0)
			temp['depth'] = "%.1f m" %float(r)
			treatments = ['A','B']
			temp['treatment'] = random.choice(treatments)

		if args.metadata:
			current = ".".join(q.split(".")[:-1])
#			print current
			if not metadata.has_key(current):
				print "sample %s has no metadata available\n"
			else:
				for meta in metadata[current].keys():
					temp[meta] = metadata[current][meta]
		sample_metadata.append(temp)
		
#	print "sample_metadata:\n%s" %sample_metadata
#	print len(sample_metadata)
	
#	print "observ_ids:\n%s" %observ_ids
#	print len(observ_ids)

	observation_metadata=[]

	Taxonomy=defaultdict(dict)
	syn = {'kingdom': 'k__', 'phylum': 'p__', 'class': 'c__', 'order': 'o__', 'family': 'f__', 'genus':'g__', 'species': 's__'}
	levels = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
	level_indices = []
        for level in levels:
        	for i in range(len(tax_dict['tax_id'])):
                	if level == tax_dict['tax_id'][i]:
                        	level_indices.append(i)

	global_taxids_hit = list(set(global_taxids_hit))
	for tid in global_taxids_hit:
		ind_taxonomy=[]
		if tax_dict.has_key(tid):
			for i in range(len(levels)):
				if tax_dict[tid][level_indices[i]]:
					ind_taxonomy.append('%s%s' %(syn[levels[i]],tax_dict[tax_dict[tid][level_indices[i]]][2].replace(' ','_')))
				else:
					ind_taxonomy.append('%sunknown' %syn[levels[i]])
		else:
			"taxid %s has not been seen before\n" %tid

		for j in reversed(range(len(levels))):
			if '_unknown' in ind_taxonomy[j]:
				del ind_taxonomy[j]

		Taxonomy[tax_dict[tid][2]]['taxonomy'] = ind_taxonomy

##	print Taxonomy

	if observ_ids[-1] == 'unassigned':
		ind_taxonomy = []
		for lev in levels:
			ind_taxonomy.append('%sunassigned' %syn[lev])

		Taxonomy['unassigned']['taxonomy'] = ind_taxonomy

	for taxon in observ_ids:
#		print "%s: %s" %(taxon, Taxonomy[taxon])
		observation_metadata.append(Taxonomy[taxon])
	
#	print "observation metadata:\n%s" %observation_metadata
#	print len(observation_metadata)

### double check all inputs for table
#	print "data: %s" %data
#	print len(data)
#	print "cluster data: %s" %clust_data
#	print len(clust_data)
#	print "observed_ids: %s" %observ_ids
#	print len(observ_ids)
#	print "sample ids: %s" %sample_ids
#	print len(sample_ids)
#	print "observation metadata: %s" %observation_metadata
#	print len(observation_metadata)
#	print "sample_metadata: %s" %sample_metadata
#	print len(sample_metadata)

#	print "\n"

	table = Table(data, observ_ids, sample_ids, observation_metadata, sample_metadata, table_id='Read count Table')
	clust_table = Table(clust_data, observ_ids, sample_ids, observation_metadata, sample_metadata, table_id='Cluster count Table')


#	print tables
	out=open(args.output_prefix+".biom","w")
	table.to_json('metaBEAT v.'+VERSION, direct_io=out)
	out.close()

	out=open(args.output_prefix+".tsv","w")
	out.write(table.to_tsv(header_key='taxonomy', header_value='taxomomy')) #to_json('generaged by test', direct_io=out)
	out.close()

	out=open(args.output_prefix+"_clusters.biom","w")
	clust_table.to_json('metaBEAT v.'+VERSION, direct_io=out)
	out.close()

	out=open(args.output_prefix+"_clusters.tsv","w")
	out.write(clust_table.to_tsv(header_key='taxonomy', header_value='taxomomy')) #to_json('generaged by test', direct_io=out)
	out.close()
print "\n##### DONE! #####\n"
print '\n'+time.strftime("%c")+'\n'

#print "remove read-pair id from extended reads. \n output overall run summary, i.e. per sample: raw reads, trimmed reads, merged reads, clusters, etc. \n make OTU table output standard"
#
#add functionality to iterate over parameters, e.g. cluster_coverage 10,20,30,40,50 and treat each as different sample in the biom file
