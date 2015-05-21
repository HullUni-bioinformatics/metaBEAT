#! /usr/bin/python

from Bio import SeqIO
from Bio import Entrez
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast import NCBIXML
from Bio.Alphabet import generic_dna

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

informats = {'gb': 'gb', 'genbank': 'gb', 'fasta': 'fasta', 'fa': 'fasta', 'fastq': 'fastq'}
all_seqs = []
skipped_ref_seqs = [] #defaultdict(list)
references = {}
queries = defaultdict(list)
seq_info = ['"seqname","accession","tax_id","species_name","is_type"']
records = {}
reference_taxa = {}
taxids = defaultdict(int)
date = time.strftime("%d-%b-%Y").upper() 

parser = argparse.ArgumentParser(description='metaBEAT - metaBarcoding and Environmental DNA Analyses tool')
#usage = "%prog [options] REFlist"
#parser = argparse.ArgumentParser(usage='%(prog)s [options] REFlist', formatter_class=RawTextHelpFormatter)

parser.add_argument("-Q", "--querylist", help="file containing a list of query files", metavar="<FILE>", action="store")
parser.add_argument("-v", "--verbose", help="turn verbose output on", action="store_true")
parser.add_argument("-s", "--seqinfo", help="write out seq_info.csv file", action="store_true")
parser.add_argument("-f", "--fasta", help="write out ref.fasta file", action="store_true")
parser.add_argument("-p", "--phyloplace", help="perform phylogenetic placement", action="store_true")
parser.add_argument("-t", "--taxids", help="write out taxid.txt file", action="store_true")
parser.add_argument("-b", "--blast", help="compile local blast db and blast queries", action="store_true")
parser.add_argument("-m", "--marker", help="marker ID (default: marker)", metavar="<string>", action="store", default="marker")
parser.add_argument("-n", "--n_threads", help="Number of threads (default: 1)", type=int, metavar="<INT>", action="store", default="1")
parser.add_argument("-E", "--extract_centroid_reads", help="extract centroid reads to files", action="store_true")
parser.add_argument("-e", "--extract_all_reads", help="extract reads to files", action="store_true")
query_group = parser.add_argument_group('Query preprocessing', 'The parameters in this group affect how the query sequences are processed')
query_group.add_argument("--trim_adapter", help="trim adapters provided in file", metavar="<FILE>", action="store")
query_group.add_argument("--trim_qual", help="minimum phred quality score (default: 30)", metavar="<INT>", type=int, action="store", default=30)
query_group.add_argument("--trim_window", help="sliding window size (default: 5) for trimming; if average quality drops below the specified minimum quality all subsequent bases are removed from the reads", metavar="<INT>", type=int, action="store", default=5)
query_group.add_argument("--trim_minlength", help="minimum length of reads to be retained after trimming (default: 50)", metavar="<INT>", type=int, action="store", default=50)
query_group.add_argument("--merge", help="attempt to merge paired-end reads", action="store_true")
query_group.add_argument("--product_length", help="estimated length of PCR product (default: 100)", metavar="<INT>", type=int, action="store", default=100)
query_group.add_argument("--phred", help="phred quality score offset - 33 or 64 (default: 33)", metavar="<INT>", type=int, action="store", default=33)
reference_group = parser.add_argument_group('Reference', 'The parameters in this group affect the reference to be used in the analyses')
reference_group.add_argument("-R", "--REFlist", help="file containing a list of files to be used as reference sequences", metavar="<FILE>", action="store")
reference_group.add_argument("--gb_out", help="output the corrected gb file", metavar="<FILE>", action="store", default="")
reference_group.add_argument("--rec_check", help="check records to be used as reference", action="store_true")
cluster_group = parser.add_argument_group('Query clustering options', 'The parameters in this group affect read clustering')
cluster_group .add_argument("--clust_match", help="identity threshold for clustering in percent (default: 1)", type=float, metavar="<FLOAT>", action="store", default="1")
cluster_group .add_argument("--clust_cov", help="minimum number of records in cluster (default: 1)", type=int, metavar="<INT>", action="store", default="1")
blast_group = parser.add_argument_group('BLAST search', 'The parameters in this group affect BLAST search and BLAST based taxonomic assignment')
blast_group.add_argument("--www", help="perform online BLAST search against nt database", action="store_true")
blast_group.add_argument("--min_ident", help="minimum identity threshold in percent (default: 0.95)", type=float, metavar="<FLOAT>", action="store", default="0.95")
blast_group.add_argument("--min_bit", help="minimum bitscore (default: 80)", type=int, metavar="<INT>", action="store", default="80")
parser.add_argument("--version", action="version", version='%(prog)s v.0.5')
args = parser.parse_args()

####START OF MAIN PROGRAM
print '\n'+time.strftime("%c")+'\n'
print "%s\n" % (' '.join(sys.argv))

if args.trim_adapter:
	args.trim_adapter = os.path.abspath(args.trim_adapter)
	if not os.path.isfile(args.trim_adapter):
		print "adapter file is not a valid file"
		sys.exit(0)

if not os.path.isfile(args.REFlist):
	print "no valid reference file supplied\n"
	sys.exit(0)
else:
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
			if not informats.has_key(querydata[1]):
				print "query file format for sample %s is invalid" % querydata[0]
				sys.exit(0) 
			for f in querydata[2:]:
				if not os.path.isfile(f):
					print "%s is not a valid file" % f
					sys.exit(0)
			
			queries[querydata[0]] = querydata[1:]
	

#print references
for reffile, refformat in references.items():
	seqs=list(SeqIO.parse(reffile,refformat, alphabet=generic_dna))
	print "processing %s (containing %i records)" % (reffile,len(seqs))
	for i in reversed(range(len(seqs))):
#		print seqs[i].description
		if args.rec_check:
			if not seqs[i].id in records:	#records list keeps track of all records.
				records[seqs[i].id] = 1
			else:				#if a record id occurs again, say in the case of custom fasta sequences the record id is the first name in the fasta header, e.g. the Genus name
				records[seqs[i].id] += 1
				seqs[i].id += "_%s" % records[seqs[i].id]

			if not seqs[i].features:		#if the record does not have a source feature yet. Will happen for records originating from custom fasta files
#				print "the record has no features\n"	
				seqs[i].features.append(SeqFeature(FeatureLocation(0,len(seqs[i].seq),strand=1), type="source"))	#create source feature
				seqs[i].features[0].qualifiers['original_desc'] = seqs[i].description
			
			if not seqs[i].features[0].qualifiers.has_key('organism'):	#if the record does not have an organism key. Again will happen for custom fasta sequences
				seqs[i].features[0].qualifiers['organism'] = ["%s %s" % (seqs[i].description.split(" ")[0], seqs[i].description.split(" ")[1])]	#add an organism key and make it the first 2 words in the record description will are per default the first 2 words in the fasta header

			if not seqs[i].features[0].qualifiers['organism'][0] in reference_taxa:	#if the organism name has not yet been encountert the following lines will find the corresponding taxid, add it to the record and store it in the reference_taxa dictionary
#				print "seeing %s for the first time" %seqs[i].features[0].qualifiers['organism'][0]
				handle = Entrez.esearch(db="Taxonomy", term=seqs[i].features[0].qualifiers['organism'][0])	#search the taxonomy database for the taxon by organism name
				taxon = Entrez.read(handle)
				if not taxon["IdList"]:	#if the search term has not yielded any result
					print "WARNING: %s in record '%s' was not recognized as a valid species name -> record will be ommitted from further analyses.\nI'll put the sequence record in the file 'skipped_seqs.fasta' in case you want to fix it." % (seqs[i].features[0].qualifiers['organism'][0], seqs[i].name)
					skipped_ref_seqs.append(seqs.pop(i))

				else:
					taxid = taxon["IdList"][0]	#the taxid is the first element in the dictionary that is returned from Entrez
					seqs[i].features[0].qualifiers['db_xref'] = ["taxon:" + taxid]	#update the db_rxref qualifier with the correct taxid
					taxids[taxid] += 1
					current = '"%s","%s","%s","%s",0' % (seqs[i].id, seqs[i].id, taxid, seqs[i].features[0].qualifiers['organism'][0]) #producing a string for the current record for the seq_info file that is needed for the reference package required by pplacer
					seq_info.extend([current])	#add the string to the seq_info array
#					print "adding %s with taxid: %s to the dictionary" %(seqs[i].features[0].qualifiers['organism'][0], taxid)
					reference_taxa[seqs[i].features[0].qualifiers['organism'][0]] = taxid
			else:
				seqs[i].features[0].qualifiers['db_xref'] = ["taxon:" + reference_taxa[seqs[i].features[0].qualifiers['organism'][0]]]
#				print "Have seen %s before. The taxid is: %s" %(seqs[i].features[0].qualifiers['organism'][0], reference_taxa[seqs[i].features[0].qualifiers['organism'][0]])

			if not seqs[i].annotations.has_key('date'):	#set the date in the record if it does not exist
				seqs[i].annotations['date'] = date
				
		else:
			taxid = seqs[i].features[0].qualifiers['db_xref'][0].split(":")[1]
#			print taxid
			taxids[taxid] += 1

	if skipped_ref_seqs:
		SKIPPED = open("skipped_seqs.fasta","w")
		for rec in skipped_ref_seqs:
			outstring = ">%s\n%s" %(rec.description, rec.seq)
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


if args.taxids:
	print "write out taxids to taxids.txt\n"
	f = open("taxids.txt","w")
	for key in taxids.keys():
#		print key
		f.write(key + "\n")
	f.close()
	cmd = "taxit taxtable -d /home/chrishah/src/taxtastic/taxonomy_db/taxonomy.db -t taxids.txt"# -o taxa.csv"
	print "running taxit to generate reduced taxonomy table"
	print cmd

#	handle = subprocess.call(taxtable, shell=True)
	taxtable,err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	if args.phyloplace:
		f = open("taxa.csv","w")
		f.write(taxtable + "\n")
		f.close()

	tax_dict = {}
	for line in taxtable.split("\n"):
#		print line
		line = re.sub('"','',line)
#		print line
		array = line.split(',')
		if len(array) > 1:
#			print len(array)
#			print array
			key = array.pop(0)
			array.pop(-1)
#			print "key: %s - %s" % (key, array[-1])
#			print len(array)
			tax_dict[key]=array

#for key,value in tax_dict.items():
#	print "%s: %s" % (key, value)

if args.seqinfo:
	print "write out seq_info.csv\n"
	f = open("seq_info.csv","w")
	for elem in seq_info:
#		print elem
		f.write(elem + "\n")

if args.fasta:	#this bit writes out the sequences that will become the reference dataset
	print "write out reference sequences to refs.fasta\n"
	OUT = open("refs.fasta","w")
#	OUT_temp = open ("temp_refs.fasta","w")
	for record in all_seqs:
#		outstring = ">%s\n%s" % (record.features[0].qualifiers['db_xref'][0].split(":")[1], record.seq) #note that the header of the sequence will be the taxid
		outstring = ">%s|%s\n%s" % (record.id, record.features[0].qualifiers['db_xref'][0].split(":")[1], record.seq) #note that the header of the sequence will be the taxid
		
		outstring2 = ">%s\n%s" % (record.id, record.seq)
		OUT.write(outstring + "\n")
#		OUT_temp.write(outstring2 + "\n")
	OUT.close()
#	OUT_temp.close()

if args.blast:
	print "building blast db for marker %s\n" % args.marker
	makeblastdb="makeblastdb -in refs.fasta -dbtype nucl -out %s_blast_db" % args.marker
	cmdlist = shlex.split(makeblastdb)
	cmd = subprocess.Popen(cmdlist, stdout=subprocess.PIPE)
	stdout = cmd.communicate()[0]
	print stdout


print '\n'+time.strftime("%c")+'\n'
querycount = defaultdict(int)
if args.blast:
	for queryID, querydata in sorted(queries.items()): #loop through the query files and read in the current query id and the path to the files
		print "\nprocessing query ID: %s\n###############\n" % (queryID)
		species_count = defaultdict(list)
		taxonomy_count = defaultdict(dict)
		nohit_count = defaultdict(list)
		if not os.path.exists(queryID):
			os.makedirs(queryID)

		os.chdir(queryID)
		if (querydata[0]=="fastq"):
			if len(querydata)==3:
				trimmomatic_path="java -jar /usr/bin/trimmomatic-0.32.jar"
				print "\ntrimming PE reads with trimmomatic"
				trimmomatic_exec=trimmomatic_path+" PE -threads %i -phred%i %s %s %s_forward.paired.fastq.gz %s_forward.singletons.fastq.gz %s_reverse.paired.fastq.gz %s_reverse.singletons.fastq.gz ILLUMINACLIP:%s:5:5:5 TRAILING:%i LEADING:%i SLIDINGWINDOW:%i:%i MINLEN:%i" % (args.n_threads, args.phred, querydata[1], querydata[2], queryID, queryID, queryID, queryID, args.trim_adapter, args.trim_qual, args.trim_qual, args.trim_window, args.trim_qual, args.trim_minlength)
				trimmed_files=[queryID+'_forward.paired.fastq.gz', queryID+'_forward.singletons.fastq.gz', queryID+'_reverse.paired.fastq.gz', queryID+'_reverse.singletons.fastq.gz']
			elif len(querydata)==2:
				print "\ntrimming SE reads with trimmomatic"
				trimmomatic_exec= trimmomatic_path+" SE -threads %i -phred%i %s %s_trimmed.fastq.gz ILLUMINACLIP:%s:5:5:5 TRAILING:%i LEADING:%i SLIDINGWINDOW:%i:%i MINLEN:%i" % (args.n_threads, args.phred, querydata[1], queryID, args.trim_adapter, args.trim_qual, args.trim_qual, args.trim_window, args.trim_qual, args.trim_minlength)
				trimmed_files=[queryID+'_trimmed.fastq.gz']

			trimmomatic="%s" % trimmomatic_exec
			print trimmomatic
			cmdlist = shlex.split(trimmomatic)
			cmd = subprocess.Popen(cmdlist, stdout=subprocess.PIPE, stderr=subprocess.PIPE)	
			stdout,stderr = cmd.communicate() #[0]
			print stdout
			print stderr

			if len(trimmed_files)==4:
				if args.merge:
					print "merging paired-end reads with flash\n"
					cmd="flash %s %s -M %s -t %i -p %i -o %s -z" % ( trimmed_files[0], trimmed_files[2], args.product_length, args.n_threads, args.phred, queryID)
					print cmd
					cmdlist = shlex.split(cmd)
					cmd = subprocess.Popen(cmdlist, stdout=subprocess.PIPE, stderr=subprocess.PIPE)	
					stdout,stderr = cmd.communicate() #[0]
					print stdout
					print stderr
					trimmed_files[0]= queryID+'.notCombined_1.fastq.gz'
					trimmed_files[2]= queryID+'.notCombined_2.fastq.gz'
					trimmed_files.insert(0, queryID+'.extendedFrags.fastq.gz')

				files = " ".join(trimmed_files[-2:])
				cmd="zcat %s | fastq_to_fasta | fastx_reverse_complement > temp2.fasta" % (files) # trimmed_files[2], trimmed_files[3])
				print cmd
				cmdlist = shlex.split(cmd)
				cmd = subprocess.call(cmd, shell=True)
				
				files = " ".join(trimmed_files[:-2])
				cmd="zcat %s | fastq_to_fasta > temp1.fasta" % (files) #trimmed_files[0], trimmed_files[1])
				print cmd
				cmdlist = shlex.split(cmd)
				cmd = subprocess.call(cmd, shell=True)

				cmd="cat temp1.fasta temp2.fasta > temp_trimmed.fasta"
				print cmd
				cmdlist = shlex.split(cmd)
				cmd = subprocess.call(cmd, shell=True)
				
				os.remove("temp1.fasta")
				os.remove("temp2.fasta")
			else:
				cmd="zcat %s | fastq_to_fasta > temp_trimmed.fasta" % trimmed_files[0]
				print cmd
				cmdlist = shlex.split(cmd)
				cmd = subprocess.call(cmd, shell=True)

			fin=open("temp_trimmed.fasta","r")
			f_out=open(queryID+'_trimmed.fasta',"w")
			for line in fin:
				f_out.write(line.replace(" ","_"))
				
			fin.close()
			f_out.close()
			os.remove("temp_trimmed.fasta")

#			cmd="zcat %s | fastx_clipper -a CTAGAGGAGCCTGTTCTA | fastq_to_fasta > temp1.fasta" % (querydata[1])
#			print cmd
#			cmdlist = shlex.split(cmd)
#                       cmd = subprocess.call(cmd, shell=True)

#			cmd="zcat %s | fastx_clipper -a GGGGTATCTAATCCCAGT | fastq_to_fasta | fastx_reverse_complement > temp2.fasta" % (querydata[2])
#			print cmd
#			cmdlist = shlex.split(cmd)
#                       cmd = subprocess.call(cmd, shell=True)

#			cmd="cat temp1.fasta temp2.fasta > temp_concat.fasta"
#			print cmd
#			cmdlist = shlex.split(cmd)
#                       cmd = subprocess.call(cmd, shell=True)


			queryfile=queryID+'_trimmed.fasta'

		unknown_seqs_dict = SeqIO.to_dict(SeqIO.parse(queryfile,'fasta'))
		unknown_seqs=list(SeqIO.parse(queryfile,'fasta'))	#read in query sequences, atm only fasta format is supported. later I will check at this stage if there are already sequences in memory from prior quality filtering
		querycount[queryID] += len(unknown_seqs)

		##running clustering
		print "\nclustering using vsearch"
		cmd = "vsearch --cluster_fast %s --id %.2f --threads %s --centroids %s_centroids.fasta --uc %s.uc" % (queryfile, args.clust_match, args.n_threads, queryID, queryID )
		print cmd
		cmdlist = shlex.split(cmd)
	        stdout, stderr = subprocess.Popen(cmdlist, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate() # , stdout=subprocess.PIPE).communicate()
		if stdout:
			print stdout

		#read in the results from the clustering, i.e. number of reads per cluster, later I could read in the actual read ids to allow for retrievel of reads assigned to a given taxon
		cluster_counts = {}
		cluster_reads = defaultdict(list)
		all_clust_count = int(0)
		f=open(queryID+".uc","r") #read the file
		for line in [l.strip() for l in f]: #loop through the file one line at a time, stripping off any newline characters
			if line.startswith("C"): #process only lines that start with a "C"
				all_clust_count+=1
				elem = line.split("\t")	#split the lines at tab
				if int(elem[2]) >= args.clust_cov:
					cluster_counts[elem[8]] = int(elem[2]) #write the counts to dictionary with key being the id of the centroid read
			if args.extract_all_reads:
				if line.startswith("H"):
					elem = line.split("\t")
					if cluster_reads.has_key(elem[9]):
						cluster_reads[elem[9]].append(elem[8]) #add the new read id to the centroid cluster
					else:
						cluster_reads[elem[9]] = [elem[9]] #create a new key for the centroid id and add the centroid id as the first element into the list
						cluster_reads[elem[9]].append(elem[8]) #add the new read id to the centroid cluster
		f.close()
		if args.clust_cov>1:
			os.rename(queryID+'_centroids.fasta', queryID+'_centroids_backup.fasta')
			f=open(queryID+'_centroids.fasta',"w")
			seqs=list(SeqIO.parse(queryID+'_centroids_backup.fasta','fasta'))
			for record in seqs:
				if cluster_counts.has_key(record.id):
					outstring=">%s\n%s\n" % (record.id, record.seq)
					f.write(outstring)

			querycount[queryID] -= all_clust_count - len(cluster_counts)
			f.close()
		
		print "vsearch identified %i clusters (clustering threshold %.2f) - %i clusters (minimum of %i records per cluster) are used in subsequent analyses\n" % (all_clust_count, float(args.clust_match), len(cluster_counts), args.clust_cov)

		#running blast search against previously build database
		queryfile = "%s_centroids.fasta" % queryID
		blast_db = "../%s_blast_db" % args.marker
		blast_out = "%s_%s_blastn.out.xml" % (args.marker, queryID)

		print "running blast search against local database %s" % blast_db
		blast_cmd = "blastn -query %s -db %s -evalue 0.001 -outfmt 5 -out %s -num_threads %i -max_target_seqs 50" % (queryfile, blast_db, blast_out, args.n_threads) 
		print blast_cmd #this is just for the output, the actual blast command is run using the NCBI module below 


		blast_cmd = "blastn -query %s -db %s -evalue 0.001 -out %s -num_threads %i" % (queryfile, blast_db, "blastn.standard.out", args.n_threads) 
		cmdlist = shlex.split(blast_cmd)
		stdout, stderr = subprocess.Popen(cmdlist, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate() 


		blast_handle = NcbiblastxCommandline(cmd='blastn', query=queryfile, db=blast_db, evalue=0.001, outfmt=5, out=blast_out, num_threads=args.n_threads, max_target_seqs=50)		
		stdout, stderr = blast_handle()
		blast_result_handle = open(blast_out) #read in blast result (in xml format)


		#parse blast xml file
		blast_results = NCBIXML.parse(blast_result_handle) #parse blast xml file
		
		processed=[]
		for res in blast_results: #loop through the blast result query after query
#			print res.ka_params
#			print dir(res)
			print "\nquery: %s" % res.query #the current query
			if not res.alignments:	#if no alignment was found for the query
				print "no hit - done"
				processed.append(res.query)
				nohit_count['no_hit'].append(res.query) #append the query id to the dictionary that contains all the nohit reads
			elif (float(res.alignments[0].hsps[0].identities)/len(res.alignments[0].hsps[0].query) < args.min_ident) or res.alignments[0].hsps[0].bits < 80:
				print "no significant hit - done"
				processed.append(res.query)
				nohit_count['no_hit'].append(res.query) #append the query id to the dictionary that contains all the nohit reads

			else: #if a hit has been found
				max_bit_score = res.alignments[0].hsps[0].bits #record the maximum bitscore
#				print max_bit_score
#				print "%f" % (max_bit_score*0.9)
				hit_taxids = defaultdict(int)	#define dictionaries
				perfect_hit_taxids = defaultdict(int)
				ambig_hits_taxids = defaultdict(int)

				for alignment in res.alignments: #for each hit of the current query
					for hsp in alignment.hsps:	#for each alignment the parser creates a list through which we loop here
						perc_ident = 100*hsp.identities/float(res.query_letters)	#calculate percent identity
						if (res.alignments[0].hsps[0].identities == len(res.alignments[0].hsps[0].query) and (res.query_letters==hsp.identities)):	#if the query length is equal to the number of hit bases, we have a perfect match
#							perfect_hit_taxids[alignment.title.split(" ")[1]] += 1	#count the number of perfect hits to a given taxon, the name (which is the taxid of the taxon, because I have build the database using taxids as sequence headers) is recorded in the alignment.title after the first space
							perfect_hit_taxids[alignment.title.split(" ")[1].split("|")[1]] += 1
#							print "perfect hit recorded "
						if hsp.bits > (max_bit_score*0.9): # and (float(hsp.identities)/len(hsp.query)*100) > 95:	#if a hit has a bitscore that falls within the top 90 % of the bitscores recorded
#							hit_taxids[alignment.title.split(" ")[1]] += 1	#count the number of top90 hits to a given taxon
							hit_taxids[alignment.title.split(" ")[1].split("|")[1]] += 1	#count the number of top90 hits to a given taxon
#							print "top 90 hit recorded"
#						else:
#							print "not recorded because a hit of %f is below the threshold of %f" % (hsp.bits, max_bit_score*0.9)
			
				if len(perfect_hit_taxids)==1:	#if the perfect hit dictionary contains only one key, that means that only one taxon has been hit perfectly
					print "1 perfect match - done:\n%s" % tax_dict[perfect_hit_taxids.keys()[0]][2]
					processed.append(res.query)
#					print "query %s (cluster countains %i sequences) assigned to %s (%s; %s)" % (res.query, cluster_counts[res.query], perfect_hit_taxids.keys()[0], tax_dict[perfect_hit_taxids.keys()[0]][2], tax_dict[perfect_hit_taxids.keys()[0]][1])
#					print "this should be species: %s" % tax_dict[perfect_hit_taxids.keys()[0]][2]
					species_count[tax_dict[perfect_hit_taxids.keys()[0]][2]].append(res.query) #add the perfect hit, which is bound to be a species to a dictionary. The key is the species name (which is fetched from the tax_dict dictionary based on the taxid) and the id of the query that had a perfect hit will be added to a list in the value of the dictionary
				else:	#if the perfect hit dictionary contains zero or more than one keys, i.e. more than one species in the database had perfect hits
#					print "reached ambiguous hits"
					if len(perfect_hit_taxids)>1:	#if more than one taxids had perfect hits
#						print "%i perfect matches" % len(perfect_hit_taxids)
#						print perfect_hit_taxids
						ambig_hits_taxids.update(perfect_hit_taxids)	#make the list of perfectly hitting taxids the current list to be fed into the LCA part below
					else:	#if no perfect hit had been found
#						print "no perfect matches"
#						print hit_taxids
						if len(hit_taxids)==1:	#if only one top 90% hit has been found
							print "1 top90 match - done:\n%s" % tax_dict[hit_taxids.keys()[0]][2]
							processed.append(res.query)
#							print "query %s (cluster contains %i sequences) assigned to %s (%s; %s)" % (res.query, cluster_counts[res.query], hit_taxids.keys()[0], tax_dict[hit_taxids.keys()[0]][2], tax_dict[hit_taxids.keys()[0]][1])
#							print "this should be species: %s" % tax_dict[hit_taxids.keys()[0]][2]
							species_count[tax_dict[hit_taxids.keys()[0]][2]].append(res.query) #add the one top90 hit to a dictionary. The key is the species name (which is fetched from the tax_dict dictionary based on the taxid) and the id of the query that had a perfect hit will be added to a list in the value of the dictionary
						else:	#if more than one top90 hits were found
#							print "%i top90 matching species" % len(hit_taxids)
							ambig_hits_taxids.update(hit_taxids)	#make this the list of taxids to be fed to the LCA part below

					print "number of amgibuous hits: %i" % len(ambig_hits_taxids)
					if len(ambig_hits_taxids)>0: #if the list contains more than one taxids, we attempt to identify the Lowest common ancestor
#						print "identifying LCA"
						tax_list = []	#defining variables
						parent_count = defaultdict(int)
						counter=1
						no_species_hit = defaultdict(int)
					
						for counter in range(len(tax_dict["tax_id"])):
#						while counter <= len(tax_dict["tax_id"]): #loop through the taxonomic levels encountert (see tax_dict has been generated with taxtastic above). The key "tax_id" contains a list of taxonomic ranks
							counter += 1
#							print counter
							index=counter*-1
#							print "current index: %i" % index
							for hit,count in ambig_hits_taxids.items(): #generate a list of taxids that have been hit and were retained by the above criteria
#								print "%s: %s" % (tax_dict[hit][2], count)
								tax_list.extend([hit])
	
#							print "list of taxids to be included in LCA analysis:\n%s" % tax_list

							taxok=int(0)
							for tax in tax_list: #for every taxid in the list
#								print "taxon: %s; taxid at level %i: %s" % ( tax, counter, tax_dict[tax][index])
								if tax_dict[tax][index]: 
								#we access the tax_dict dictionary using the taxid as the key, the value is a list, which contains taxids 
								#for the parents of the current taxid. The last element in the list is a taxid at species level, the next to last
								#contains a taxid for the corresponding genus level, the next to that the taxid for family and so on
								#as you go to higher taxonomic levels some levels will be empty because in one lineage there is subfamily rank
								#and in the other this field is just empty.
								#Here I just check if there is valid taxids at the current taxonomic level (via the index which reads through the list backwards
								#for each of the taxids in the list
									taxok += 1
									parent_count[tax_dict[tax][index]] += 1
#								else:
#									print "empty"

#							print "number of valid unique taxids: %i" % len(parent_count)
							if taxok == len(tax_list):
								if len(parent_count)==1:
#									print "query %s was assigned to LCA %s (%s; %s)" % (res.query, parent_count.keys()[0], tax_dict[parent_count.keys()[0]][2], tax_dict[parent_count.keys()[0]][1])
									print "assigned to LCA - done\n%s" % tax_dict[parent_count.keys()[0]][2]
									processed.append(res.query)
#
									if not taxonomy_count.has_key(tax_dict[parent_count.keys()[0]][1]): #if the taxonomic rank has not been encountert before:
#										print "taxonomic rank %s has not been encountert before" % tax_dict[parent_count.keys()[0]][1]
#										print "ading new empty dictionary %s" % tax_dict[parent_count.keys()[0]][1]
										
										taxonomy_count[tax_dict[parent_count.keys()[0]][1]] = defaultdict(list) #empty dictionary
										taxonomy_count[tax_dict[parent_count.keys()[0]][1]][tax_dict[parent_count.keys()[0]][2]].append(res.query)
										
#This worked for counts								taxonomy_count[tax_dict[parent_count.keys()[0]][1]] = defaultdict(int) #empty dictionary
#This worked for counts								taxonomy_count[tax_dict[parent_count.keys()[0]][1]][tax_dict[parent_count.keys()[0]][2]] += cluster_counts[res.query]

#										print "level 1: %s" % taxonomy_count.keys()
									#	for rank in taxonomy_count.keys():
#											print "level 2 (key: %s): %s" % (rank, taxonomy_count[rank].keys())
									#		for count in taxonomy_count[rank].keys():
#												print "level 3 (key %s): %s" % (count, taxonomy_count[rank][count])
									else:
#										print "I have seen rank %s before" % ( tax_dict[parent_count.keys()[0]][1] )
#										print "currently the rank %s contains the following keys: %s" % ( tax_dict[parent_count.keys()[0]][1], taxonomy_count[tax_dict[parent_count.keys()[0]][1]].keys())
#										print "check if the rank %s already contains the key: %s" % (tax_dict[parent_count.keys()[0]][1], tax_dict[parent_count.keys()[0]][2] )
										if not taxonomy_count[tax_dict[parent_count.keys()[0]][1]].has_key(tax_dict[parent_count.keys()[0]][2]):
#											print "not found"
#											print "adding new key %s to rank %s" % ( tax_dict[parent_count.keys()[0]][2], tax_dict[parent_count.keys()[0]][1])
#This worked for counts									taxonomy_count[tax_dict[parent_count.keys()[0]][1]][tax_dict[parent_count.keys()[0]][2]] += cluster_counts[res.query]
											taxonomy_count[tax_dict[parent_count.keys()[0]][1]][tax_dict[parent_count.keys()[0]][2]].append(res.query)
										else:
#											print "found"
#											print "current count for rank %s, key %s is: %s" % ( tax_dict[parent_count.keys()[0]][1], tax_dict[parent_count.keys()[0]][2], taxonomy_count[tax_dict[parent_count.keys()[0]][1]][tax_dict[parent_count.keys()[0]][2]] )
#											print "add %i to count" % cluster_counts[res.query]
#This worked for counts									taxonomy_count[tax_dict[parent_count.keys()[0]][1]][tax_dict[parent_count.keys()[0]][2]] += cluster_counts[res.query]
											taxonomy_count[tax_dict[parent_count.keys()[0]][1]][tax_dict[parent_count.keys()[0]][2]].append(res.query)
#										print taxonomy_count

									break #if LCA is found break out of loop
#								else:
#									print "There was at least one invalid taxid encountert at this level"
#							counter += 1
#							print "no LCA found at level %i" % counter
							tax_list = []
							parent_count = defaultdict(int)
							if counter == len(tax_dict["tax_id"]):
								print "no LCA could be found for query: %s" % res.query

		print "number of clusters processed: %i (of %i)" % ( len(processed), len(cluster_counts))
		if len(processed)!=len(cluster_counts):	#final check
			print "not all clusters were properly processed"
			print "%i of %i" % (len(processed), len(cluster_counts))

		taxonomy_count['species'] = species_count
		taxonomy_count['nohit'] = nohit_count
#		print "\nThe final dictionary"
#		print taxonomy_count
#		print
		tax_dict["tax_id"].append("species")
		tax_dict["tax_id"].insert(0,'nohit')

		out=open(queryID+"-results.txt","w")
		outstring="\nThe sample %s contained %i valid query sequences:\n" % (queryID, querycount[queryID])
		print outstring
		out.write(outstring+"\n")
		
		for tax_rank in reversed(tax_dict["tax_id"]):
#			print tax_rank
			if taxonomy_count.has_key(tax_rank):
				total_per_rank_count=int(0)
				output=[]
				for hit in sorted(taxonomy_count[tax_rank].keys()):
#					print "\t%s: %i" % (hit, len(taxonomy_count[tax_rank][hit])) #This is the count of unique clusters that were assigned
					current_count=int(0)
					current_reads=[]
					for read in taxonomy_count[tax_rank][hit]:
#						print "centroid: %s" % read
#						print "read count in cluster: %i" % cluster_counts[read]
						current_count+=cluster_counts[read]	#sum up the number of actual reads in the clusters assigned to this taxon
						if args.extract_centroid_reads:
							current_reads.append(read)

						elif args.extract_all_reads:
							current_reads.extend(cluster_reads[read])

						total_per_rank_count+=cluster_counts[read]
					output.append("\t%s: %i (%.2f %%)" % (hit, current_count, 100*float(current_count)/querycount[queryID]))
#					print "\t%s: %i (%.2f %%)" % (hit, current_count, 100*float(current_count)/querycount[queryID])
					
					if current_reads: #This list is only non empty if either -e or -E was specified
#						print current_reads
						fas_out = open(hit.replace(" ", "_")+'.fasta',"w")
						for ID in current_reads:
							if args.extract_centroid_reads:
								unknown_seqs_dict[ID].id += "|%i" % cluster_counts[unknown_seqs_dict[ID].id] 
								unknown_seqs_dict[ID].description = unknown_seqs_dict[ID].id
#							print unknown_seqs_dict[ID]
							fas_out.write(unknown_seqs_dict[ID].format("fasta"))
#							fas_out.write(">%s|%i\n%s\n" % ( unknown_seqs_dict[ID].id, cluster_counts[unknown_seqs_dict[ID].id], unknown_seqs_dict[ID].seq ))
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
		os.chdir("../")
		del tax_dict["tax_id"][0] #remove the first element, i.e. 'nohit'
		del tax_dict["tax_id"][-1] #remove the last element, i.e. 'species'
		
		print '\n'+time.strftime("%c")+'\n'
print "Allow to optionally retrieve all read ids (or even reads) that were assigned to a given taxononomic level, i.e. Salmonidae - need to make a dictionary with cetroid ids as keys and a list of reads as values"
