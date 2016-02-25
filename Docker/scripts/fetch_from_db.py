#! /usr/bin/python

import sys, warnings
import os
from Bio import Entrez
import argparse
import shlex, subprocess
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
import urllib
import time
from Bio.Alphabet import generic_dna

### define variables
VERSION = '0.3'
Entrez.email = ""
date = time.strftime("%d-%b-%Y").upper()
geo_syn = {'Europe':["Albania","Andorra","Armenia","Austria","Azerbaijan","Belarus","Belgium","Bosnia and Herzegovina","Bulgaria","Croatia","Cyprus","Czech Republic","Denmark","Estonia","Finland","France","Georgia","Germany","Greece","Hungary","Iceland","Ireland","Italy","Kazakhstan","Kosovo","Latvia","Liechtenstein","Lithuania","Luxembourg","Macedonia","Malta","Moldova","Monaco","Montenegro","Netherlands","Norway","Poland","Portugal","Romania","Russia","San Marino","Serbia","Slovakia","Slovenia","Spain","Sweden","Switzerland","Turkey","Ukraine","United Kingdom"]}
marker_syn = {'COI':['MT-CO1','cox1','COX1','COI','coi','cytochrome oxidase 1','cytochrome oxidase I',
             'cytochrome oxidase subunit 1','cytochrome oxidase subunit I',
             'cox I','coxI', 'cytochrome c oxidase subunit I'],
             'COII':['MT-CO2','cox2','COX2','COII','coii','cytochrome oxidase 2','cytochrome oxidase II',
             'cytochrome oxidase subunit 2','cytochrome oxidase subunit II','cox II',
             'cytochrome c oxidase subunit II'],
             'COIII':['MT-CO3','cox3','COX3','COIII','coiii','cytochrome oxidase 3',
             'cytochrome oxidase III','cytochrome oxidase subunit 3','cytochrome oxidase subunit III'
             ,'cox III','cytochrome c oxidase subunit III'],
             '18S':['18s','18S','SSU rRNA','18S ribosomal RNA','small subunit 18S ribosomal RNA', '18S rRNA'],
             '28S':['28s','28S','LSU rRNA','28S ribosomal RNA','28S large subunit ribosomal RNA'],
             'ATP6':['MT-ATP6','atp6','ATP6','atpase6','atpase 6','ATPase6','ATPASE6'],
             'ATP8':['MT-ATP8','atp8','ATP8','atpase8','atpase 8','ATPase8','ATPASE8'],
             'ATP9':['ATP9','atp9'],
             'CYB':['CYB','MT-CYB','cytochrome-b','cytb','CYTB','Cytb', 'cytochrome B','cytochrome b'],
             'ef1a':['ef1a','Ef1a', 'elongation factor 1-alpha','elongation factor-1 alpha',
             'ef1-alpha','elongation factor 1 alpha'],
             'ND1':['MT-ND1','nd1','ND1','nadh1','NADH1','nadh 1','NADH 1',
             'NADH dehydrogenase subunit 1','NADH dehydrogenase subunit I','nad1','nadh1'],
             'ND2':['MT-ND2','nd2','ND2','nadh2','NADH2','nadh 2','NADH 2','NADH dehydrogenase subunit 2', 
             'NADH dehydrogenase subunit II','nad2','nadh2'],
             'ND3':['MT-ND3','nd3','ND3','nadh3','NADH3','nadh 3','NADH 3','NADH dehydrogenase subunit 3',
             'NADH dehydrogenase subunit III','nad3','nadh3'],
             'ND4':['MT-ND4','nd4','ND4','nadh4','NADH4','nadh 4','NADH 4','NADH dehydrogenase subunit 4',
             'NADH dehydrogenase subunit IV','nad4','nadh4'],
             'ND4L':['MT-ND4L','nd4l','ND4L','nad4l','nadh4l', 'nadh4L', 'NADH dehydrogenase subunit 4l',
             'NADH dehydrogenase subunit 4L'],
             'ND5':['MT-ND5','nd5','ND5','nad5', 'nadh5','NADH5','nadh 5','NADH 5','NADH dehydrogenase subunit 5',
             'NADH dehydrogenase subunit V'],
             'ND6':['MT-ND6','nd6','ND6','nad6','nadh6','NADH6','nadh 6','NADH 6','NADH dehydrogenase subunit 6',
             'NADH dehydrogenase subunit VI'],
             '12S':['12S','rrnS','12srRNA','12s rRNA','12S ribosomal RNA','12S rRNA'],
             '16S':['16S','rrnL','16srRNA','16s rRNA','16S ribosomal RNA','16S rRNA'],
             'C_mos':['C_mos', 'C-mos','c-mos','C-MOS'],
             'GAPDH':['GAPDH','gapdh'],
             'RNApol2':['RNApol2','RNApolII','RNA polymerase II']}
non_gene = ['12S', '16S', '18S', '28S']
genes=[]
gene_search_term = ''
taxids = []
#search_taxa = {}
record_dict = {}
batch_size=1000
BOLD_records = []
gb_accessions = []
####################
#define functions
def get_accession(record):
	"""Given a SeqRecord get the accession"""
	accession = record.name
	return accession
	
def check_email(mail):
        """
        The function checks that you provide an email address to Entrez
        """
	
	from datetime import datetime

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



##########
parser = argparse.ArgumentParser(description='Fetch all available DNA sequences for a list of taxa and a given marker gene', prog='fetch_from_db.py')

parser.add_argument("-t", "--taxlist", help="text file containing list of taxa to search for on Genbank", metavar="<FILE>", action="store")
parser.add_argument("-m", "--marker", help="name of gene to be searched for (put in \"\" if the search term consists of more than one word", metavar="<string>", action="store")
parser.add_argument("-G","--Genbank", help="search Genbank", action="store_true")
parser.add_argument("-B","--BOLD", help="search BOLD", action="store_true")
parser.add_argument("--geo", help="limit BOLD search to countries/continents. for more than one write as follows: \"Austria|UK\"", metavar="<string>", action="store")
parser.add_argument("-o","--out", help="prefix for output files (default=\"out\")", metavar="<string>", action="store", default="out")
parser.add_argument("-@", "--email", help='provide your email address for identification to NCBI', metavar='<email-address>', action="store", default="")

parser.add_argument("--version", action="version", version='%(prog)s v.'+VERSION)
args = parser.parse_args()


if len(sys.argv) < 2:   #if the script is called without any arguments display the usage
	parser.print_usage()
	sys.exit(1)

if args.email:
	Entrez.email = args.email
Entrez.email = check_email(mail = Entrez.email)

if not args.Genbank and not args.BOLD:
	print "\nPlease specify a database to be searched\n"
	parser.print_usage()
	sys.exit(1)

if args.taxlist:
	fh = open(args.taxlist,'r')
	taxa = sorted(list(set([line.strip() for line in open(args.taxlist)])))	#the list(set()) around the list comprehension removes duplicates
#	print taxa

if args.Genbank:
	print "\nQUERYING GENBANK\n"
	if not args.marker:
		print "no marker specified - will fetch all available sequences for the desired taxa"
	else:
		print "check for synonyms for \"%s\" (this is relevant only for Genbank searches)" %args.marker
		for g in marker_syn.keys():
#			print "checking %s" %g
			if args.marker in marker_syn[g]:
#				print "found: %s: %s" %(g, marker_syn[g])
				if g in non_gene:
					gene_search_term = " AND ((%s))" %") OR (".join(marker_syn[g])
				else:
					gene_search_term = " AND ((%s[gene]))" %"[gene]) OR (".join(marker_syn[g])
				break

		if not gene_search_term:
			gene_search_term = " AND (%s[gene])" %args.marker

	print "\nfetching accessions ..\n"
	for tax in taxa:
		search = "("+tax+"[orgn])"+gene_search_term
#		print search
		handle = Entrez.esearch(db='nuccore', term=search, retmax=10000)
		records = Entrez.read(handle)
#		print records.keys()
		if len(records['IdList']) < records['Count']: #['IdList'] contains the number of records downloaded, while ['Count'] holds the number of actually existing records
			handle = Entrez.esearch(db='nuccore', term=search, retmax=records['Count'])
			records = Entrez.read(handle)
		print "%s\t%s" %(tax, records['Count'])
#		search_taxa[tax]=records['Count']
		taxids.extend(records['IdList'])
#		print "total: %s" %len(taxids)
#		search_taxa[tax]={'Genbank': records['IdList']}

#	print "\ntotal: %s\n" %len(taxids)

	taxids = list(set(taxids))	
#	print "\ntotal number of accessions fetched: %s\n" %len(list(set(taxids)))
	print "\ntotal number of accessions fetched: %s\n" %len(taxids)
#	print search_taxa

	#now download the records
	print "\ndownloading the records .. processing %s accessions per batch\n" %batch_size
#	for start in range(0,len(list(set(taxids))),batch_size):
	for start in range(0,len(taxids),batch_size):
#		end = min(len(list(set(taxids))), start+batch_size)
		end = min(len(taxids), start+batch_size)
#		handle = Entrez.efetch(db='nuccore', id=sorted(list(set(taxids)))[start:end], rettype='gb',retmax=batch_size)
		handle = Entrez.efetch(db='nuccore', id=taxids[start:end], rettype='gb',retmax=batch_size)
		recs=SeqIO.parse(handle,'gb')
		for rec in recs:
			if not rec.id in record_dict.keys():
				record_dict[rec.id]=rec
		print "dowloaded\t%i unique records" %len(record_dict)	

if args.BOLD:
	print "\nQUERYING BOLD - atm experimental\n"
	if args.geo:
		if geo_syn.has_key(args.geo):
			geo_search ='&geo='+"|".join(geo_syn[args.geo])
		else:
			geo_search = '&geo='+args.geo
	else:
		geo_search = ''
	for tax in taxa:
		url = "http://www.boldsystems.org/index.php/API_Public/sequence?taxon=%s%s" %(tax, geo_search)
#		url = "http://www.boldsystems.org/index.php/API_Public/specimen?taxon=%s%s&format=tsv" %(tax, geo_search)
#		print url
		urllib.urlretrieve(url, filename='temp.fasta')
		temp_BOLD_records = list(SeqIO.parse(open('temp.fasta','r'),'fasta'))
		print "%s\t%s" %(tax, len(temp_BOLD_records))
		BOLD_records.extend(temp_BOLD_records)
		os.remove('temp.fasta')

	print "\nquerying NCBI for taxonomic ids (taxids) associated with the BOLD records\n"
	for rec in BOLD_records:
#		print "%s (%s)" %(rec.description, len(rec.description.split("|")))
		if len(rec.description.split("|"))==4: #if four elements exist then there is a valid genbank accession.
#			print "record has genbank accession: %s" %rec.description.split("|")[3]
			gb_accessions.append(rec.description.split("|")[3]) #collect the genbank accessions and download those records directly from genbank below
		else: #these seqeuences do not exist on Genbank, so I will create a seqrecord for those
			rec.id = rec.description.split("|")[0]
			rec.name = rec.id
			rec.seq = Seq(str(rec.seq).replace("-",""), generic_dna) #remove any gaps from the sequence
			rec.features.append(SeqFeature(FeatureLocation(0,len(rec.seq),strand=1), type="source"))        #create source feature
			rec.features[0].qualifiers['original_desc'] = rec.description
			rec.features[0].qualifiers['organism'] = [rec.description.split("|")[1]]
			rec.features[0].qualifiers['comment'] = ['downloaded from BOLD']
			rec.annotations['date'] = date
		
			if len(rec.features[0].qualifiers['organism'][0].split(" ")) > 1:
#				print "check if a taxid exists for the taxon %s\n" %rec.features[0].qualifiers['organism'][0]
				handle = Entrez.esearch(db="Taxonomy", term=rec.features[0].qualifiers['organism'][0])
				taxon = Entrez.read(handle)
				if taxon['IdList']: #if the search yielded a result
					rec.features[0].qualifiers['db_xref'] = taxon['IdList'][0]
#					print "found %s" %taxon['IdList']

				else:
#					print "check if a taxid exists for the taxon %s\n" %rec.features[0].qualifiers['organism'][0].split(" ")[0]
					handle = Entrez.esearch(db="Taxonomy", term=rec.features[0].qualifiers['organism'][0].split(" ")[0])
					taxon = Entrez.read(handle)
					if taxon['IdList']:
#						print "found %s" %taxon['IdList']
						pass
					else:
						print "NO VALID TAXID FOUND - WHAT TO DO NOW? %s" %rec.features[0].qualifiers['organism'][0]		

			elif len(rec.features[0].qualifiers['organism'][0].split(" ")) == 1:
				rec.features[0].qualifiers['organism'] = [rec.features[0].qualifiers['organism'][0]+" sp."]
#				print "check if a taxid exists for the taxon %s\n" %rec.features[0].qualifiers['organism'][0].split(" ")[0]
				handle = Entrez.esearch(db="Taxonomy", term=rec.features[0].qualifiers['organism'][0].split(" ")[0])
				taxon = Entrez.read(handle)
				if taxon['IdList']:
					pass
#					print "found %s" %taxon['IdList']
				else:
					print "NO VALID TAXID FOUND - WHAT TO DO NOW? %s" %rec.features[0].qualifiers['organism'][0]		
			
#			print rec
#			print rec.features[0]
			record_dict[rec.id] = rec	#add record to final dictionary

	print "total: %s" %len(BOLD_records)
	gb_accessions = list(set(gb_accessions))
	print "with genbank accession: %s" %len(gb_accessions)

	for i in reversed(range(len(gb_accessions))):	#remove any accessions from the list that are already in the dictionary
		if record_dict.has_key(gb_accessions[i]):
			del gb_accession[i]

	print "\ndownloading %s novel records from Genbank.. %s records per batch\n" %(len(list(set(gb_accessions))), batch_size)
	

	for start in range(0,len(gb_accessions),batch_size):
#		end = min(len(list(set(taxids))), start+batch_size)
		end = min(len(taxids), start+batch_size)
#		handle = Entrez.efetch(db='nuccore', id=sorted(list(set(taxids)))[start:end], rettype='gb',retmax=batch_size)
		handle = Entrez.efetch(db='nuccore', id=gb_accesssions[start:end], rettype='gb',retmax=batch_size)
		recs=SeqIO.parse(handle,'gb')
		for rec in recs:
			if not rec.id in record_dict.keys():
				record_dict[rec.id]=rec
		print "dowloaded\t%i unique records" %len(record_dict)	


#print record_dict
print "\nwriting the records to file: %s.gb\n" %args.out
output_handle = open(args.out+'.gb','w')
for key in record_dict.keys():
	SeqIO.write(record_dict[key], output_handle, "gb")

output_handle.close()
	
#wget -O- "http://www.boldsystems.org/index.php/API_Public/sequence?taxon=Chordata&geo=Florida|Germany|France" > test.fasta

#http://www.boldsystems.org/index.php/API_Public/specimen?taxon=Aleiodes%20rossicus&format=tsv
