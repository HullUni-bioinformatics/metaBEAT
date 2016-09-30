#! /usr/bin/python

VERSION="0.97.3-global"

def vsearch_cluster_full_length(infile, cluster_match, threads, sampleID, verbose=0):
    """
    The function runs vsearch
    """
    import shlex, subprocess

    cmd = "vsearch --cluster_fast %s --id %.2f --strand both --threads %s --centroids %s_centroids.fasta --uc %s.uc" % (infile, cluster_match, threads, sampleID, sampleID)
    if verbose:
        print cmd
    cmdlist = shlex.split(cmd)
    stdout, stderr = subprocess.Popen(cmdlist, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate() # , stdout=subprocess.PIPE).communicate()
    if verbose:
        if stdout:
            print stdout
        if stderr:
            print stderr

def find_most_abundant_seq_from_uc(uc):
    """
    The function finds most abundant haplotype
    """
    
    centroids={}
    
    uc_in=open(uc, 'r')
    for line in uc_in:
        if line.startswith('C'):
            if not line.split("\t")[2] in centroids:
                centroids[int(line.split("\t")[2])]=[]
            centroids[int(line.split("\t")[2])].append(line.split("\t")[8])
            
    return centroids[sorted(centroids)[-1]][0]



def find_target_OTUs_by_taxonomy(BIOM, target, level='all'):

    """
    Find all OTUs that have been assigned to target taxon
    """
    from biom.table import Table
    
    return_OTUs = {}
    levels = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    levels_valid = ['all','unassigned']
    levels_valid.extend(levels)
    search_levels = []
    
    print "SEARCH TERM: '%s'" %target
    
    #check if the first OTU has 'taxonomy' metadata attached, if yes assume all others have too and resume
    if not 'taxonomy' in BIOM.metadata(axis='observation')[0]:
        raise KeyError('The BIOM table you are trying to screen does not have taxonomy metadata attached to it')
    else:
        print "Found taxonomy metadata with OTUs - ok!"
    
    if not level in levels_valid:
        raise KeyError("The taxonomic level you are trying to search: '%s', does not exist" %level)
    else:
        if level == 'all':
            search_levels = []
        
        elif level == 'unassigned':
            search_levels = [0]
        else:
            search_levels.append(int(levels.index(level)))
            
        print "Screening at taxonomic level: '%s'\n" %level
#        print search_levels
   
    OTUs=BIOM.ids(axis='observation')
    per_OTU_metadata = BIOM.metadata(axis='observation')
    
    
    for i in range(len(per_OTU_metadata)):

        if level == 'all':
            search_levels = range(len(per_OTU_metadata[i]['taxonomy']))
        
        if len(per_OTU_metadata[i]['taxonomy']) > search_levels[-1]:
            for l in search_levels:
                if target in per_OTU_metadata[i]['taxonomy'][l]:
                    return_OTUs[OTUs[i]] = per_OTU_metadata[i]['taxonomy']
                    break
    
    print "\nIdentified %i OTU(s) assigned to '%s'." %(len(return_OTUs), target)
    
    return return_OTUs

def filter_BIOM_by_per_sample_read_prop(BIOM, min_prop=0.01):
    """
    Filter OTU table by mininimum reads per sample
    """

    import numpy as np
    from biom.table import Table

    print "\nFiltering at level: %s %%\n" %(min_prop*100)
    
#    print "input table:\n"
#    print BIOM
#    print "\n"
    
    sample_ids = BIOM.ids(axis='sample')
    observation_ids = BIOM.ids(axis='observation')
    data_to_biom = []
    sample_metadata = BIOM.metadata(axis='sample')
    observation_metadata = BIOM.metadata(axis='observation')
    sums = BIOM.sum(axis='sample')

    for OTU in observation_ids:
        orig=BIOM.data(OTU, axis='observation')
        for i in range(len(orig)):
            if not int(orig[i]) == 0:
                if not int(orig[i]) >= sums[i]*min_prop:
                    orig[i] = '0.0'
        data_to_biom.append(orig)
    
    data = np.asarray(data_to_biom)

    #construct adjusted table
    table = Table(data, observation_ids, sample_ids, table_id='OTU table', sample_metadata=sample_metadata, observation_metadata=observation_metadata)

    #Filter OTUs with sum = '0'
    to_exclude = []
    observation_sums = table.sum(axis='observation')
    for i in range(len(observation_sums)):
        if int(observation_sums[i]) == 0:
            to_exclude.append(observation_ids[i])
    
    print "Removing %i OTUs for lack of support\n" %len(to_exclude)
    table.filter(to_exclude, invert=True, axis='observation',inplace=True)
    
#    print table
    return table

def collapse_biom_by_taxonomy(in_table):
	"""
	collapse OTU table observations by taxonomy informatition
	"""
	
	bin_f = lambda id_, x: x['taxonomy'][-1].split("__")[1]
	biom_out_collapsed = in_table.collapse(bin_f, norm=False, min_group_size=1, axis='observation')

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
	biom_out_collapsed.add_metadata(md=collapsed_meta, axis='observation')
	biom_out_collapsed_sorted = biom_out_collapsed.sort_order(collapsed_order, axis='observation')
#	print biom_out_collapsed_sorted.to_tsv(header_key='taxonomy', header_value='taxomomy')
#	print "\n\n"

	return biom_out_collapsed_sorted

def BIOM_tsv_to_R_transpose(in_tsv, out_csv):
    """
    Parse a biom table in tsv format and transpose it for input into R
    """
    
    from biom import Table
    
    tsv = open(in_tsv)
    #in_tsv = open('COI-trim30min100-merge-c3-id97-OTU-taxonomy.kraken.tsv')
    func = lambda x : x
    intable = Table.from_tsv(tsv,obs_mapping=None, sample_mapping=None, process_func=func)
    outtable = intable.transpose()
    out=open("transposed.tsv","w")
    out.write(outtable.to_tsv(header_key=None, header_value=None))
    out.close()

    #refine
    intable = open('transposed.tsv','r')
    temp = intable.next()

    out=''
    for line in intable:
        if line.startswith('#'):
            if line.strip().endswith('taxomomy'):
                print "Removing taxonomy"
                line = ",".join(line.strip().split("\t")[:-1]).replace('#OTU ID','Sample').replace('\t',',')+'\n'
            line = line.replace('#OTU ID','Sample').replace('\t',',')
            out+=line
        else:
            line = line.replace('\t',',')
            out+=line

    outtable = open(out_csv,'w')
    outtable.write(out)
    outtable.close()

def load_BIOM(table):
    """
    load a BIOM table from BIOM format
    """
    from biom.table import Table
    import json

    with open(table) as data_file:
        data = json.load(data_file)
    t = Table.from_json(data)

    return t


def plot_BIOM_as_heatmap(BIOM, width_in=5, target_png=None):
    """
    plot a simple heatmap from a BIOM table object
    """

    import matplotlib.pyplot as plt
    import pandas as pd
    import numpy as np
#    %matplotlib inline
    #%pylab inline

    column_labels = BIOM.ids('observation')
    row_lables = BIOM.ids('sample')
    data = []

    for o in column_labels:

         data.append(BIOM.data(o,axis='observation'))

    data = np.asarray(data)

    # Plot it out
    fig, ax = plt.subplots()
    heatmap = ax.pcolor(data, cmap=plt.cm.Reds, alpha=0.8)

    # Format
    fig = plt.gcf()
    fig.set_size_inches(width_in,8)

    # turn off the frame
    ax.set_frame_on(False)

    # put the major ticks at the middle of each cell
    ax.set_yticks(np.arange(data.shape[0]) + 0.5, minor=False)
    ax.set_xticks(np.arange(data.shape[1]) + 0.5, minor=False)

    # want a more natural, table-like display
    ax.invert_yaxis()
    ax.xaxis.tick_top()

    ax.set_xticklabels(row_lables, minor=False)
    ax.set_yticklabels(column_labels, minor=False)

    # rotate the
    plt.xticks(rotation=90)

    ax.grid(False)

    # Turn off all the ticks
    ax = plt.gca()

    for t in ax.xaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    for t in ax.yaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False

   #write to file
    if target_png:
        fig.savefig(target_png, dpi=100)

def BIOM_return_by_tax_level(taxlevel, BIOM):
    """
    Returns a new BIOM object containing only OTUS that were assigned at least to a given taxonomic level
    """
    
    from biom.table import Table
    
    return_OTUs = {}
    levels = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'unassigned']
    search_level=''
    
    if not taxlevel in levels:
        raise KeyError("The taxonomic level you are trying to search: '%s', is not valid" %level)
    
    #check if the first OTU has 'taxonomy' metadata attached, if yes assume all others have too and resume
    if not 'taxonomy' in BIOM.metadata(axis='observation')[0]:
        raise KeyError('The BIOM table you are trying to screen does not have taxonomy metadata attached to it')
    else:
        print "Found taxonomy metadata with OTUs - ok!"
    
    if taxlevel == 'unassigned':
        search_level = 0;
    else:
        search_level = int(levels.index(taxlevel))
                           
    OTUs=BIOM.ids(axis='observation')
    per_OTU_metadata = BIOM.metadata(axis='observation')
                           
    for i in range(len(per_OTU_metadata)):
        if len(per_OTU_metadata[i]['taxonomy']) > search_level:
            return_OTUs[OTUs[i]] = per_OTU_metadata[i]['taxonomy']
                           
    outBIOM = BIOM.filter(return_OTUs.keys(), invert=False, axis='observation',inplace=False)
    
    return outBIOM

def BIOM_return_clipped_taxonomy(taxlevel, BIOM):
    """
    Returns a BIOM table for which the taxonomy has been clipped at a certain level
    """
    
    from biom.table import Table
    import numpy as np
    
    return_OTUs = {}
    levels = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'unassigned']
    clip_level=''
    to_drop=[]
    
    if not taxlevel in levels:
        raise KeyError("The taxonomic level you are trying to search: '%s', is not valid" %level)
    
    clip_level = int(levels.index(taxlevel))+1
    #check if the first OTU has 'taxonomy' metadata attached, if yes assume all others have too and resume
    if not 'taxonomy' in BIOM.metadata(axis='observation')[0]:
        raise KeyError('The BIOM table you are trying to screen does not have taxonomy metadata attached to it')
    else:
        print "Found taxonomy metadata with OTUs - ok!"
        
    sample_ids = BIOM.ids(axis='sample')
    observation_ids = BIOM.ids(axis='observation')
    data_to_biom = []
    sample_metadata = BIOM.metadata(axis='sample')
    observation_metadata = BIOM.metadata(axis='observation')
    
    for OTU in observation_ids:
        orig=BIOM.data(OTU, axis='observation')
        data_to_biom.append(orig)
        
    data = np.asarray(data_to_biom)
    
    for i in range(len(observation_metadata)):
        if len(observation_metadata[i]['taxonomy']) > clip_level:
            observation_metadata[i]['taxonomy'] = observation_metadata[i]['taxonomy'][:clip_level]
        if 'unknown' in observation_metadata[i]['taxonomy'][-1]:
            print "fishy: %s" %observation_metadata[i]['taxonomy']
            to_drop.append(observation_ids[i])
#        print observation_metadata[i]['taxonomy']
        
    #construct adjusted table
    outtable = Table(data, observation_ids, sample_ids, table_id='OTU table', sample_metadata=sample_metadata, observation_metadata=observation_metadata)
    
    if to_drop:
        outtable.filter(to_drop, invert=True, axis='observation',inplace=True)
        
    
    return outtable

def per_taxon_abundance(BIOM,tsv):
    """
    Return a table (tsv)"""
    
    #Get observation ids
    obs_ids = BIOM.ids(axis='observation')
    
    #Determine the total number of samples
    total = len(BIOM.ids(axis='sample'))

    
    by_count = {}
    #For each OTU determine the number of samples it is missing from 
    for o in obs_ids:
        count=0
        temp = BIOM.data(o, axis='observation')
        
        for i in reversed(range(len(temp))):
            if temp[i] == 0:
                count+=1
    
        if not by_count.has_key(int(len(temp)-count)):
            by_count[int(len(temp)-count)] = []
        by_count[int(len(temp)-count)].append(o)
    
    #could be sorted by abundance
    out = "#taxon\tsample_count\tsample_proportion\n"
    out_fh = open(tsv,'w')
    for c in sorted(by_count):
        for t in sorted(by_count[c]):
            out+="%s\t%i\t%.2f\n" %(t,c,float(c)/total)

    out_fh.write(out)
    out_fh.close()


def write_BIOM(BIOM, target_file, BIOM_out=True, tsv_out=True, program='metaBEAT v.'+VERSION):
    """
    Write BIOM object to file. Will create two files - 1. in table in BIOM format; 2. table in tsv format
    """
    from biom.table import Table

    if BIOM_out:
        out=open(target_file+".biom","w")
        BIOM.to_json(program, direct_io=out)
        out.close()

    if tsv_out:
            out=open(target_file+".tsv","w")
            out.write(BIOM.to_tsv(header_key='taxonomy', header_value='taxomomy'))
            out.close()

def adjust_metadata_sample_ids(intable, metadata):
    """
    Adjust metadata sample ids to match the sample ids and the order in table
    """
    
    
    #get sample ids in table
    ids = []
    table = open(intable,'r')
    table.next()
    
    for line in table:
        ids.append(line.split(",")[0])
        
    table.close()
#    print ids
    
    #parse metadata
    meta_in = open(metadata,'r')
    out = meta_in.next()

    meta_ids = {}
    meta_data = {}
    for line in meta_in:
        l = line.strip().split(",")
#        print l[0]
        for i in range(len(ids)):
            if ids[i].startswith(l[0]):
                if not l[0] in meta_ids:
                    meta_ids[l[0]] = []
                meta_ids[l[0]].append(ids[i])
                
                if len(meta_ids[l[0]]) == 1:
#                    print "Found matches for %s: %s" %(l[0],meta_ids[l[0]] )
                    l[0] = meta_ids[l[0]][0]
                    meta_data[i] = l
                else:
                    raise KeyError("ambiguous match for %s" %(l[0]))
                
        else:
            "Did not find a match for %s" %l[0]
    
    meta_in.close()
    
    if len(meta_ids) == len(ids):
        to_correct = 1 #always sort 
        for s in meta_ids:
            if not s == meta_ids[s][0]:
                print "Sample ids in '%s' match those in '%s'" %(metadata, intable)
                to_correct+=1
#            else:
#                print "%s is identical to %s" %(s,meta_ids[s][0])

        if to_correct:
            print "ok - updating file %s" %metadata
            #convert
            outfh = open(metadata,'w')
            for s in sorted(meta_data):
                out+=",".join(meta_data[s])+'\n'
        
            outfh.write(out)
            outfh.close()
        
        else:
            print "Nothing to do here - sample ids in '%s' match those in '%s'" %(metadata, intable)
        
        
    else:
        print "Did not find metadata for all samples. "


def find_target(BIOM, target):

    """
    Find all samples and the proportion of reads a taxon was detected in
    """

    from biom.table import Table

    samples = BIOM.ids('sample')

    if not target in BIOM.ids(axis='observation'):
        print "The taxon you are looking for '%s' was not detected" %target
    else:
        positives = BIOM.norm(axis='sample').data(target,'observation')

        index=0
        for o in positives:
            if str(o) != '0.0':
                print "%s\t(%.4f %%)" %(samples[index],o*100)
            index+=1



