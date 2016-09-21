#! /usr/bin/python

VERSION="0.97-global"

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



