#! /usr/bin/python

VERSION="0.96.2-global"

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

    if not target in samples:
        sys.exit("The taxon you are looking for '%s' was not detected" %target)

    positives = BIOM.norm(axis='sample').data(target,'observation')

    index=0
    for o in positives:
#        print o
        if str(o) != '0.0':
            print "%s\t(%.4f %%)" %(samples[index],o)
        index+=1



