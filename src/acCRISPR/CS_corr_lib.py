import csv
from collections import defaultdict
import numpy as np


def generate_CS_corr_lib(CS_FS_file, cutoff, prefix):
    """ Generating a CS-corrected sgRNA library"""
    f1 = CS_FS_file # File containing sgRNA ID, gene name, CS & FS
    sgRNA = []; gene = []; CS = []; FS = [] # Lists for storing sgRNA ID, gene name, CS & FS of sgRNA
    with open(f1,'r') as f:
        reader = csv.reader(f, delimiter = '\t')
        for row in reader:
            sgRNA.append(row[0])
            gene.append(row[1])
            CS.append(row[2])
            FS.append(row[3])
    sgRNA.remove(sgRNA[0]); gene.remove(gene[0]); CS.remove(CS[0]); FS.remove(FS[0])
    guide_details = tuple(zip(sgRNA,gene,CS,FS))
    # guide_details is a tuple, each element of which is a tuple of the form (sgRNA ID, gene, CS, FS) 
    gene_set = list(set(gene)) # List of genes covered by the GW-library
    
    gene_guides = defaultdict(list) # A dictionary for storing gene name as key and a list of attributes (i.e., ID, gene, CS, FS) of all sgRNA targeting that gene as values
    for g in gene_set:
        for i in range(len(guide_details)):
            if g in guide_details[i][1]: gene_guides[g].append(guide_details[i])
            
    val = []; # List for storing details of guides included in the CS-corrected library
        
    for g in gene_guides:
        c = 0; # Counts the number of good cutters (CS > cutoff) targeting a gene
        y = np.array(gene_guides[g])
        for i in range(y.shape[0]):
            if float(y[i,2]) > cutoff:
                c += 1
                val.append(y[i])
        if c == 0: # if no good cutters target a gene, the best of the existing cutters is picked
            p = y[:,2].astype(float) # Converting CS values to float from str
            x = np.where(p == max(p)) # Finding the array position of guide having max. CS (best cutter)
            val.append(y[x[0][0]])
            
    # Writing results in a new .TAB file
    g = prefix + '_corrected_guide_CS_FS.tab'
    with open(g,'w',newline = '') as f:
        writer = csv.writer(f, delimiter = '\t')
        writer.writerow(['Guide ID','gene','CS','FS'])
        for i in val:
            row = []
            row.append(i[0])
            row.append(i[1])
            row.append(float(i[2]))
            row.append(float(i[3]))
            writer.writerow(row)
