import csv
from collections import defaultdict
import numpy as np
import math

def counts_to_log2fc(count_file, rep_file, prefix):
    """ Calculating CS & FS of guides using raw read counts as input """
    f1 = rep_file # File containing information on replicates
    col_1 = []; col_2 = []; # Lists for storing information in columns 1 & 2 of the replicate file
    with open(f1,'r') as f:
        reader = csv.reader(f, delimiter = '\t')
        for row in reader:
            col_1.append(row[0])
            col_2.append(row[1])
    col_1.remove(col_1[0]); col_2.remove(col_2[0])
        
    # If control is common between the 2 treatments:
    if len(set(col_1)) == 3:
        replicate_map = {'Control':[], 'Treatment_CS':[], 'Treatment_FS':[]} # Dictionary for storing replicate information for control & treatment conditions
        for i in range(len(col_1)):
            if 'Control' in col_1[i]: replicate_map['Control'].append(col_2[i])
            elif 'Treatment_CS' in col_1[i]: replicate_map['Treatment_CS'].append(col_2[i])
            else: replicate_map['Treatment_FS'].append(col_2[i])
                
        n_ctrl_rep = len(replicate_map['Control']) # No. of replicates for control condition
        n_CS_rep = len(replicate_map['Treatment_CS']) # No. of replicates for CS treatment condition
        n_FS_rep = len(replicate_map['Treatment_FS']) # No. of replicates for FS treatment condition
        n_samples = n_ctrl_rep + n_CS_rep + n_FS_rep; # Total number of samples
            
        counts = defaultdict(list) # Dictionary for storing sample names as keys and list of raw sgRNA counts as values
            
    # If control is different for CS & FS:
    elif len(set(col_1)) == 4:
        replicate_map = {'Control_CS':[], 'Control_FS':[], 'Treatment_CS':[], 'Treatment_FS':[]} # Dictionary for storing replicate information for control & treatment conditions
        for i in range(len(col_1)):
            if 'Control_CS' in col_1[i]: replicate_map['Control_CS'].append(col_2[i])
            elif 'Control_FS' in col_1[i]: replicate_map['Control_FS'].append(col_2[i])
            elif 'Treatment_CS' in col_1[i]: replicate_map['Treatment_CS'].append(col_2[i])
            else: replicate_map['Treatment_FS'].append(col_2[i])
                
        n_ctrl_CS_rep = len(replicate_map['Control_CS']) # No. of replicates for CS control condition
        n_ctrl_FS_rep = len(replicate_map['Control_FS']) # No. of replicates for FS control condition
        n_CS_rep = len(replicate_map['Treatment_CS']) # No. of replicates for CS treatment condition
        n_FS_rep = len(replicate_map['Treatment_FS']) # No. of replicates for FS treatment condition
        n_samples = n_ctrl_CS_rep + n_ctrl_FS_rep + n_CS_rep + n_FS_rep; # Total number of samples
            
        counts = defaultdict(list) # Dictionary for storing sample names as keys and list of raw sgRNA counts as values
        
    f2 = count_file # File containing sgRNA ID, associated gene name & raw counts for each sample
    sgRNA = []; gene = []; # Lists storing sgRNA IDs and associated gene names
    j = 0
    with open(f2,'r') as f:
        reader = csv.reader(f, delimiter = '\t')
        for row in reader:
            j += 1
            if j > 1: # The first line is ignored since the column titles do NOT have to be included in the lists
                sgRNA.append(row[0])
                gene.append(row[1])
                for i in range(2,n_samples+2): counts[col_2[i-2]].append(int(row[i])+1) # Storing (count + 1) of each guide for each sample in the dictionary
            
    for k in counts:
        s = sum(counts[k])*pow(10,-6) # Calculating total of sgRNA counts per million for each replicate
        for i in range(len(counts[k])): counts[k][i] = counts[k][i]/s # Total normalizing counts
        
    # If control is common between the 2 treatments:
    if len(set(col_1)) == 3:
        # Defining lists for storing avg. abundances (across replicates) of each guide in control & treatment conditions
        ctrl_abundance = list(np.zeros(len(sgRNA)));
        tr_CS_abundance = list(np.zeros(len(sgRNA)));
        tr_FS_abundance = list(np.zeros(len(sgRNA)));
        
        for k in counts:
            if k in replicate_map['Control']: # If a replicate belongs to control, the counts for that replicate get added to the total control abundance of each sgRNA
                for i in range(len(counts[k])): ctrl_abundance[i] += counts[k][i]
            elif k in replicate_map['Treatment_CS']: # If a replicate belongs to CS treatment, the counts for that replicate get added to the total CS treatment abundance of each sgRNA
                for i in range(len(counts[k])): tr_CS_abundance[i] += counts[k][i]
            else: # If a replicate belongs to FS treatment, the counts for that replicate get added to the total FS treatment abundance of each sgRNA
                for i in range(len(counts[k])): tr_FS_abundance[i] += counts[k][i]
                
        CS = []; FS = []; # Lists for storing CS & FS of each sgRNA
        for i in range(len(ctrl_abundance)):
            # Converting total abundance to average abundance (by dividing total abundances by no. of replicates) for all controls and treatments
            ctrl_abundance[i] = ctrl_abundance[i]/n_ctrl_rep
            tr_CS_abundance[i] = tr_CS_abundance[i]/n_CS_rep
            tr_FS_abundance[i] = tr_FS_abundance[i]/n_FS_rep
            # Calculation of CS & FS of each sgRNA & storing the values in respective lists
            CS.append(math.log2(ctrl_abundance[i]/tr_CS_abundance[i]))
            FS.append(math.log2(tr_FS_abundance[i]/ctrl_abundance[i]))
            
    # If control is different for CS & FS:
    elif len(set(col_1)) == 4:
        # Defining lists for storing avg. abundances (across replicates) of each guide in control & treatment conditions
        ctrl_CS_abundance = list(np.zeros(len(sgRNA)));
        ctrl_FS_abundance = list(np.zeros(len(sgRNA)));
        tr_CS_abundance = list(np.zeros(len(sgRNA)));
        tr_FS_abundance = list(np.zeros(len(sgRNA)));
        
        for k in counts:
            if k in replicate_map['Control_CS']: # If a replicate belongs to CS control, the counts for that replicate get added to the total CS control abundance of each sgRNA
                for i in range(len(counts[k])): ctrl_CS_abundance[i] += counts[k][i]
            elif k in replicate_map['Control_FS']: # If a replicate belongs to FS control, the counts for that replicate get added to the total FS control abundance of each sgRNA
                for i in range(len(counts[k])): ctrl_FS_abundance[i] += counts[k][i]
            elif k in replicate_map['Treatment_CS']: # If a replicate belongs to CS treatment, the counts for that replicate get added to the total CS treatment abundance of each sgRNA
                for i in range(len(counts[k])): tr_CS_abundance[i] += counts[k][i]
            else: # If a replicate belongs to FS treatment, the counts for that replicate get added to the total FS treatment abundance of each sgRNA
                for i in range(len(counts[k])): tr_FS_abundance[i] += counts[k][i]
                
        CS = []; FS = []; # Lists for storing CS & FS of each sgRNA
        for i in range(len(ctrl_CS_abundance)):
            # Converting total abundance to average abundance (by dividing total abundances by no. of replicates) for all controls and treatments
            ctrl_CS_abundance[i] = ctrl_CS_abundance[i]/n_ctrl_CS_rep
            ctrl_FS_abundance[i] = ctrl_FS_abundance[i]/n_ctrl_FS_rep
            tr_CS_abundance[i] = tr_CS_abundance[i]/n_CS_rep
            tr_FS_abundance[i] = tr_FS_abundance[i]/n_FS_rep
            # Calculation of CS & FS of each sgRNA & storing the values in respective lists
            CS.append(math.log2(ctrl_CS_abundance[i]/tr_CS_abundance[i]))
            FS.append(math.log2(tr_FS_abundance[i]/ctrl_FS_abundance[i]))
            
    # Writing results in a new .TAB file
    g = prefix + '_guide_CS_FS.tab'
    with open(g,'w',newline = '') as f:
        writer = csv.writer(f, delimiter = '\t')
        writer.writerow(['Guide ID','gene','CS','FS'])
        for i in range(len(sgRNA)):
            row = []
            row.append(sgRNA[i])
            row.append(gene[i])
            row.append(float(CS[i]))
            row.append(float(FS[i]))
            writer.writerow(row)
