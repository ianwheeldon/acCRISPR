import csv
import numpy as np
import random
import matplotlib.pyplot as plt
import scipy.stats
import statistics
import scipy.integrate
import math
import statsmodels.stats.multitest

def std_norm(x):
    """ Defining the standard normal distribution function """
    F = (1/np.sqrt(2*math.pi))*np.exp(-pow(x,2)/2)
    return F

def p_value_calc(corr_CS_FS_file, FS_gene_file, guides_per_gene, n_noness, significance, prefix):
    """ Calculating p-values for genes """
    
    f1 = corr_CS_FS_file # File containing sgRNA ID, gene name, CS & FS
    sgRNA = []; sg_gene = []; FS_guide = [] # Lists for storing sgRNA ID, gene name & FS of sgRNA
    with open(f1,'r') as f:
        reader = csv.reader(f, delimiter = '\t')
        for row in reader:
            sgRNA.append(row[0])
            sg_gene.append(row[1])
            FS_guide.append(row[3])
    sgRNA.remove(sgRNA[0]); sg_gene.remove(sg_gene[0]); FS_guide.remove(FS_guide[0])
    
    f2 = FS_gene_file # File containing gene names and associated FS
    gene = []; FS_gene = [] # Lists for storing gene names and fitness scores of genes
    with open(f2,'r') as f:
        reader = csv.reader(f, delimiter = '\t')
        for row in reader:
            gene.append(row[0])
            FS_gene.append(row[1])
    gene.remove(gene[0]); FS_gene.remove(FS_gene[0])
    FS_gene = list(map(float, FS_gene))
    
    coverage = len(sgRNA)/len(gene) # Average coverage of the CS-corrected library
    
    mu_null = statistics.median(FS_gene) # The mean of the null distribution is equal to the median of FS_genes
    non_ess = []; # List for storing the set of putative non-essential genes
    
    if significance == '1-tailed':
        for j in range(len(gene)):
            if FS_gene[j] > np.percentile(FS_gene,20): # Genes having FS > 20th percentile are considered as putatively non-essential
                non_ess.append(gene[j])
                
    elif significance == '2-tailed':
        for j in range(len(gene)):
            if FS_gene[j] > np.percentile(FS_gene,2.5) and FS_gene[j] < np.percentile(FS_gene,97.5):
        # Genes having (2.5th percentile < FS < 97.5th percentile) are considered as putatively non-essential
                non_ess.append(gene[j])
        
    sd = []; # List for storing S.D. of FS_pseudogenes from each iteration
    n_pseudo = []; # List for storing no. of sgRNA per pseudogene in each iteration
    
    for k in range(50):
        random_noness = random.sample(non_ess,n_noness) # Random sampling of putative non-essential genes
        guide_pool = {} # Dictionary that stores guides targeting non-essential genes as keys and their FS as values
        for i in random_noness:
            for j in range(len(sgRNA)):
                if i in sg_gene[j]: guide_pool[sgRNA[j]] = FS_guide[j]
                    
        FS_pseudo = []; # List that will store FS values of 1000 pseudogenes
        # Each guide in the pool should be picked no more than twice on average. This
        # means that the max. ratio of total number of guides picked to the pool size
        # should be 2:1. For n guides/pseudogene (i.e. 1000*n guides for 1000 pseudogenes),
        # the pool size has to be at least (1000*n/2) (or 2*pool_size >= 1000*n).
        # No. of guides to be picked per pseudogene (n):
        if 2*len(guide_pool) >= 1000*guides_per_gene: n = guides_per_gene
        # If the pool size is smaller than (1000*n/2), we can pick no more than 2*pool_size
        # guides in total. In that case, no. of guides per pseudogene would be (2*pool_size)/1000
        else: n = int(2*len(guide_pool)/1000)
        
        n_pseudo.append(n)
        
        for i in range(1000): # Generating 1000 pseudogenes
            c = 0
            for j in range(n): # 'n' guides targeting each pseudogene
                x = random.choice(list(guide_pool)) # picking a guide at random
                fs = guide_pool[x] # Finding FS of that guide
                c += float(fs)
            FS_pseudo.append(c/n)
        mu, sigma = scipy.stats.norm.fit(FS_pseudo) # Fitting the FS data to a normal distribution
        sd.append(sigma)
            
    sd_null = np.mean(sd) # S.D. of null dist. = Avg. value of sigma from 50 iterations
    x = np.linspace(1.5*min(FS_gene),1.5*max(FS_gene),1000)
    fitted_data = scipy.stats.norm.pdf(x, mu_null, sd_null)
        
    # Plotting the null distribution
    plt.figure()
    plt.plot(x,fitted_data*1000,'#0E5474')
    plt.title('Mean = %.3f,  S.D. = %.3f, %d guides/pseudogene' % (mu_null, sd_null, round(sum(n_pseudo)/50)),fontsize = 14)
    plt.xlabel('FS_pseudogene', fontsize = 14)
    plt.ylabel('Frequency', fontsize = 14)
    plt.savefig(prefix + '_Null distribution.png', dpi = 300)
    
    # Performing a z-test of significance for every gene in the Cas9 library
    z = []; p_val = [] # Lists for storing z-scores & p-values of gene FS
    if significance == '1-tailed':
        for i in range(len(FS_gene)):
            z.append((FS_gene[i]-mu_null)/sd_null) # Calculation of z-score
            p_val.append(scipy.integrate.quad(std_norm,-np.inf,z[i])[0]) # Calculation of p-value
        
        p_corr = statsmodels.stats.multitest.fdrcorrection(p_val, alpha = 0.05, method = 'indep', is_sorted = False)
        # 'p_corr' is a tuple consisting of 2 arrays. The second array contains corrected p-values of genes
        
    elif significance == '2-tailed':
        for i in range(len(FS_gene)):
            z.append((FS_gene[i]-mu_null)/sd_null) # Calculation of z-score
            p_val.append(1-abs(scipy.integrate.quad(std_norm,-z[i],z[i])[0])) # Calculation of p-value
            
        p_corr = statsmodels.stats.multitest.fdrcorrection(p_val, alpha = 0.05, method = 'indep', is_sorted = False)
        # 'p_corr' is a tuple consisting of 2 arrays. The second array contains corrected p-values of genes
    
    g = prefix + '_p-values_FS_gene.tab'
    with open(g,'w',newline = '') as f:
        writer = csv.writer(f, delimiter = '\t')
        writer.writerow(['gene','FS','Raw p-value','Corrected p-value'])
        for i in range(len(gene)):
            row = []
            row.append(gene[i])
            row.append(float(FS_gene[i]))
            row.append(float(p_val[i]))
            row.append(float(p_corr[1][i]))
            writer.writerow(row)
    
    return coverage, p_corr[1]
