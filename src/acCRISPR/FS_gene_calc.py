import csv

def FS_gene(corr_CS_FS_file, prefix):
    """ Calculating fitness scores of genes using sgRNA FS data """
    f1 = corr_CS_FS_file
    guide = []; gene_name = []; FS_guide = [] # Lists containing all guides in the library, associated gene names and FS
    with open(f1,'r') as f:
        reader = csv.reader(f, delimiter = '\t')
        for row in reader:
            guide.append(row[0])
            gene_name.append(row[1])
            FS_guide.append(row[3])
    guide.remove(guide[0]); gene_name.remove(gene_name[0]); FS_guide.remove(FS_guide[0])
    gene = list(set(gene_name)) # List containing all genes covered by the GW library
    FS_gene = [] # List storing fitness scores of genes
    for i in gene:
        total = 0; c = 0 # variable c counts the no. of guides targeting gene i
        for j in range(len(guide)):
            if i in gene_name[j]:
                total += float(FS_guide[j])
                c += 1
        FS_gene.append(total/c) # FS of a gene is the arithmetic mean of FS of all sgRNA targeting that gene
    
    # Writing results in a new .TAB file
    g = prefix + '_FS_genes.tab'
    with open(g,'w',newline = '') as f:
        writer = csv.writer(f, delimiter = '\t')
        writer.writerow(['gene','FS_gene'])
        for i in range(len(gene)):
            row = []
            row.append(gene[i])
            row.append(float(FS_gene[i]))
            writer.writerow(row)
