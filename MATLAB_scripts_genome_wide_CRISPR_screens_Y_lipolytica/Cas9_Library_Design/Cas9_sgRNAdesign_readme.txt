crispr_analyzer ・Determines type of comparison (amino acid or nucleotide), call genome using gene_import・ and analyzes sequences and designs sgRNAs with gene_analysis・

Gene_Import ・Imports fasta file of genes (mRNA-encoding DNA sequence) and truncates genes for sgRNA design

Gene_Analysis ・Imports gene and reference chromosome, divides sequence into 30mers, calculates on target cutting efficiency, ranks scores per gene, and BLASTs against genome to ensure uniqueness

This was the main, overarching file I used.  First, it made you say whether you were conducting an amino or nucleotide test, and then called the genome based on what you input using the gene_import or amino_import functions.  Once the files were imported, it analyzed each imported file using the Gene_analysis or gene_analysis_amino functions, again, depending on what was selected earlier.  These functions output the string sequences and corresponding blast scores for the sgRNA sequences that met all criteria.  It then took those values and output them into excel using the sgRNA_to_Excel function.  

Files used in gene_Analysis・・divider_30mer, on_target_score_calculator, score_ranking, BLAST_Redux, FormatOutput, and sgRNA_to_Excel

Cory Schwartz
PhD Candidate, Wheeldon Lab
UC Riverside
