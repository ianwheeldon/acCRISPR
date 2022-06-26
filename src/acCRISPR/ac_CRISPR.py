from CS_FS_guide import counts_to_log2fc
from CS_corr_lib import generate_CS_corr_lib
from FS_gene_calc import FS_gene
from FS_to_p_val import p_value_calc
import argparse
import logging

logging.basicConfig(level = logging.INFO, format = '%(name)-5s %(asctime)s: %(message)s', datefmt = '%m/%d/%Y %I:%M:%S %p')
LOG = logging.getLogger('acCRISPR')
LOG.setLevel(logging.INFO)

def acCRISPR_parser():
    """ Creating command-line parser """
    parser = argparse.ArgumentParser(formatter_class = argparse.HelpFormatter, description = 'DESCRIPTION: Essential gene identification using acCRISPR')
    parser.add_argument('--counts',
                        type = str,
                        default = '',
                        metavar = '',
                        help = 'A TAB delimited file containing raw sgRNA counts for each sample')
    parser.add_argument('--replicate_info',
                        type = str,
                        default = '',
                        metavar = '',
                        help = 'A TAB delimited file containing replicate information for each sample')
    parser.add_argument('--cov',
                        type = int,
                        metavar = '',
                        help = 'Average no. of sgRNA per gene in the original library rounded off to the nearest integer')
    parser.add_argument('--cutoff',
                        type = float,
                        default = None,
                        metavar = '',
                        help = 'CS cutoff to be used to generate a CS-corrected library')
    parser.add_argument('--skip_log2fc_calc',
                        action = 'store_true',
                        default = False,
                        help = 'Skip sgRNA CS & FS calculation from raw counts')
    parser.add_argument('--CS_FS_file',
                        type = str,
                        default = '',
                        metavar = '',
                        help = 'File containing sgRNA ID, associated gene name, CS & FS (required only if log2fc is to be skipped)')
    parser.add_argument('--n_noness',
                        type = int,
                        default = 1000,
                        metavar = '',
                        help = 'No. of genes to be randomly sampled from the set of putative non-essential genes (default: 1000)')
    parser.add_argument('--p_val',
                        type = float,
                        default = 0.05,
                        metavar = '',
                        help = 'FDR corrected p-value cutoff for significance (default: 0.05)')
    parser.add_argument('--significance',
                        type = str,
                        default = '1-tailed',
                        choices = ['1-tailed', '2-tailed'],
                        help = 'Method of significance testing: 1-tailed or 2-tailed (default: 1-tailed)')
    parser.add_argument('--output_prefix',
                        type = str,
                        default = '',
                        metavar = '',
                        help = 'Prefix for output files')
    return parser

def run_ac_CRISPR():
    """ Running the tool and providing output """
    parser = acCRISPR_parser()
    args = parser.parse_args()
    LOG.info('       #####       Running acCRISPR       #####')
    if args.skip_log2fc_calc == False:
        LOG.info('Calculating sgRNA CS & FS')
        counts_to_log2fc(args.counts, args.replicate_info, args.output_prefix)
    
    if args.cutoff != None:
        if args.skip_log2fc_calc == True: CS_FS_file = args.CS_FS_file # CS_FS_file to be provided by the user
        else: CS_FS_file = args.output_prefix + '_guide_CS_FS.tab'
        LOG.info('Generating CS-corrected library')
        generate_CS_corr_lib(CS_FS_file, args.cutoff, args.output_prefix)
        LOG.info('Calculating FS of genes')
        FS_gene(args.output_prefix + '_corrected_guide_CS_FS.tab', args.output_prefix)
        LOG.info('Computing p-values for genes')
        CS_corr_coverage, p_corr = p_value_calc(args.output_prefix + '_corrected_guide_CS_FS.tab', args.output_prefix + '_FS_genes.tab', args.cov, args.n_noness, args.significance, args.output_prefix)
        
    if args.cutoff == None: 
        if args.skip_log2fc_calc == True: CS_FS_file = args.CS_FS_file # CS_FS_file to be provided by the user
        else: CS_FS_file = args.output_prefix + '_guide_CS_FS.tab'
        LOG.info('Calculating FS of genes')
        FS_gene(CS_FS_file, args.output_prefix)
        LOG.info('Computing p-values for genes')
        CS_corr_coverage, p_corr = p_value_calc(CS_FS_file, args.output_prefix + '_FS_genes.tab', args.cov, args.n_noness, args.significance, args.output_prefix)
    
    signif_genes = 0 # Variable for storing number of significant genes
    for i in range(len(p_corr)):
        if p_corr[i] < args.p_val: signif_genes += 1
    result = [];
    result.append('Significant genes = ' + str(signif_genes) + '\n')
    if args.cutoff != None:
        result.append('Average coverage of the CS-corrected library = ' + str(CS_corr_coverage) + '\n')
        result.append('ac-coefficient = ' + str(CS_corr_coverage*args.cutoff))
        
    else:
        result.append('Average coverage of the library = ' + str(CS_corr_coverage) + '\n')
        result.append('ac-coefficient = N/A')
    
    LOG.info('Writing final results')
    g = open(args.output_prefix + '_signif_genes_results.txt','w')
    g.writelines(result)
    g.close()
    LOG.info('Done!')
