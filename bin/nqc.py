from __future__ import print_function
import argparse
import numpy as np
import datetime
import pandas as pd
import chartify
import glob
import sys
import os
import concurrent.futures
import warnings
import subprocess
from operator import itemgetter
from argparse import RawTextHelpFormatter
from multiprocessing import Pool
from io import StringIO

def get_regression_data(params):
    
    (bedgraph,genome) = params
    sample_name=os.path.splitext(os.path.basename(bedgraph))[0]  
    
    coverage_file = pd.read_csv(bedgraph, names=['chr','start','stop','coverage'], sep='\t' \
                                 ).assign(sample=sample_name)
    coverage_file['log2_coverage'] = np.log2(abs(coverage_file.coverage))
    coverage_variance = coverage_file.log2_coverage.var(ddof=1)    
    
    complement = subprocess.Popen(['complementBed', '-i', bedgraph, '-g', genome], 
               stdout=subprocess.PIPE, 
               stderr=subprocess.STDOUT)
    compout,comperr = complement.communicate()
    compdata = StringIO(str(compout,'utf-8')) 
    gaps=pd.read_csv(compdata, sep='\t', header=None, \
                  names = ['chr','start','end'])    
    gaps['gap_size'] = gaps.end - gaps.start
    gap_std = np.log2(gaps.gap_size.std(ddof=1))
    gap_var = np.log2(gaps.gap_size.var(ddof=1))
    gap_median_size = np.log2(gaps['gap_size'].median())
    number_of_1bp_gaps = np.log2((gaps.gap_size == 1).sum()) 
    
    coverage_stats = { 'sample': sample_name, \
                         'gap_std': float(gap_std), \
                         'gap_var': float(gap_var), \
                         'gap_median_size': float(gap_median_size), \
                         'coverage_variance': float(coverage_variance), \
                         'number_of_1bp_gaps': float(number_of_1bp_gaps) }  
    
    print("Done with sample %s." % sample_name)
    
    return coverage_stats


def main():
    parser = argparse.ArgumentParser(description='Nascent QC\n\n===============================================================================================================\n\nAn algorithm which provides a metric for nascent data quality.', epilog='@Dowell Lab, Margaret Gruca, margaret.gruca@colorado.edu\nFor questions and issues, see https://github.com/Dowell-Lab/nqc', usage='%(prog)s --bedgraph /my/sample/dir/*.bedGraph --output /my/out/file.bed', formatter_class=RawTextHelpFormatter)
    
    required = parser.add_argument_group('Required Arguments')
    optional = parser.add_argument_group('Optional Arguments')
    
    required.add_argument('-b', '--bedgraph', dest='bedgraphdir', \
                        help='Path to bedgraph file(s) of interest. Wildcard may be used to process multiple bedGraphs simultaneously with one output report.', required=True)
    
    required.add_argument('-o', '--output', dest='output', metavar='SAMPLE.BED', \
                        help='Path to where output BED file will be saved.', required=False, default='nqc_stats.txt')
    
    required.add_argument('-g', '--genome', dest='genome', \
                        help='Chromosome sizes file -- bedGraph file sort must match genome file. See README for details on generating this file for your specific genome.', \
                        default='', required=False, type=str)    
   
    optional.add_argument('-t', '--threads', type=int, dest='threads', metavar='<THREADS>', \
                        help='Number of threads for multi-processing. Default=1', default=1, required=False)
    
    args = parser.parse_args()
    
    genome = args.genome
    warnings.simplefilter(action='ignore', category=(FutureWarning,RuntimeWarning))
    
    threads = args.threads
    file_list = glob.glob(args.bedgraphdir + '*.bedGraph')   
    
    print("Starting to process bedGraphs...\n" + str(datetime.datetime.now()))
    
    with concurrent.futures.ProcessPoolExecutor(threads) as executor:
        jobs = [executor.submit(get_regression_data, [file,genome])
                for file in file_list]
        regression_stats = [r.result() for r in jobs]
        
    sorted_regression_stats = sorted(regression_stats, key=itemgetter('gap_std'))

    regression_stats_file = open(args.output, 'w')
    for stat in sorted_regression_stats:
        regression_stats_file.write("%s\t%s\t%s\t%s\t%s\t%s\n" % \
                        (stat['sample'], stat['gap_std'], stat['gap_var'], \
                        stat['gap_median_size'], stat['coverage_variance'], stat['number_of_1bp_gaps']))
    regression_stats_file.close()
    
    print("Done.\n" + str(datetime.datetime.now()))
    
    sys.exit(0)
   
if __name__ == '__main__':
    main()