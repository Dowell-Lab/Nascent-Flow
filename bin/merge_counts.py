from __future__ import print_function
import argparse
import numpy as np
import datetime
import pandas as pd
import chartify
import glob
import sys
import os
from functools import reduce
from pybedtools import BedTool
from argparse import RawTextHelpFormatter
from multiprocessing import Pool

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Merge Count\n\n===============================================================================================================\n\nA short and simple python script to merge count files generated from multicov.', epilog='@Dowell Lab, Margaret Gruca', usage='%(prog)s --beds /my/sample/dir/ --output /my/out/file.bed', formatter_class=RawTextHelpFormatter)
    
    required = parser.add_argument_group('Required Arguments')
    optional = parser.add_argument_group('Optional Arguments')
    
    required.add_argument('-b', '--bed', dest='beddir', \
                       help='Path to count file(s) of interest.', required=True)
    
    
    required.add_argument('-o', '--output', dest='output', metavar='SAMPLE.BED', \
                        help='Path to where output BED file will be saved.', required=True)
    
    optional.add_argument('-t', '--threads', type=int, dest='threads', metavar='<THREADS>', help='Number of threads for multi-processing. Default=1', default=1, required=False)
    
    args = parser.parse_args()

def reader(f):
    sample = str(os.path.splitext(os.path.basename(f))[0])
    return pd.read_csv(f, sep="\t", usecols=[0,1,2,3,12], \
                         names=['chr', 'start', 'stop', 'id', sample])


def main():
    pool = Pool(args.threads) # number of cores you want to use
    file_list = sorted(glob.glob(args.beddir + '*.bed'))
    df_list = pool.map(reader, file_list) #creates a list of the loaded df's
    df_merged = reduce(lambda  left,right: pd.merge(left,right,on=['chr', 'start', 'stop','id'], how='outer'), df_list)
    df_merged.to_csv((args.output), sep="\t", index=False)

if __name__ == '__main__':
    main()