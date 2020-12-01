"""Wrapper for ancestry_pipeline.py


Usage:  python ancestry_pipeline_wrapper.py [meta] [locations] [main_dir]

[meta]: Tab separated file with at least two columns: 1) sample name, 2) cell type
[locations]: Tab-separated file with two columns: 1) sample name, 2) .fastq file location
[main_dir]: Base directory for all commands


This pipeline is based on Barral-Arca et al. 2019, "Ancestry patterns inferred from massive RNA-seq data" (RNA)

"""

import os
import pandas as pd
import sys



#Read in arguments
meta = sys.argv[1]
locations = sys.argv[2]
main_dir = sys.argv[2]

meta = pd.read_csv(meta, sep  = "\t",index_col = "sample") 
locations = pd.read_csv(locations, sep = "\t")  

#Run pipeline on files that correspond to fibroblasts

os.chdir(main_dir)

dictionary = {}

for i in locations.index:
    
    sample = locations.loc[i,"Sample"]
    filename = locations.loc[i,"File"]
    
    if meta.loc[sample,"cell_type"].strip()=="Fib":
        if sample not in dictionary.keys():
            dictionary[sample] = [filename]
        else:
            dictionary[sample].append(filename)

for individual in list(dictionary.keys())[0:1]:
    
    fastq_1 = dictionary[individual][0]
    fastq_2 = dictionary[individual][1]
    
    os.system("bsub python ancestry_pipeline.py {0} {1} {2}".format(individual,fastq_1, fastq_2))
        