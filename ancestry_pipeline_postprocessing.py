"""Perform PLINK and ADMIXTURE analysis for filtered VCF files


Usage:  python ancestry_pipeline_postprocessing.py  [meta] [locations] [main_dir]

[meta]: Tab separated file with at least two columns: 1) sample name, 2) cell type
[locations]: Tab-separated file with two columns: 1) sample name, 2) .fastq file location
[main_dir]: Base directory for all commands


This pipeline is based on Barral-Arca et al. 2019, "Ancestry patterns inferred from massive RNA-seq data" (RNA)

"""


import os
import pandas as pd
import sys


"""Part I:  Get files to merge"""

meta = sys.argv[1]
locations = sys.argv[2]
main_dir = sys.argv[2]

meta = pd.read_csv(meta, sep  = "\t",index_col = "sample") 
locations = pd.read_csv(locations, sep = "\t")  

#Run fastqc on files that correspond to fibroblasts

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


"""Part II:  Merge files"""

variant_files = ["{0}/platypus_variants_filtered.vcf.gz".format(individual) for individual in list(dictionary.keys())]
number_files = len(variant_files)
variant_files = " ".join(variant_files)


#Find SNPs shared in both files       #SLOW:  Takes around 5 hr
os.system("bcftools isec -p {0} -n={1} 1000g_vcf/1000g.vcf.gz {2}".format(main_dir, number_files+1,variant_files))  #takes around 5 hr -f HapScore,PASS


#Zip and index files with shared SNPs
def get_filename(i):
    if i<10:
        string = "0{0}".format(str(i))
    else:
        string = str(i)
    filename = "00{0}.vcf".format(string)
    return(filename)

for i in range(0,number_files+1):
    filename = get_filename(i)
    os.system("bgzip {0}".format(filename))
    os.system("tabix -p vcf {0}.gz".format(filename))
    

#Merge all samples together into single file with shared SNPs          #MEDIUM (maybe 30 min)
vcf_files = ["{0}.gz".format(get_filename(i)) for i in range(0,number_files+1)]
vcf_files = " ".join(vcf_files)
os.system("bcftools merge {0} -Oz -o merged.vcf.gz".format(vcf_files)) #-f HapScore,PASS



"""Part III:  Run PLINK"""

#Filter to remove SNPs in LD (r2 > 0.75)  #FAST
os.system("plink --vcf merged.vcf.gz --recode vcf bgz --set-missing-var-ids @:# --out merged")  #-const-fid  #output:  merged.vcf.gzary
os.system("plink --vcf merged.vcf.gzary --indep-pairwise 50 5 0.75  --out merged_LD")
os.system("plink --vcf merged.vcf.gzary --extract merged_LD.prune.in --recode vcf bgz --out merged_LD")

os.system("plink --distance square 1-ibs --vcf merged_LD.vcf.gzary") #FAST



"""Part IV: Run Admixture"""

os.system("plink --vcf merged_LD.vcf.gzary --make-bed --out merged_LD") #2min
os.system("admixture merged_LD.bed 4") #FAST

