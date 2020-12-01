"""Generate filtered VCF file for each RNA-Seq sample


Usage:  python pipeline_with_filtering.py [individual] [fastq_1] [fastq_2] [main_dir]

individual:  Sample for which you are generating filtered VCF
[fastq_1]:  Full path to fastq_1
[fastq_2]: Full path to fastq_2
[main_dir]: Base directory for all commands


ENSURE THESE ARE TRUE BEFORE USE:
- hg38_chr.fa must be in main_dir
- Homo_sapiens.GRCh38.101.chr.gtf must be in main_dir
- Opposum.py must be in directory named Opposum which is accessed from main_dir
- Platypus.py must be in Platypus/bin which is accessed from main_dir


This pipeline is based on Barral-Arca et al. 2019, "Ancestry patterns inferred from massive RNA-seq data" (RNA)

"""

import os
import sys


#Read in arguments
individual = sys.argv[1]
fastq_1 = sys.argv[2]
fastq_2 = sys.argv[3]
main_dir = sys.argv[4]


#Make STAR index (only needs to be done once)
if not os.isdir("STAR_99"):
    os.system("mkdir STAR_99")
    os.system("STAR --runThreadN 8 --runMode genomeGenerate --genomeDir STAR_99 --genomeFastaFiles " + main_dir + "hg38_chr.fa --sjdbGTFfile " + main_dir + "Homo_sapiens.GRCh38.101.chr.gtf --sjdbOverhang 99")


#Run STAR alignment
os.system("mkdir " + individual)   
os.system("STAR --runThreadN 8 --genomeDir STAR_99 --readFilesIn " + fastq_1 + " "  + fastq_2 + " --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix " + individual + "/")

#Run OPOSSUM
os.system("samtools calmd -b " + individual + "/Aligned.sortedByCoord.out.bam " + main_dir + "hg38_chr.fa > " + individual+ "/Aligned.sortedByCoord.out.md.bam")

os.system("python2 " + main_dir + "Opossum/Opossum.py --BamFile=" + individual + "/Aligned.sortedByCoord.out.md.bam --SoftClipsExist=True  --OutFile=" + individual + "/opossum_output.bam")

#Run PLATYPUS
os.system("python2 " + main_dir + "Platypus/bin/Platypus.py callVariants --bamFiles " + individual + "/opossum_output.bam --refFile " + main_dir + "hg38_chr.fa --filterDuplicates 0 --minMapQual 0 --minFlank 0 --minGoodQualBases 10 --minBaseQual 20 -o " + individual + "/platypus_variants.vcf --source " + main_dir + "1000g_vcf/1000g.vcf.gz --minPosterior 0 --getVariantsFromBAMs 0")


#Filter VCF file to remove low-quality variants

with open("{0}/platypus_variants.vcf".format(individual),"r") as old:
    with open("{0}/platypus_variants_filtered.vcf".format(individual),"w") as new:
        for line in old:
            if line[0]=="#":
                new.write(line)
            else:
                fields = line.split("\t")
                if fields[-1].split(":")[0]=="0/0":
                    if "MQ" not in fields[6].split(";"):
                        new.write(line)
                elif fields[-1].split(":")[0]!="./.":
                    if ("Q20" not in fields[6].split(";")):
                        new.write(line)

#Replace sample name
new_name = individual.split("_")[0] + individual.split("_")[1]
os.system("sed -i 's/opossum_output/{0}/g' {1}/platypus_variants_filtered.vcf".format(new_name, individual))  
                   
#Zip VCF files
os.system("bgzip {0}/platypus_variants.vcf".format(individual))
os.system("bgzip {0}/platypus_variants_filtered.vcf".format(individual))

#Index filtered VCF for downstream steps analysis 
os.system("tabix -p vcf {0}/platypus_variants_filtered.vcf.gz".format(individual))