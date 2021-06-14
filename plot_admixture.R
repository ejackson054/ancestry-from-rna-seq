
library(RColorBrewer)
library(hash)

#Visualization based on this paper: https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_4

#Read in file:  Results from Admixture program (.4.Q)
table = read.table("/Users/emilyjackson/Documents/Ancestry_project/merged_LD.4.Q")

#Read in the IDs for the plink distance table
ids = read.table("/Users/emilyjackson/Documents/Ancestry_project/plink.mdist.id",stringsAsFactors=FALSE)

#Read in 1000 genomes metadata
#Note:  There are some lines with the comment mark ("#") in them
genomes <- read.table("/Users/emilyjackson/Documents/1000genomes.sequence.index", sep = "\t", comment = "#", stringsAsFactors = FALSE, fill = TRUE)
genomes <- genomes[!duplicated(genomes$V10), ]
rownames(genomes) <- genomes$V10

#Re-organize by ancestry
#Make dictionary of populations and superpops (https://www.internationalgenome.org/category/population/)
pop_to_super = hash("CHB"="EAS","JPT"="EAS","CHS"="EAS","CDX"="EAS","KHV"="EAS",
                    "CEU"="EUR","TSI"="EUR","FIN"="EUR","GBR"="EUR","IBS"="EUR",
                    "YRI"="AFR","LWK"="AFR","GWD"="AFR","MSL"="AFR","ESN"="AFR","ASW"="AFR","ACB"="AFR",
                    "MXL"="AMR","PUR"="AMR","CLM"="AMR","PEL"="AMR",
                    "GIH"="SAS","PJL"="SAS","BEB"="SAS","STU"="SAS","ITU"="SAS")

#Get super-population for each individual in the plink distance table
super_pops = c()
for (i in 1:length(ids$V1)) {
  ind = ids$V1[i]
  if (substring(ind,1,2) %in% c("NA","HG")) {
    pop = genomes[[ind,"V11"]]
    super_pop = pop_to_super[[pop]]
    super_pops = c(super_pops, super_pop)
  }
  else {
    super_pops = c(super_pops,"test_data")
  }
}

table$population = super_pops
table = table[c("V1","V3","V4","V2","population")]

table_SAS = subset(table,population=="SAS")
table_AFR = subset(table,population=="AFR")
table_EUR = subset(table,population=="EUR")
table_EAS = subset(table,population=="EAS")
table_AMR = subset(table,population=="AMR")
table_test = subset(table,population=="test_data")

total <- rbind(table_AFR,table_EUR,table_SAS,table_EAS,table_AMR,table_test)

#Plot separately for each ancestry group

#Americas
par(mar=c(1.5,4,2.5,2),cex.lab=0.75,cex.axis=0.6)
barplot(t(as.matrix(subset(table_AMR,select=-c(population)))),
        col=brewer.pal(4,"Set1"),ylab="Anc.Proportions",border=NA,space=0)

#African
par(mar=c(1.5,4,2.5,2),cex.lab=0.75,cex.axis=0.6)
barplot(t(as.matrix(subset(table_AFR,select=-c(population)))),
        col=brewer.pal(4,"Set1"),ylab="Anc.Proportions",border=NA,space=0)

#East Asian
par(mar=c(1.5,4,2.5,2),cex.lab=0.75,cex.axis=0.6)
barplot(t(as.matrix(subset(table_EAS,select=-c(population)))),
        col=brewer.pal(4,"Set1"),ylab="Anc.Proportions",border=NA,space=0)

#European
par(mar=c(1.5,4,2.5,2),cex.lab=0.75,cex.axis=0.6)
barplot(t(as.matrix(subset(table_EUR,select=-c(population)))),
        col=brewer.pal(4,"Set1"),ylab="Anc.Proportions",border=NA,space=0)

#South Asian
par(mar=c(1.5,4,2.5,2),cex.lab=0.75,cex.axis=0.6)
barplot(t(as.matrix(subset(table_SAS,select=-c(population)))),
      col=brewer.pal(4,"Set1"),ylab="Anc.Proportions",border=NA,space=0)

#Individuals collected by our lab
par(mar=c(1.5,4,2.5,2),cex.lab=0.75,cex.axis=0.6)
barplot(t(as.matrix(subset(table_test,select=-c(population)))),
        col=brewer.pal(4,"Set1"),ylab="Anc.Proportions",border=NA,space=0)

#All individuals
par(mar=c(1.5,4,2.5,2),cex.lab=0.75,cex.axis=0.6)
barplot(t(as.matrix(subset(total,select=-c(population)))),
        col=brewer.pal(4,"Set1"),ylab="Anc.Proportions",border=NA,space=0)
