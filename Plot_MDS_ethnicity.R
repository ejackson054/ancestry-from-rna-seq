library(graphics)
library(hash)
library(ggplot2)

#Read in plink distance table
table = read.table("/Users/emilyjackson/Documents/Ancestry_project/plink.mdist")
head(table)

#Run multi-dimensional scaling
fit <- cmdscale(table, k=2)
x <- -fit[,1]
y <- fit[,2]

#Read in the IDs for the plink distance table
ids = read.table("/Users/emilyjackson/Documents/Ancestry_project/plink.mdist.id",stringsAsFactors=FALSE)

#Read in 1000 genomes metadata
#Note:  There are some lines with the comment mark ("#") in them
genomes <- read.table("/Users/emilyjackson/Documents/1000genomes.sequence.index", sep = "\t", comment = "#", stringsAsFactors = FALSE, fill = TRUE)
genomes <- genomes[!duplicated(genomes$V10), ]
rownames(genomes) <- genomes$V10

#Make dictionary of populations and superpops (https://www.internationalgenome.org/category/population/)
pop_to_super = hash("CHB"="EAS","JPT"="EAS","CHS"="EAS","CDX"="EAS","KHV"="EAS",
                    "CEU"="EUR","TSI"="EUR","FIN"="EUR","GBR"="EUR","IBS"="EUR",
                    "YRI"="AFR","LWK"="AFR","GWD"="AFR","MSL"="AFR","ESN"="AFR","ASW"="AFR","ACB"="AFR",
                    "MXL"="AMR","PUR"="AMR","CLM"="AMR","PEL"="AMR",
                    "GIH"="SAS","PJL"="SAS","BEB"="SAS","STU"="SAS","ITU"="SAS")

#Get super-population for each individual in the plink distance table
super_pops = c()
inds = c()
pops = c()
for (i in 1:length(ids$V1)) {
  ind = ids$V1[i]
  inds = c(inds,ind)
  if (substring(ind,1,2) %in% c("NA","HG")) {
    pop = genomes[[ind,"V11"]]
    super_pop = pop_to_super[[pop]]
    super_pops = c(super_pops, super_pop)
    pops = c(pops, pop)
  }
  else {
    super_pops = c(super_pops,"SCA_samples")
    pops = c(pops, "SCA_samples")
  }
}

#Assemble dataframe (x, y, plus population)
df = data.frame(cbind(x, y, super_pops,inds,pops), stringsAsFactors=FALSE)
df = transform(df, x = as.numeric(x))
df = transform(df, y = as.numeric(y))

non_white = c("7179CX","7202CXXY","7247DXX","7348CXXY",
              "7165GXY","7207DXiXq","7169CXY","7184DX")
non_white_df = subset(df,inds %in% non_white)

ggplot(df, aes(x=x,y=y,color = super_pops)) + geom_point() +
  geom_text(data=non_white_df,aes(x=x,y=y,label=inds),hjust=0,vjust=0) +
  scale_color_manual(values=c("#F8766D","#B79F00","#00BA38","#00BFC4",
                              "#619CFF","#000000"))

#Plot by JUST POPULATION (not super population)
ggplot(df, aes(x=x,y=y,color = pops)) + geom_point() +
  geom_text(data=non_white_df,aes(x=x,y=y,label=inds),hjust=0,vjust=0)
    


#Repeat all of this using JUST 1000 GENOMES
table = head(table,-50)
n = dim(table)[2]
table = table[,1:(n-50)]

ids = head(ids,-50)


#Run multi-dimensional scaling
fit <- cmdscale(table, k=2)
x <- fit[,1]
y <- fit[,2]


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
    super_pops = c(super_pops,"SCA_samples")
  }
}

#Assemble dataframe (x, y, plus population)
df = data.frame(cbind(x, y, super_pops), stringsAsFactors=FALSE)
df = transform(df, x = as.numeric(x))
df = transform(df, y = as.numeric(y))

ggplot(df, aes(x=x,y=y,color = super_pops)) + geom_point() + 
  scale_color_manual(values=c("#F8766D","#B79F00","#00BA38","#00BFC4",
                              "#619CFF","#000000"))




