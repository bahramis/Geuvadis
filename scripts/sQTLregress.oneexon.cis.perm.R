#!/cm/shared/apps_chop/R/3.5.1/bin/Rscript

#############
# sQTL analysis for SNPs within 200kb window of one single exon


## Arguments needed for program listed at the command line
# exoninforfile
# alljunctionfile
# IJfile
# genoplinkfile
#chromosome   "chr#"
# exonindex
# exonindex_end
# permutation id
args = commandArgs(TRUE)

exoninforfile = args[1]
alljunctionfile = args[2]
IJfile= args[3]
genoplinkfile= args[4]
chrom=args[5]
exonindex = as.numeric(args[6])
exonindex_end = as.numeric(args[7])
perm_id=as.numeric(args[8])
phenofile.IJ = IJfile
phenofile = alljunctionfile
genofile = genoplinkfile

#print (c(exonindex, exonindex_end))

################################
#  read in the exon information

exon.infor <- read.csv(exoninforfile,header=T,sep="\t")
exonIDs <- as.character(exon.infor[,1])
chrstr <- exon.infor[, 4]
chr<- sub("chr","",chrstr)
SNPstartpos=exon.infor[,6]-200000
SNPendpos=exon.infor[,7]+200000
exonnames <- exonIDs 

# select the right se
selectchr=exon.infor[chrstr==chrom,]
if (exonindex_end > nrow(selectchr)) exonindex_end=nrow(selectchr)
selected.se=as.vector(as.matrix(selectchr[exonindex:exonindex_end,1]))
selected.index=c(1:nrow(exon.infor))[exonIDs %in% selected.se]

TMPASSODIR <-paste("../Glimmps_each_exon_cis_perm",perm_id, sep="")
library(lme4) 

source("../../GLiMMPS_functions.R")

###############################################################################################################################################
# read in phenotype expression levels in plink format
allreads.data <- read.csv( phenofile ,header=T,sep="\t")
nsample <- dim(allreads.data)[1] 
nexons <-dim(allreads.data)[2] -2

allreads.matrix <- allreads.data[,seq(1,nexons)+2 ]
IDs.pheno <- as.character(allreads.data[,1])      									# I made change here

IJ.data <- read.csv( phenofile.IJ ,header=T,sep="\t")
IJ.matrix <- IJ.data[,seq(1,nexons)+2 ]
### create output directory ##
if (!file.exists( TMPASSODIR )) {system(paste ("mkdir",TMPASSODIR)) }

  
# do the association for all SNPs near the target exon

cmmd <- paste("mkdir ",TMPASSODIR,"/",chrom,sep="")
if (!file.exists( paste(TMPASSODIR,"/",chrom,sep="")) ) {system(cmmd)}

##############
# read in genotype information

map<-read.csv(paste("../genotype/",chrom,".",genofile, ".map", sep=""), header=F, sep="\t")       #  I changed here 
names(map) <- c("Chr","SNPID","cM","Pos")

genoplink <- read.csv(paste("../genotype/perm", perm_id, "/", chrom, ".", genofile, "_tposed.raw", sep=""), header=T,sep=" ",  na.strings=c("NA"))     			 
genoplink=genoplink[,-1] 
 
IDs.geno <- names(genoplink)           #  I changed here 
IDs.common <- IDs.geno[is.element(IDs.geno, IDs.pheno)]
nsnps <- dim(genoplink)[1] -5                              #  I changed here 
sub.geno <- match(IDs.common , IDs.geno)
geno <- as.matrix(genoplink[seq(1, nsnps)+5, sub.geno])    #  I changed here 

### only take those individuals with genotypes, and sort the phenotype to be the same order as in the genotype matrix


sub <- match(IDs.common , IDs.pheno)
#IDs.pheno[sub]


###############
#  start association for psi~ SNP for every SNP. 


for (i in selected.index){


    n<- allreads.matrix[sub,i]
    y <- IJ.matrix[sub,i]
    
 	snp.filter=map[,1]==chr[i] & map[,4]>SNPstartpos[i] & map[,4]<SNPendpos[i]
 	if (length(snp.filter[snp.filter])==0){print (paste("nocis", exonIDs[i], sep=" ")); next;}
 	print (paste("start", exonIDs[i], sep=" "))
 	pvals.glmm <- rep(NA,nsnps)
    betas <-  rep(NA,nsnps) 
for ( gi in  c(1: nsnps)[snp.filter]) {         #701:720){ #
 #cat(paste("SNP ", gi,"\n"))
  ### data association
  onedata <- list(n=n,y=y,SNP=  as.vector(geno[gi,]) )
  #results.glm <- glm.sQTL ( onedata )
  
  #pvals.glm[gi] <- results.glm$pval

  #results.quasi <- glmquasi.sQTL ( onedata )
  #pvals.glmquasi[gi] <- results.quasi$pval

  ###############
  # GLiMMPS method 
  results.glmm <- glmm.sQTL ( onedata )
  pvals.glmm[gi] <- results.glmm$pval
  betas[gi] <- results.glmm$betas[2]
  ############
  #results.lm <- lm.sQTL ( onedata )
  #pvals.lm[gi] <- results.lm$pval
  
  #results.glmmWald <- glmmWald.sQTL ( onedata )
  #pvals.glmmWald[gi] <- results.glmmWald$pval

 }

tmpout <- cbind(map[snp.filter,c(1,2,4)], formatC(pvals.glmm[snp.filter],format="e",digits=3), round(betas[snp.filter],3))

colnames(tmpout) <- c("Chr","SNPID","Pos","pvals.glmm","Beta")
write.table(tmpout  , paste(TMPASSODIR,"/",chrom,"/",exonIDs[i],".asso",sep=""),row.names=F,col.names=T, quote=F,sep="\t" )

}
######################
cat("Finished association!\n")

