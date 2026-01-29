##### SETBP1 genotyping ######

# sudo apt-get install bcftools
# sudo cp /nfs/users/nfs_m/my4/bin/alleleCounter /usr/local/bin/alleleCounter
# sudo chmod +x /usr/local/bin/alleleCounter

setwd("/lustre/scratch117/casm/team274/mt22/SETBP1/")


#############
# Libraries #
#############
# Load libraries
library(tidyverse)
library(alleleIntegrator)

##############################################################################
# Genotyping to check if the samples are from the same or different patients #
##############################################################################

# We don't have DNA data so just look at RNA only

refGenome10X = '/nfs/users/nfs_m/my4/Projects/FetalCN/Data/scRNA/genomeRNA.fa'
liftChain = '/nfs/users/nfs_m/my4/Projects/FetalCN/Data/hg19ToHg38_noChr.over.chain'
nParallel=24


outDir = file.path('Results/0_SETBP1_genotypeCheck')
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}


#### Get list of RNA and DNA BAM files ####

#----- RNA BAMs
bams10X = paste0(c(list.files('/lustre/scratch117/casm/team274/mt22/SETBP1/Data/SETBP1/',full.names = T),
                    list.files('/lustre/scratch117/casm/team274/mt22/Data/SETBP1/',pattern = 'GRCh38-2020-A',full.names = T)),
                  '/possorted_genome_bam.bam')

names(bams10X) = gsub('_GRCh38-2020-A.*$','',gsub('^.*_SB_','',bams10X))
bams10X = bams10X[file.exists(bams10X)]

#############################
# Check genotype consistency
#Are all the BAMs you're going to use from the same individual?  Check before you start

genoCheck = matchBAMs(BAMs = c(bams10X),
                      refGenomes = rep(c(refGenome10X),c(length(bams10X))),
                      outputs = file.path(outDir,paste0(c(names(bams10X)),'_genotypeCheck.tsv')),
                      liftOvers=rep(c(NA),c(length(bams10X))),
                      is10X=rep(c(TRUE),c(length(bams10X))),
                      nParallel=nParallel)
#If anything is less than 0.8 and you should be concerned...
message(sprintf("The minimum similarity found was %g",min(genoCheck$ibs$ibs)))




