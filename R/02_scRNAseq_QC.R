# QC SETBP1 single-cell RNASeq dataset

outDir = 'Results/02_SETBP1_scQC'
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}

##----------------##
##   Libraries  ####
##----------------##
library(tidyverse)
library(Seurat)
#source("R/utils/misc.R")
source("R/utils/sc_basicQC.R")

##----------------------------##
##   Set Global parameters  ####
##----------------------------##

seed = 2397


maxMT = 30
minGenes = 200
minUMIs = 300
maxBadFrac = 0.5
numPCs = 75
clusteringRes = 10
skipScrub = F
skipSoup=F
scrubScoreMax = 0.5
scrubPath='./scrubletScores.tsv'
scPath="./strainedCounts"
doPlot=T
verbose = T
skipIfExists=F
keepMTCells=T
  
outDir_sub = file.path(outDir,'2505')  
plotDir = file.path(outDir_sub,'SETBP1_SGS_')
outPath = file.path(outDir_sub,paste0('SETBP1_SGS'))

cleanSrat_fp = ifelse(keepMTCells,paste0(outPath,'_clean_withMTCells.RDS'),paste0(outPath,'_clean_noMTCells.RDS'))
if(file.exists(cleanSrat_fp) & skipIfExists){
    cleanSrat = readRDS(cleanSrat_fp)  
}else{
    if(!dir.exists(outDir_sub)){
        message(sprintf('Creating output directory'))
        dir.create(outDir_sub,recursive = T)
    }
    
    
    # please get the rest of the script template from aneuploidy/scripts/sampleQC_2.R
    dataDirs = file.path(list.files('Data/inhouse_SETBP1_10X',full.names = T,pattern = 'GRCh38-2020-A'),'/filtered_feature_bc_matrix')
    # Remove MDS-SETBP1 sample
    dataDirs = dataDirs[!grepl('SB_Leuk13645525',dataDirs)]
    # Remove BJ9-solidTumour snRNA-seq samples
    dataDirs = dataDirs[!grepl('MY_200531_13839752|MY_200531_13839753|MY_200531_13839754',dataDirs)]    

    names(dataDirs) = gsub('^.*_CG_|_GRCh38-2020-A.*$','',dataDirs)    
    names(dataDirs) = gsub('^.*_MY_','MY_',names(dataDirs))    
    names(dataDirs) = gsub('^.*_SB_','SB_',names(dataDirs))
    names(dataDirs) = gsub('_','.',names(dataDirs))    
    
  
    dataDirs=dataDirs[file.exists(dataDirs)]
    dataDirs=dataDirs[sapply(dataDirs, function(x){length(list.files(x))>0})]
    print(n_distinct(dataDirs))
  
    cleanCountDir = dataDirs
    metadata = NULL
    matchBy = NULL
    is10X=TRUE
    
    # Run basicQC
    message('\nPerforming scRNAseq QC...')
    
    QC.output = basicQC(dataDirs = dataDirs,maxMT = maxMT,
                    minGenes=minGenes,minUMIs=minUMIs,maxBadFrac=maxBadFrac,numPCs=numPCs,clusteringRes=clusteringRes,
                    cleanCountDir=cleanCountDir,
                    skipScrub=skipScrub,skipSoup=skipSoup,scrubScoreMax=scrubScoreMax,scrubPath=scrubPath,
                    metadata=metadata,matchBy=matchBy,scPath=scPath,outPath=outPath,skipIfExists=skipIfExists,
                    doPlot=doPlot,plotDir=plotDir,verbose=verbose,is10X=is10X)
  
    cleanSrat = QC.output[[1]]
    df.out = QC.output[[2]]
    
    write.csv(df.out,paste0(outPath,'qc_summary.csv'))
}


