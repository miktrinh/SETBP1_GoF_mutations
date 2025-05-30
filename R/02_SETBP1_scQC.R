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

#     qc.summary = QC.output[[2]]
  
  
#   ## add metadata to srat object
#   srat = QC.output[[1]]
#   srat$channelID = srat$orig.ident
#   srat$cellID = rownames(srat@meta.data)
  
#   srat$donorID = 'NA'
#   srat$donorID[srat$channelID %in% c('c09','c10','c11','c75','c76','c77','c78')] = 'GOSH84'
#   srat$donorID[srat$channelID %in% c('c34','c35','c36','c37')] = 'BJ9'
#   srat$batch = 'NA'
#   srat$batch[srat$channelID %in% c('c34','c35','c36','c37')] = 'BJ9_1'
#   srat$batch[srat$channelID %in% c('c75','c76','c77','c78')] = 'GOSH84_1'
#   srat$batch[srat$channelID %in% c('c09','c10','c11')] = 'GOSH84_2'
#   srat$sex = 'NA'
#   srat$sex[srat$donorID == 'GOSH84'] = 'XX'
#   srat$sex[srat$donorID == 'BJ9'] = 'XY'
#   srat$age = 'NA'
#   srat$age[srat$donorID == 'GOSH84'] = '10'
#   srat$age[srat$donorID == 'BJ9'] = '2'
  
#   saveRDS(srat,'/lustre/scratch117/casm/team274/mt22/SETBP1/Results/1_setbp1QC/jan23/SETBP1_QC_clean_noMTCells.RDS')
  
#   ## Add mdat to highMT cells srat too
#   srat_highMT = readRDS('/lustre/scratch117/casm/team274/mt22/SETBP1/Results/1_setbp1QC/jan23/SETBP1_QC_clean_withMTCells.RDS')
#   srat_highMT$cellID = rownames(srat_highMT@meta.data)
#   m = match(srat_highMT$orig.ident,srat$orig.ident)
#   sum(is.na(m))
  
#   srat_highMT$channelID = as.character(srat$channelID[m])
#   srat_highMT$donorID = as.character(srat$donorID[m])
#   srat_highMT$batch = as.character(srat$batch[m])
#   srat_highMT$sex = as.character(srat$sex[m])
#   srat_highMT$age = as.character(srat$age[m])
  
#   saveRDS(srat_highMT,'/lustre/scratch117/casm/team274/mt22/SETBP1/Results/1_setbp1QC/jan23/SETBP1_QC_clean_withMTCells.RDS')
  
# }

