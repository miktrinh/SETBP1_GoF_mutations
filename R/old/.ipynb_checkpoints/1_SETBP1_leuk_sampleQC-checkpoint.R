## Preprocessing of L061 - leukaemia sample with SETBP1 mutation
setwd('~/SETBP1_GoF_mutations/')

outDir = 'Results/04_public_scRNAseq/2601/L061'
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}

##----------------##
##   Libraries  ####
##----------------##
library(tidyverse)
library(Seurat)
source("R/utils/sc_basicQC.R")

# library(readxl)
# library(ggVennDiagram)
# library(SoupX)
# library(RColorBrewer)


##----------------------------##
##   Set Global parameters  ####
##----------------------------##

seed = 2397

maxMT = 30
minGenes = 300
minUMIs = 500
maxBadFrac = 0.5
numPCs = 75
clusteringRes = 10
skipScrub = F
skipSoup = F
scrubScoreMax = 0.5
scrubPath='scrubletScores.tsv'
scPath="strainedCounts"
doPlot=T
verbose = T
skipIfExists=T
keepMTCells=T
rho_max_limit=0.02

#cleanCountDir = '~/lustre_mt22/SETBP1/Results/1_setbp1QC/jun23/leukPBMC_cleanCount'
# Define dataDir
dataDirs = 'Data/L061_pleukPBMC/cellranger612_count_47089_SB_Leuk13645525_GRCh38-2020-A/filtered_feature_bc_matrix'
names(dataDirs) = gsub('.*SB_|_GR.*$','',basename(dirname(dataDirs)))

projMan = readxl::read_excel('~/lustre_mt22/projectManifest.xlsx',sheet = 'SETBP1')
projMan = projMan[!is.na(projMan$assay) & grepl('scRNA',projMan$assay) & 
                    projMan$DonorID == 'L061',]
metadata = projMan[,c('DonorID','externalID','Tissue','Sex','assay','sangerSampleID')]
colnames(metadata) = c('donorID','externalID','tissue','sex','assay','sangerSampleID')
metadata$channelID = gsub('^SB_','',metadata$sangerSampleID)

# Run basicQC
plotDir = file.path(outDir,paste0('L061_'))
outPath = file.path(outDir,paste0('L061_'))

cleanCountDir = NULL # Specify NULL to save scrublet and soupX results in the same dataDir folder
matchBy = c('orig.ident','chanelID')
is10X=TRUE

QC.output = basicQC(dataDirs = dataDirs,maxMT = maxMT, 
                    minGenes=minGenes,minUMIs=minUMIs,maxBadFrac=maxBadFrac,numPCs=numPCs,
                    clusteringRes=clusteringRes,cleanCountDir=cleanCountDir,
                    skipScrub=skipScrub,skipSoup=skipSoup,scrubScoreMax=scrubScoreMax,scrubPath=scrubPath,
                    metadata=metadata,matchBy=matchBy,scPath=scPath,rho_max_limit=rho_max_limit,
                    outPath=outPath,skipIfExists=skipIfExists,
                    doPlot=doPlot,plotDir=plotDir,verbose=verbose,is10X=is10X)

cleanSrat = QC.output[[1]]

df.out = QC.output[[2]]

#cleanSrat = readRDS('~/lustre_mt22/SETBP1/Results/1_setbp1QC/jun23/leukPBMC_clean_noMTCells.RDS')

# Process L067 GATA1 -----------------------------------------------------------
outDir = 'Results/04_public_scRNAseq/2601/L067'
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}

dataDirs = 'Data/L067_pleukPBMC/cellranger700_count_47089_SB_Leuk13645530_GRCh38-2020-A/filtered_feature_bc_matrix'
names(dataDirs) = gsub('.*SB_|_GR.*$','',basename(dirname(dataDirs)))

projMan = readxl::read_excel('~/lustre_mt22/projectManifest.xlsx',sheet = 'GOSH_others_included')
projMan = projMan[!is.na(projMan$assay) & !grepl('WGS',projMan$assay) & 
                    projMan$donorID == 'L067',]
metadata = projMan[,c('donorID','Tissue','Sex','assay','sangerSampleID')]
colnames(metadata) = c('donorID','tissue','sex','assay','sangerSampleID')
metadata$channelID = gsub('^SB_','',metadata$sangerSampleID)

# Run basicQC
plotDir = file.path(outDir,paste0('L067_'))
outPath = file.path(outDir,paste0('L067_'))

cleanCountDir = NULL # Specify NULL to save scrublet and soupX results in the same dataDir folder
matchBy = c('orig.ident','chanelID')
is10X=TRUE

QC.output = basicQC(dataDirs = dataDirs,maxMT = maxMT, 
                    minGenes=minGenes,minUMIs=minUMIs,maxBadFrac=maxBadFrac,numPCs=numPCs,
                    clusteringRes=clusteringRes,cleanCountDir=cleanCountDir,
                    skipScrub=skipScrub,skipSoup=skipSoup,scrubScoreMax=scrubScoreMax,scrubPath=scrubPath,
                    metadata=metadata,matchBy=matchBy,scPath=scPath,rho_max_limit=rho_max_limit,
                    outPath=outPath,skipIfExists=skipIfExists,
                    doPlot=doPlot,plotDir=plotDir,verbose=verbose,is10X=is10X)

cleanSrat = QC.output[[1]]

df.out = QC.output[[2]]


### Merge with combSrat
## SETBP1
combSrat = readRDS('Results/03_SETBP1_annotation/2505/SETBP1_SGS_clean_noMTCells_annot_2505.RDS')
combSrat_L061 = merge_seurat_objects(combSrat,leuk_srat,keepAllGenes = F,genomeVersions = c('v38','v38'))
combSrat_L061 = standard_clustering(combSrat_L061)

#combSrat = readRDS(('~/lustre_mt22/SETBP1/Results/2_setbp1_HCA_annotation/jan23/SETBP1_REF24_mergedProcessed_annotated.RDS'))
#combSrat = readRDS(('~/lustre_mt22/SETBP1/Results/2_setbp1_HCA_annotation/jan23/SETBP1_REF24_mergedProcessed_annotated_harm_jun23.RDS')) # harmonized version
annot = read.csv('~/lustre_mt22/SETBP1/Results/2_setbp1_HCA_annotation/jan23/SETBP1_REF24_mergedProcessed_annotated_jun23_mdat.csv') # generate by 9_reAnalysis.R
combSrat$annot = annot$finalAnn_jun23[match(combSrat$cellID,annot$cellID)]
combSrat$annot[is.na(combSrat$annot)] = 'NA'
combSrat$donorID[combSrat$dataset == 'HCA'] = combSrat$sampleID[combSrat$dataset == 'HCA']
combSrat$finalAnn_jun23 = combSrat$annot

combSrat$dataset2 = ifelse(combSrat$donorID == 'GOSH84','GOSH84',ifelse(combSrat$donorID == 'BJ9','BJ9','HCA'))
#combSrat$dataset2 = ifelse(combSrat$donorID == 'GOSH84','GOSH84',combSrat$dataset2)
combSrat$ann = ifelse(combSrat$finalAnn_jun23 %in% c('DC.precursor','DC1','HSC_MEMP','Mast_cell','unknown','Plasma.cell'),'unknown',
                      paste0(combSrat$finalAnn_jun23,':',combSrat$dataset2))




leuk_srat = readRDS('~/lustre_mt22/SETBP1/Results/1_setbp1QC/jun23/leukPBMC_clean_noMTCells.RDS')
leuk_srat = standard_clustering(leuk_srat)

# add metadata
leuk_srat@meta.data$donorID = metadata$donorID[match(leuk_srat$orig.ident,metadata$channelID)]
leuk_srat@meta.data$externalID = metadata$externalID[match(leuk_srat$orig.ident,metadata$channelID)]
leuk_srat@meta.data$tissue = metadata$tissue[match(leuk_srat$orig.ident,metadata$channelID)]
leuk_srat@meta.data$sex = metadata$sex[match(leuk_srat$orig.ident,metadata$channelID)]
leuk_srat@meta.data$assay = metadata$assay[match(leuk_srat$orig.ident,metadata$channelID)]
leuk_srat@meta.data$sangerSampleID = metadata$sangerSampleID[match(leuk_srat$orig.ident,metadata$channelID)]

#merge(leuk_srat@meta.data,metadata,by.x='orig.ident',by.y='channelID',all=T)
#rownames(leuk_srat@meta.data) = leuk_srat@meta.data$cellID
#annot = read.csv('/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/BALL_leuk_anno.csv')
#leuk_srat$externalID_2 = leuk_srat$externalID
#leuk_srat$externalID_2[leuk_srat$externalID_2 == 'L007_BALL_RelapseD0_Blood'] = 'L007_BALL_Relapse_Blood'
#m = match(leuk_srat$externalID,annot$sample_ID)
leuk_srat_markers = quickMarkers(leuk_srat@assays$RNA@counts,leuk_srat$seurat_clusters)
DimPlot(big.srat.sub_harm,cells.highlight = leuk_srat$cellID[leuk_srat$seurat_clusters==15])
FeaturePlot(leuk_srat,'CD14')
big.srat.sub_harm
table(leuk_srat$orig.ident[leuk_srat$seurat_clusters==10])
DimPlot(leuk_srat,group.by = 'seurat_clusters',label = T,label.box = T) + NoLegend()
DimPlot(leuk_srat,group.by = 'donorID',label = T,label.box = T) + NoLegend()
DimPlot(leuk_srat,group.by = 'annot',label = T,label.box = T) + NoLegend()

leuk_srat$seurat_clusters = as.character(leuk_srat$seurat_clusters)
leuk_srat$annot = as.character(leuk_srat$seurat_clusters)
leuk_srat$annot[leuk_srat$seurat_clusters %in% c('3','16','18')] = 'BALL'
leuk_srat$annot[leuk_srat$seurat_clusters %in% c('7','4','36','34','25','35','30')] = 'AML'
leuk_srat$annot[leuk_srat$seurat_clusters %in% c(5,10,19,15)] = as.character(leuk_srat$seurat_clusters[leuk_srat$seurat_clusters %in% c(5,10,19,15)])
leuk_srat$annot[leuk_srat$seurat_clusters %in% c(23,33,31)] = 'B.cell'
leuk_srat$annot[leuk_srat$seurat_clusters %in% c(21,22)] = 'normalT'
leuk_srat$annot[leuk_srat$seurat_clusters %in% c(20,32)] = 'BALL'
leuk_srat$annot[leuk_srat$seurat_clusters %in% c(12,11)] = '?'
leuk_srat$annot[leuk_srat$seurat_clusters %in% c(0,1,2,9,8,26,14,13,17,6)] = 'TALL'
leuk_srat$annot[leuk_srat$seurat_clusters %in% c(29)] = 'doublets'
leuk_srat$annot[leuk_srat$seurat_clusters %in% c(40)] = '?'

Idents(big.srat.sub_harm) = big.srat.sub_harm$annot2
DotPlot(big.srat.sub_harm,group.by = 'seurat_clusters',idents = 'BALL',features = unique(c('CD34','CD38','SPINK2','MLLT3','PRSS57', # HSC_MPP
                                                                  'SERPINB1', 'GATA1',	'GATA2', 'TESPA1',	'CTNNBL1',#MEMP
                                                                  'FCER1A', 'ITGA2B', 'HBD','KLF1','GYPA','PLEK', # MEP
                                                                  'ZBTB16','LTB', 'CD52',# ILC precursor
                                                                  'IL7R','CD3D','GZMA','CD4','CD8A', #Early.lymphoid_T.lymphocyte
                                                                  'TRDV2','TRGV9', # gamma delta T-cell
                                                                  'SLC4A10','TRAV1-2', #MAIT t-cell
                                                                  'PRF1', # effector T cells
                                                                  'FOXP3',	'CDH1', # regulatory T cells
                                                                  'NKG7','KLRD1', #NK
                                                                  
                                                                  'IGLL1','CD99', # Pre-pro
                                                                  'DNTT','CD79B','VPREB1','EBF1','CD19','RAG1',# pro-b-cells
                                                                  # naive_b_cell = CD19+ CD27-
                                                                  # memory_b_cell = CD19+ and CD27 + 
                                                                  'MME','CD79A',# pre-b-cells
                                                                  'TCL1A','MME','RAG1','MS4A1',  # B-cells
                                                                  
                                                                  'PF4','ITGA2B', #Megakaryocyte
                                                                  'CSF2RB','HDC','SERPINB1','TPSAB1','KIT', # Mast.cell
                                                                  'NRGN','PPBP','PF4', # Platelets
                                                                  
                                                                  'GATA1','KLF1', # early Erythroid
                                                                  'ALAS2', # mid.erythroid
                                                                  'HBA1','BPGM', # late.erythroid
                                                                  
                                                                  'DEFA3','DEFA4', # pro-myelocytes
                                                                  'CAMP','LCN2', #myelocytes
                                                                  'CXCR2','CSF3R','FCGR3B','FUT4', # Neutrophil
                                                                  
                                                                  'CTSG',	'PRTN3', # CMP
                                                                  'AZU1','MPO','FLT3','PTPRC', # GMP
                                                                  
                                                                  'ITGB2','SELL','ITGAM','CD14','FCGR3A',# Monocytes
                                                                  'CD68','MSR1','FTL','SELENOP','APOE','HMOX1','CD5L',#Macrophages
                                                                  
                                                                  'ITGAX','CD1C','MME', #DCs
                                                                  'IRF8',	'CLEC10A', #DC.precursor
                                                                  'CD1C','CLEC9A',#DC1
                                                                  'IL3RA', # DC2
                                                                  'CLEC4C', #pDC
                                                                  
                                                                  'CD27' # Plasma cells
                                                                  
                                                                  
)))+RotatedAxis() + ggtitle('Celltype marker genes')

leuk_srat$dataset2 = ifelse(leuk_srat$donorID == 'L061','L061','leuk')




#### Import Bram's leuk annotation
annot = read.csv('/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/BALL_leuk_anno_soupX_scrublet_comments.csv')
annot = annot[grepl('Blood',annot$sample_ID),]
# umap clusters definition
s.list = list()
ann = data.frame()
for(sample in list.files('/lustre/scratch126/casm/team274sb/bl10/B-ALL/Intermediate/16_UMAP_clustering',pattern = 'Blood',full.names = T)){
  s = readRDS(sample)
  a = annot[annot$sample_ID == gsub('\\.rds$','',basename(sample)),]
  s$annot = a$leuk_anno[match(as.character(s$seurat_clusters),a$seurat_clusters)]
  DimPlot(s,group.by = 'annot',label = T,label.box = T,repel = T,cols=col25[-6]) + NoLegend()
  s.list[[gsub('\\.rds$','',basename(sample))]] = s
  tmp = data.frame(cellID = s$cell_ID,sampleID = gsub('\\.rds$','',basename(sample)),annot=s$annot)
  
  if(nrow(ann) == 0){
    ann = tmp
  }else{
    ann = rbind(ann,tmp)
  }
}


ann$donorID = sapply(strsplit(ann$cellID,'_'),'[',1)
ann$channelID = leuk_srat$orig.ident[match(ann$donorID,leuk_srat$donorID)]
ann$cellID2 =  paste0(ann$channelID,'_',sapply(strsplit(ann$cellID,'::'),'[',2))




## Import my object
leuk_srat = readRDS('~/lustre_mt22/SETBP1/Results/1_setbp1QC/jul23/leuk_processed_annotated_july23.RDS')
## Add L029 and L058
leuk_srat.new = merge_seurat_objects(leuk_srat,s.list[['L029_TALL_D0_Blood']],keepAllGenes = F,genomeVersions = c('v38','v38'))
leuk_srat.new = merge_seurat_objects(leuk_srat.new,s.list[['L058_AML_D0_Blood']],keepAllGenes = F,genomeVersions = c('v38','v38'))
leuk_srat.new = standard_clustering(leuk_srat.new)


leuk_srat.new$donorID[is.na(leuk_srat.new$donorID)] = ann$donorID[match(leuk_srat.new$cellID[is.na(leuk_srat.new$donorID)],ann$cellID)]

ann$cellID3 = ifelse(ann$donorID %in% c('L029','L058'),ann$cellID,ann$cellID2)
ann$inLeuk = (ann$cellID3 %in% leuk_srat.new$cellID)
table(ann$cellID3 %in% leuk_srat.new$cellID)
table(leuk_srat.new$cellID %in% ann$cellID3,leuk_srat.new$donorID)

table(ann$annot[ann$inLeuk == F],ann$sampleID[ann$inLeuk == F])



leuk_srat.new$brams_annot = leuk_srat.new$finalAnn_july23
leuk_srat.new$brams_annot[leuk_srat.new$cellID %in% ann$cellID3] =  ann$annot[leuk_srat.new$cellID %in% ann$cellID3]
leuk_srat.new$finalAnn_july23[is.na(leuk_srat.new$finalAnn_july23)] = leuk_srat.new$brams_annot[is.na(leuk_srat.new$finalAnn_july23)]
tmp = as.data.frame(table(leuk_srat.new$brams_annot,leuk_srat.new$finalAnn_july23,leuk_srat.new$donorID))
tmp = tmp[tmp$Freq >0,]

leuk_srat.new$brams_annot[is.na(leuk_srat.new$brams_annot)] = 'NA'
DimPlot(leuk_srat.new,group.by = 'brams_annot',label = T,label.box = T,cols = col25[-6])+NoLegend()
table(is.na(leuk_srat.new$brams_annot))




## Combine SETBP1gl + HCA + LEUK
big.srat = merge_seurat_objects(combSrat,leuk_srat,keepAllGenes = F,genomeVersions = c('v38','v38'))
#big.srat$donorID[is.na(big.srat$donorID)] = '15M'
#big.srat$dataset2[is.na(big.srat$dataset2)] = '15M_leuk'

# Clustering with and withouth Harmony
big.srat_harm = standard_clustering(big.srat,runHarmony = T,harmonyVar = c('donorID'),doPlot = T)

big.srat = standard_clustering(big.srat,runHarmony = F)

# Some plotting
DimPlot(big.srat,cells.highlight = big.srat$cellID[big.srat$donorID == 'L061'])
DimPlot(leuk_srat,cells.highlight = big.srat.sub$cellID[big.srat.sub$seurat_clusters == 13])

DimPlot(big.srat,group.by = 'annot',label = T,label.box = T,cols=c(col25[-6],col22))
DimPlot(leuk_srat,group.by = 'donorID',label = T,label.box = T) + NoLegend()
DimPlot(leuk_srat,group.by = 'annot',label = T,label.box = T) + NoLegend()

# remove BALL, AML and TALL leukaemic cells, and doublets
cells_toRemove = leuk_srat$cellID[leuk_srat$donorID != 'L061' & leuk_srat$annot %in% c('BALL','AML','TALL','doublets')]

DimPlot(big.srat,cells.highlight = cells_toRemove)
leuk_srat$cellToRemove = (leuk_srat$cellID %in% cells_toRemove)
leuk_srat.sub = subset(leuk_srat,subset = cellID %in% leuk_srat$cellID[!leuk_srat$cellID %in% cells_toRemove])
big.srat.sub = subset(big.srat,subset = cellID %in% big.srat$cellID[!big.srat$cellID %in% cells_toRemove])

DimPlot(big.srat.sub,group.by = 'seurat_clusters',label = T,label.box = T,cols=c(col25[-6],col22,col25[-6])) + NoLegend()
big.srat.sub$dataset2[big.srat.sub$donorID == 'L061'] = 'L061'
big.srat.sub_harm$assay[big.srat.sub_harm$dataset2 %in% c('BJ9','GOSH84')] = "scRNA 10x5'v2(DUAL)"
big.srat.sub_harm$assay[big.srat.sub_harm$dataset2 %in% c('HCA')] = "sc3prime"

# Recluster srat objects
big.srat.sub = standard_clustering(big.srat.sub)
big.srat.sub_harm = standard_clustering(big.srat.sub,runHarmony = T,harmonyVar = c('dataset2'),doPlot = T) 
DimPlot(big.srat.sub_harm,group.by = 'seurat_clusters',label = T,label.box = T,cols=c(col25[-6],col22,col25[-6]),repel = T) 
DimPlot(big.srat.sub_harm,group.by = 'finalAnn_july23',label = T,label.box = T,cols=c(col25[-6],col22,col25[-6]),repel = T,label.size = 3) + NoLegend()
big.srat.sub$annot2 = big.srat.sub_harm$annot2[match(big.srat.sub$cellID,big.srat.sub_harm$cellID)]

big.srat.sub_harm$annot[big.srat.sub_harm$annot == '24'] = 'T_CD8'


DimPlot(big.srat.sub_harm,cells.highlight = big.srat.sub_harm$cellID[big.srat.sub_harm$seurat_clusters %in% c(7,37,0) & big.srat.sub_harm$annot2 == 'neutrophil'])
DimPlot(big.srat.sub_harm,cells.highlight = big.srat.sub_harm$cellID[big.srat.sub_harm$annot2 == 'BALL'])
DimPlot(leuk_srat,cells.highlight = big.srat.sub_harm$cellID[big.srat.sub_harm$seurat_clusters %in% c(25)])
table(big.srat.sub_harm$donorID[big.srat.sub_harm$seurat_clusters == 25])

big.srat.sub_harm$annot2 = big.srat.sub_harm$annot
big.srat.sub_harm$annot2[big.srat.sub_harm$annot2 %in% c('28','29','31','normalT')] = 'T_CD4'
big.srat.sub_harm$annot2[big.srat.sub_harm$annot2 %in% c('32')] = 'BALL'
big.srat.sub_harm$annot2[big.srat.sub_harm$annot2 %in% c('40')] = 'AML'
big.srat.sub_harm$annot2[big.srat.sub_harm$annot2 %in% c('33','23','3','41')] = 'B.cell'
big.srat.sub_harm$annot2[big.srat.sub_harm$annot2 %in% c('38')] = 'Mono_CD16'
big.srat.sub_harm$annot2[big.srat.sub_harm$annot2 %in% c('39') & big.srat.sub_harm$seurat_clusters == 33] = 'Plasma.cell'
big.srat.sub_harm$annot2[big.srat.sub_harm$annot2 %in% c('B?')] = as.character(big.srat.sub_harm$seurat_clusters[big.srat.sub_harm$annot2 %in% c('B?')])
big.srat.sub_harm$annot2[big.srat.sub_harm$annot2 %in% c('1','10','11','6','9','5','36','21','15','35','32','25')] = 'T_CD4'
big.srat.sub_harm$annot2[big.srat.sub_harm$annot2 %in% c('18','2','8')] = 'T_CD8'
big.srat.sub_harm$annot2[big.srat.sub_harm$annot2 %in% c('20')] = 'T_CD4_reg'
big.srat.sub_harm$annot2[big.srat.sub_harm$annot2 %in% c('37')] = 'pDC'
big.srat.sub_harm$annot2[big.srat.sub_harm$annot2 %in% c('30','17','4','40')] = 'NK'
big.srat.sub_harm$annot2[big.srat.sub_harm$annot2 %in% c('14')] = 'MAIT_Tcell'
big.srat.sub_harm$annot2[big.srat.sub_harm$annot2 %in% c('22')] = 'MK'
big.srat.sub_harm$annot2[big.srat.sub_harm$annot2 %in% c('39')] = 'Plasma.cell'
big.srat.sub_harm$annot2[big.srat.sub_harm$annot2 %in% c('38')] = 'unknown'
big.srat.sub_harm$annot2[big.srat.sub_harm$seurat_clusters %in% c(34,44)] = 'neutrophil'
big.srat.sub_harm$annot2[big.srat.sub_harm$seurat_clusters %in% c(7)] = 'Mono_CD14'
big.srat.sub_harm$annot2[big.srat.sub_harm$seurat_clusters %in% c(28)] = 'TALL?'
big.srat.sub_harm$annot2[big.srat.sub_harm$annot2 %in% c('TALL','27')] = 'unknown'
big.srat.sub_harm$annot2[big.srat.sub_harm$seurat_clusters == 40] = 'neutrophil'
big.srat.sub_harm$annot2[big.srat.sub_harm$seurat_clusters == 34] = 'c34'
big.srat.sub_harm$annot2[big.srat.sub_harm$seurat_clusters == 17] = 'Mono_CD14'

big.srat.sub_harm$annot2[big.srat.sub_harm$seurat_clusters %in% c(7,37,0) & big.srat.sub_harm$annot2 == 'neutrophil'] = 'Mono_CD14'
big.srat.sub_harm$annot2[big.srat.sub_harm$seurat_clusters %in% c(7,37,0) & big.srat.sub_harm$annot2 == 'BALL'] = 'Mono_CD14'
big.srat.sub_harm$annot2[big.srat.sub_harm$seurat_clusters %in% c(31,25)] = 'T_CD4'
big.srat.sub_harm$annot2[big.srat.sub_harm$seurat_clusters %in% c(3)] = 'B.cell'
big.srat.sub_harm$annot2[big.srat.sub_harm$seurat_clusters %in% c(35) & big.srat.sub_harm$annot2 == 'MK'] = 'Mono_CD14'
big.srat.sub_harm$annot2[big.srat.sub_harm$seurat_clusters %in% c(35) & big.srat.sub_harm$seurat_clusters %in% c(7,0) == 'BALL'] = 'Mono_CD14'
big.srat.sub_harm$annot2[big.srat.sub_harm$seurat_clusters %in% c(28) & big.srat.sub_harm$annot2 == 'BALL'] = 'pDC'
big.srat.sub_harm$annot2[big.srat.sub_harm$seurat_clusters %in% c(13) & big.srat.sub_harm$annot2 == 'BALL'] = 'Mono_CD16'

DimPlot(big.srat.sub_harm,cells.highlight = big.srat.sub_harm$cellID[big.srat.sub_harm$seurat_clusters == 35 & big.srat.sub_harm$annot2 == 'MK'])
big.srat.sub_harm$annot2[big.srat.sub_harm$seurat_clusters == 24 & big.srat.sub_harm$annot2 == 'BALL'] = 'DC2'
big.srat.sub_harm$annot2[big.srat.sub_harm$seurat_clusters == 24 & big.srat.sub_harm$annot2 == 'T_CD4'] = 'unknown'
big.srat.sub_harm$annot2[big.srat.sub_harm$annot2 == 'BALL' & big.srat.sub_harm$cellID %in% leuk_srat$cellID[leuk_srat$seurat_clusters == 11]] = 'MDS'
table(big.srat.sub_harm$donorID[big.srat.sub_harm$annot2 %in% c('MDS')]) 


big.srat.sub_harm$annot2[big.srat.sub_harm$annot2 %in% c('12') & big.srat.sub_harm$cellID %in% WhichCells(big.srat.sub_harm,expression = TRDV2 > 0.55)] = 'T_gd'
big.srat.sub_harm$annot2[big.srat.sub_harm$annot2 %in% c('12') & !big.srat.sub_harm$cellID %in% WhichCells(big.srat.sub_harm,expression = TRDV2 > 0.55)] = 'T_CD8'


library(SoupX)
markers = quickMarkers(toc=big.srat.sub_harm@assays$RNA@counts,clusters = big.srat.sub_harm$annot2)
markers_sratClusters = quickMarkers(toc=big.srat.sub_harm@assays$RNA@counts,clusters = big.srat.sub_harm$seurat_clusters)

## Finalize annotation
DimPlot(leuk_srat,group.by = 'finalAnn_july23',label = T,label.box = T,repel = T,cols=c(col25[-6],col22))+NoLegend()
DimPlot(big.srat.sub,cells.highlight = big.srat.sub_harm$cellID[big.srat.sub_harm$annot2 == 'BALL' & big.srat.sub_harm$cellID %in% leuk_srat$cellID[leuk_srat$seurat_clusters == 11]])
DimPlot(big.srat.sub_harm,cells.highlight = big.srat.sub_harm$cellID[big.srat.sub_harm$annot2 == '29'])


View(table(big.srat.sub_harm$annot,big.srat.sub_harm$annot2))
big.srat.sub_harm$finalAnn_july23 = big.srat.sub_harm$annot2
big.srat.sub_harm$finalAnn_july23[big.srat.sub_harm$finalAnn_july23 %in% c('26','29','31','AML','BALL','c34','DC.precursor','HSC_MEMP','Mast_cell','neutrophil','unknown','DC1')] = 'others'
big.srat.sub_harm$ann_july23 = paste0(big.srat.sub_harm$finalAnn_july23,':',big.srat.sub_harm$dataset2)
big.srat.sub_harm@meta.data = big.srat.sub_harm@meta.data[,colnames(big.srat.sub_harm@meta.data) != 'ann2']

outDir = '~/lustre_mt22/SETBP1/Results/1_setbp1QC/jul23'
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}
saveRDS(big.srat.sub_harm,'~/lustre_mt22/SETBP1/Results/1_setbp1QC/jul23/setbp1_HCA_leuk_processed_annotated_harm_july23.RDS')

big.srat.sub$finalAnn_july23 = big.srat.sub_harm$finalAnn_july23[match(big.srat.sub$cellID,big.srat.sub_harm$cellID)]
leuk_srat$finalAnn_july23 = leuk_srat$cellToRemove
leuk_srat$finalAnn_july23[leuk_srat$cellToRemove == F] = big.srat.sub$finalAnn_july23[match(leuk_srat$cellID[leuk_srat$cellToRemove == F],big.srat.sub_harm$cellID)]

# In this big.srat.sub object, ideally, only SETBP1-germlineMUT cells + HCA (adult normal) + LEUK_L061 (normal cells + mutated leuk) + LEUK_others (normal cells only)
saveRDS(big.srat.sub,'~/lustre_mt22/SETBP1/Results/1_setbp1QC/jul23/setbp1_HCA_leuk_processed_annotated_july23.RDS')
saveRDS(leuk_srat,'~/lustre_mt22/SETBP1/Results/1_setbp1QC/jul23/leuk_processed_annotated_july23.RDS')




### Try clustering without HCA cells - ie. remove 3prime cells to see if batch effect is improved
srat.noHCA = subset(big.srat.sub,subset = cellID %in% big.srat.sub$cellID[big.srat.sub$dataset2 != 'HCA'])
srat.noHCA = standard_clustering(srat.noHCA)

DimPlot(srat.noHCA,group.by = 'finalAnn_july23',label = T,label.box = T,repel = T,cols=col25[-6])+NoLegend()

## Quite heavy on the batch effect.... so lets do DEGs all in one go!








Idents(big.srat.sub_harm) = big.srat.sub_harm$ann2
DotPlot(big.srat.sub_harm,idents = unique(big.srat.sub_harm$ann2[grepl('T_|Mono|BALL',big.srat.sub_harm$ann2)]),
        group.by = 'ann2',features = c('PARP1','TP53','ANXA2','BAD','MTRNR2L12','NKTR'))

DotPlot(big.srat.sub_harm,group.by = 'dataset2',features = c('PARP1','TP53','BRCA1','BRCA2','ANXA2'))
leuk_srat$ann = paste0(leuk_srat$Phase,':',leuk_srat$donorID)

FeaturePlot(big.srat.sub,'MTRNR2L12')


DimPlot(big.srat.sub_harm,cells.highlight = big.srat.sub$cellID[big.srat.sub$seurat_clusters==17])
FeaturePlot(big.srat.sub_harm,'NKTR')
DotPlot(big.srat.sub_harm,features = c('MTRNR2L12'),group.by = 'assay')
DimPlot(leuk_srat,cells.highlight = leuk_srat$cellID[leuk_srat$seurat_clusters == 13])
DimPlot(leuk_srat,cells.highlight = leuk_srat$cellID[leuk_srat$seurat_clusters == 13])





big.srat.sub$sex[big.srat.sub$donorID %in% c('L061','L057','L028','L014')] = 'M'
big.srat.sub$sex[!big.srat.sub$donorID %in% c('L061','L057','L028','L014')] = 'F'
big.srat.sub$ageGroup = ifelse(big.srat.sub$dataset2=='HCA','youngAdult','paediatric')
big.srat.sub$assay = ifelse(big.srat.sub$dataset2=='HCA','GEX3p','GEX5p')
big.srat.sub$group = big.srat.sub$dataset2
big.srat.sub$group[!big.srat.sub$group %in% c('BJ9','GOSH84')] = 'healthyControl'




#------------------------------------
#### Identify differentially expressed genes #####

outDir = '~/lustre_mt22/SETBP1/Results/3_pseudobulk_DESeq2/july23/'
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}
setwd(outDir)

plotDir = file.path(outDir,'plots')
if(!dir.exists(plotDir)){
  dir.create(plotDir,recursive = T)
}


srat = big.srat.sub
srat$finalAnn_broad = srat$finalAnn_july23
table(srat$finalAnn_broad)

# Remove cycling cells
srat = subset(srat,subset = cellID %in% srat$cellID[srat$Phase == 'G1'])
table(srat$dataset2,srat$finalAnn_broad)

###### compute pseudobulk for each cell type #####
# Loop over cell-type

out = list()
for(tgtIdx in c(1:(length(unique(srat@meta.data$finalAnn_broad))))){
  tgtCell = unique(srat@meta.data$finalAnn_broad)[tgtIdx] 
  message(sprintf("Consider cell type %s",tgtCell))
  if(tgtCell %in% c('NK?','NA','MDS','T_cell?','DC.precursor','others','unknown')){
    next
  }
  
  # Check we have a reasonable number of counts from both settings. ie Get number of cells for the relevant cell type and genenotype
  #nCellsGroup = table(factor(srat@meta.data$batch[srat$finalAnn_broad == tgtCell],levels=c('BJ9_1', 'GOSH84_1','GOSH84_2','HCA')))
  nCellsGroup = table(factor(srat@meta.data$dataset2[srat$finalAnn_broad == tgtCell],levels=c('BJ9', 'GOSH84','HCA','leuk','L061')))
  nCellsGroup = nCellsGroup[nCellsGroup >20]
  if(!all(nCellsGroup>20)){
    message(sprintf('There are too few cells of type %s.  Skipping...',tgtCell))
    next
  }
  #Check how many from individual donors
  nCells = table(srat@meta.data$donorID[srat$finalAnn_broad == tgtCell])
  if(sum(nCells>20)<3){
    message(sprintf("Too few effective replicates.  Skipping..."))
    print(nCells)
    next
  }
  message("Found the following number of high confidence cells")
  print(nCells)
  
  
  #OK, we're going ahead, create the needed objects
  
  cellID_toKeep = srat$cellID[srat$finalAnn_broad == tgtCell]
  cnt_mtx = as.matrix(srat@assays$RNA@counts[,colnames(srat@assays$RNA@counts) %in% cellID_toKeep])
  rownames(cnt_mtx) = geneMap$ensID[match(rownames(cnt_mtx),geneMap$geneSym)]
  
  mDat = srat@meta.data[colnames(cnt_mtx),c('cellID','group','assay','sex','ageGroup','donorID')]
  #mDat$donorID[mDat$donorID == 'GOSH84'] = mDat$batch[mDat$donorID == 'GOSH84']
  mDat$group = factor(mDat$group,levels = c('healthyControl','BJ9','GOSH84'))
  table(mDat$donorID)
  
  pdf(file.path(plotDir,paste0(tgtCell,'.pdf')))
  #Shrink the LFC estimates
  res = tryCatch({compareCell_simplified(toc = cnt_mtx,
                                                    mDat = mDat,
                                                    coords = gns[rownames(cnt_mtx)[rownames(cnt_mtx) %in% names(gns)]],
                                                    cellTypeName=tgtCell,
                                                    formula='~ %s + assay + sex + age',
                                                    donorID='donorID',groupID='group',
                                                    geneMap = geneMap)
  }, error = function(e){
    print(e)
    message('\nTrying removing age')
    res = compareCell_simplified(toc = cnt_mtx,
                                 mDat = mDat,
                                 coords = gns[rownames(cnt_mtx)[rownames(cnt_mtx) %in% names(gns)]],
                                 cellTypeName=tgtCell,
                                 formula='~ %s + assay + sex',
                                 donorID='donorID',groupID='group',
                                 geneMap = geneMap)
    return(res)
  }, error = function(e){
    print(e)
    message('\nTrying removing sex')
    res = compareCell_simplified(toc = cnt_mtx,
                                 mDat = mDat,
                                 coords = gns[rownames(cnt_mtx)[rownames(cnt_mtx) %in% names(gns)]],
                                 cellTypeName=tgtCell,
                                 formula='~ %s + assay',
                                 donorID='donorID',groupID='group',
                                 geneMap = geneMap)
    return(res)
  },finally = {message('Successfully performed compareCell_simplified!')})
  
  out[[tgtCell]] = res
  
  # out[[tgtCell]] = compareCell_simplified(toc = cnt_mtx,
  #                                         mDat = mDat,
  #                                         coords = gns[rownames(cnt_mtx)[rownames(cnt_mtx) %in% names(gns)]],
  #                                         cellTypeName=tgtCell,
  #                                         formula='~ %s + assay',
  #                                         donorID='donorID',groupID='group',
  #                                         geneMap = geneMap)
  dev.off()
}


saveRDS(out,file.path(outDir,'DESeq2_byCT_byCondition.RDS'))


all_degs = lapply(seq(1:length(out)),FUN = function(i){
  ct = names(out)[i]
  deg = out[[i]][['all_de']]
  if(nrow(deg) > 0){
    deg$celltype = ct  
    return(deg)
  }else{
    return()
  }
  
})

all_degs = all_degs[sapply(seq(1:length(all_degs)),function(x){!is.null(all_degs[[x]])})]
all_degs = do.call(rbind,all_degs)

## Main DEGs between the conditions
all_degs.main = all_degs[all_degs$contrast %in% c('group_BJ9_vs_healthyControl', 'group_GOSH84_vs_healthyControl'),]
all_degs.main$DE = ((all_degs.main$cellFrac_g1 >= 15 | all_degs.main$cellFrac_g2 >= 15))
all_degs.main$direction = ifelse(all_degs.main$log2FoldChange > 0, 'down_HCA','up_HCA')

write.csv(all_degs.main,'~/lustre_mt22/SETBP1/Results/3_pseudobulk_DESeq2/july23/DESeq2_byCT_byCondition_allDEGs_mainContrast.csv') #####
all_degs.main = read.csv('~/lustre_mt22/SETBP1/Results/3_pseudobulk_DESeq2/july23/DESeq2_byCT_byCondition_allDEGs_mainContrast.csv')

## Technical DEGs (assay and sex)
all_degs.tech = all_degs[!all_degs$contrast %in% c('group_BJ9_vs_healthyControl', 'group_GOSH84_vs_healthyControl'),]
all_degs.tech$direction = ifelse(all_degs.tech$log2FoldChange > 0,'up','down')

write.csv(all_degs.tech,'~/lustre_mt22/SETBP1/Results/3_pseudobulk_DESeq2/july23/DESeq2_byCT_byCondition_allDEGs_techContrast.csv') #####
all_degs.tech = read.csv('~/lustre_mt22/SETBP1/Results/3_pseudobulk_DESeq2/july23/DESeq2_byCT_byCondition_allDEGs_techContrast.csv')


all_degs.tech.sum = all_degs.tech %>% group_by(contrast,geneSym,direction) %>% summarise(nCT = n_distinct(ct))


##------ Analyse main DEGs -------###
df = all_degs.main[all_degs.main$DE == T,] %>% group_by(geneSym,ct,direction) %>% summarise(nContrast = n_distinct(contrast))

df2 = all_degs.main[all_degs.main$DE == T,] %>% group_by(geneSym,contrast,direction) %>% summarise(nCT = n_distinct(ct))
df2$ensID = geneMap$ensID[match(df2$geneSym,geneMap$geneSym)]
df2 = annotateGenes(df2,geneMap = geneMap)
df2$geneSym = geneMap$geneSym[match(df2$ensID,geneMap$ensID)]
hist(df2$nCT[df2$contrast == 'dataset2_GOSH84_vs_HCA'])





##-------------------------------------------------------##
## Number of DEGs per celltype as a function of #cells ####
##-------------------------------------------------------##
# I want x = number of tgt_cells; y = number of ref_cells; dot_size = # DEGs
all_degs = all_degs.main 
nCell = as.data.frame(table(big.srat.sub$group[big.srat.sub$Phase == 'G1'],big.srat.sub$finalAnn_july23[big.srat.sub$Phase == 'G1']))
colnames(nCell) = c('dataset','celltype','nCell')

nDEG_summary = all_degs[all_degs$DE == T,] %>% group_by(celltype,direction,contrast) %>% summarise(nDEG = n_distinct(geneSym)) %>% 
  group_by(celltype,contrast) %>% mutate(nDEG_total = sum(nDEG))
ctr = nCell[nCell$dataset == 'healthyControl',]
nDEG_summary$nCell_ctrl =  ctr$nCell[match(nDEG_summary$celltype,ctr$celltype)]

gosh84 = nCell[nCell$dataset == 'GOSH84',]
nDEG_summary$nCell_gosh84 =  gosh84$nCell[match(nDEG_summary$celltype,gosh84$celltype)]

bj9 = nCell[nCell$dataset == 'BJ9',]
nDEG_summary$nCell_bj9 =  bj9$nCell[match(nDEG_summary$celltype,bj9$celltype)]


## calculate number of genes "expressed" in the cell type
# because 1 of the criteria for DEG was that it has to be expressed in >=X % of cells in at least 1 group of cells.
# and for the gene to be considered as expressed in a given cell, the raw count > 0.55
# I think it's reasonable to consider the number of genes "expressed" as also adhering to both of these criteria

nGeneExpressed_byCT = data.frame()
big.srat.sub$ann = paste0(big.srat.sub$finalAnn_july23,':',big.srat.sub$dataset2)
for(ct in unique(big.srat.sub$ann)){
  if(ct %in% c('?','DC.precursor','NA','NK?','Tcell?','others','DC1','unknown','HSC_MEMP:BJ9','others','MDS','MK:L061')){next}
  
  toc = big.srat.sub@assays$RNA@counts[,colnames(big.srat.sub@assays$RNA@counts) %in% big.srat.sub$cellID[big.srat.sub$ann == ct]]
  if(class(toc) == 'numeric'){
    next
  }
  tmp = data.frame(nCellExpressed = apply(toc,1,function(x){sum(x>0.55)}),
                   nCellTotal = ncol(toc),
                   ct = ct,
                   geneSym = rownames(toc))
  if(nrow(nGeneExpressed_byCT) == 0){
    nGeneExpressed_byCT = tmp
  }else{
    nGeneExpressed_byCT = rbind(nGeneExpressed_byCT,tmp)  
  }
}

nGeneExpressed_byCT$percExpr = 100*nGeneExpressed_byCT$nCellExpressed/ nGeneExpressed_byCT$nCellTotal
nGeneExpressed_byCT$celltype = sapply(strsplit(nGeneExpressed_byCT$ct,split = ':'),'[',1)
nGeneExpressed_byCT$dataset = sapply(strsplit(nGeneExpressed_byCT$ct,split = ':'),'[',2)
#nGeneExpressed_byCT$isExpr = ifelse(nGeneExpressed_byCT$nCellExpressed <2 | nGeneExpressed_byCT$percExpr < 15,F,T)
nGeneExpressed_byCT$isExpr = ifelse(nGeneExpressed_byCT$nCellExpressed <2 | nGeneExpressed_byCT$percExpr < 5,F,T)

nGeneExpressed_byCT_sum = nGeneExpressed_byCT %>% group_by(celltype,dataset,geneSym) %>% summarise(nExpr = sum(isExpr))
nGeneExpressed_byCT_sum_ctrl = nGeneExpressed_byCT_sum[nGeneExpressed_byCT_sum$dataset %in% c('HCA','L061','leuk'),]
nGeneExpressed_byCT_sum_ctrl = nGeneExpressed_byCT_sum_ctrl[nGeneExpressed_byCT_sum_ctrl$nExpr >0,]
nGeneExpressed_byCT_sum$isExpr_inCtrl = NA
for(ct in unique(nGeneExpressed_byCT_sum$celltype)){
  nGeneExpressed_byCT_sum$isExpr_inCtrl[nGeneExpressed_byCT_sum$celltype == ct] = (nGeneExpressed_byCT_sum$geneSym[nGeneExpressed_byCT_sum$celltype == ct] %in% nGeneExpressed_byCT_sum_ctrl$geneSym[nGeneExpressed_byCT_sum_ctrl$celltype == ct])
}
#nGeneExpressed_byCT_sum$isExpr_inHCA = nGeneExpressed_byCT_sum$geneSym %in% nGeneExpressed_byCT_sum_hca$geneSym[match(nGeneExpressed_byCT_sum$celltype,nGeneExpressed_byCT_sum_hca$celltype)]
nGeneExpressed_byCT_sum$isExpr_overall = ifelse(nGeneExpressed_byCT_sum$nExpr == 1 | nGeneExpressed_byCT_sum$isExpr_inCtrl ==T , T, F)

nGeneExpressed_byCT_sum2 = nGeneExpressed_byCT_sum %>% group_by(celltype, dataset) %>% summarise(nGeneExpr = sum(isExpr_overall))

nGeneExpressed_byCT_sum2_bj9 = nGeneExpressed_byCT_sum2[nGeneExpressed_byCT_sum2$dataset == 'BJ9',]
nGeneExpressed_byCT_sum2_gosh84 = nGeneExpressed_byCT_sum2[nGeneExpressed_byCT_sum2$dataset == 'GOSH84',]
nDEG_summary$nGeneExpressed_byCT_BJ9 = nGeneExpressed_byCT_sum2_bj9$nGeneExpr[match(nDEG_summary$celltype,nGeneExpressed_byCT_sum2_bj9$celltype)]
nDEG_summary$nGeneExpressed_byCT_GOSH84 = nGeneExpressed_byCT_sum2_gosh84$nGeneExpr[match(nDEG_summary$celltype,nGeneExpressed_byCT_sum2_gosh84$celltype)]
nDEG_summary$fracDEG = ifelse(grepl('BJ9',nDEG_summary$contrast),
                              nDEG_summary$nDEG/nDEG_summary$nGeneExpressed_byCT_BJ9,
                              nDEG_summary$nDEG/nDEG_summary$nGeneExpressed_byCT_GOSH84)
nDEG_summary$fracDEG_total = ifelse(grepl('BJ9',nDEG_summary$contrast),
                                    nDEG_summary$nDEG_total/nDEG_summary$nGeneExpressed_byCT_BJ9,
                                    nDEG_summary$nDEG_total/nDEG_summary$nGeneExpressed_byCT_GOSH84)
nDEG_summary$nCell_overall = ifelse(grepl('BJ9',nDEG_summary$contrast),
                                    nDEG_summary$nCell_ctrl + nDEG_summary$nCell_bj9,
                                    nDEG_summary$nCell_ctrl + nDEG_summary$nCell_gosh84)


ggplot(nDEG_summary[grepl('GOSH84',nDEG_summary$contrast),],
       aes(x=(nCell_gosh84 + nCell_ctrl),y=fracDEG_total))+
  geom_point(aes(col=celltype))+
  geom_smooth(method = 'lm')+
  scale_color_manual(values=col25)+
  theme_bw()


ggplot(nDEG_summary[grepl('BJ9',nDEG_summary$contrast),],
       aes(x=(nCell_bj9 + nCell_ctrl),y=fracDEG_total))+
  geom_point(aes(col=celltype))+
  geom_smooth(method = 'lm')+
  scale_color_manual(values=col25)+
  theme_bw()


ggplot(nDEG_summary,
       aes(x=nCell_overall,y=fracDEG_total,shape = contrast))+
  geom_smooth(aes(col=contrast),method = 'lm')+
  geom_point(aes(col=celltype),size=2.5)+
  scale_color_manual(values=col25)+
  theme_bw(base_size = 13)+theme(legend.title = element_blank())+
  xlab('# Cells') + ylab('Fraction of expressed genes being DEGs')


## Next part in comparing it with BEAT_AML is in script #14
