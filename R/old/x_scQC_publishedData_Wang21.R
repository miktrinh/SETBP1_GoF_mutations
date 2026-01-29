##--- Preprocess the Wang etal 2021 single-cell adult thyroid cancer dataset ---##
#     https://www.nature.com/articles/s41467-021-25771-5
setwd('~/SETBP1_GoF_mutations/')

outDir = "Results/04_public_scRNAseq"
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}

##----------------##
##   Libraries  ####
##----------------##
library(tidyverse)
library(Seurat)
source("R/utils/sc_basicQC.R")

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
skipSoup = T
scrubScoreMax = 0.5
scrubPath='./scrubletScores.tsv'
scPath="./strainedCounts"
doPlot=T
verbose = T
skipIfExists=F
keepMTCells=T
rho_max_limit = NULL


##----------------------------------##
##   Preprocessing snRNAseq data  ####
##----------------------------------##
### 1. Import cellranger output data
###    Run SoupX
###    Subset to keep only cells present in the original publications
###    Add cell labels (as published)

outDir_sub = file.path(outDir,'2601') 

if(!is.null(rho_max_limit)){
  plotDir = paste0(file.path(outDir_sub,'Wang21_rhoLim'),rho_max_limit,'_')
  outPath = paste0(file.path(outDir_sub,'Wang21__rhoLim'),rho_max_limit,'_')
}else{
  plotDir = paste0(file.path(outDir_sub,'Wang21__rhoLimNone_'))
  outPath = paste0(file.path(outDir_sub,'Wang21__rhoLimNone_'))
}

cleanSrat_fp = ifelse(keepMTCells,paste0(outPath,'_clean_withMTCells.RDS'),paste0(outPath,'_clean_noMTCells.RDS'))
if(file.exists(cleanSrat_fp) & skipIfExists){
  cleanSrat = readRDS(cleanSrat_fp)  
}else{
  if(!dir.exists(outDir_sub)){
    message(sprintf('Creating output directory'))
    dir.create(outDir_sub,recursive = T)
  }
  
  #dataDirs = list.dirs('~/lustre_mt22/Thyroid/Data/published_scRNAseq/Pu_etal_2021/GSE184362_RAW')
  #dataDirs = dataDirs[!grepl('checkpoints|GSE184362_RAW$|cleanCounts$',dataDirs)]
  dataDirs = list.dirs('Data/published_scRNAseq/Wang21')
  dataDirs = dataDirs[grepl('H\\d$',dataDirs)]
  names(dataDirs) = basename(dataDirs)
  
  dataDirs=dataDirs[file.exists(dataDirs)]
  dataDirs=dataDirs[sapply(dataDirs, function(x){length(list.files(x))>0})]
  print(n_distinct(dataDirs))
  
  cleanCountDir = NULL # Specify NULL to save scrublet and soupX results in the same dataDir folder
  metadata = NULL
  matchBy = NULL
  is10X=TRUE
  
  # Run basicQC
  message('\nPerforming scRNAseq QC...')
  
  QC.output = basicQC(dataDirs = dataDirs,maxMT = maxMT, 
                      minGenes=minGenes,minUMIs=minUMIs,maxBadFrac=maxBadFrac,numPCs=numPCs,
                      clusteringRes=clusteringRes,cleanCountDir=cleanCountDir,
                      skipScrub=skipScrub,skipSoup=skipSoup,scrubScoreMax=scrubScoreMax,scrubPath=scrubPath,
                      metadata=metadata,matchBy=matchBy,scPath=scPath,rho_max_limit=rho_max_limit,
                      outPath=outPath,skipIfExists=skipIfExists,
                      doPlot=doPlot,plotDir=plotDir,verbose=verbose,is10X=is10X)
  
  cleanSrat = QC.output[[1]]
  df.out = QC.output[[2]]
  
  write.csv(df.out,paste0(outPath,'qc_summary.csv'))
}

wang21 = readRDS(file.path(outDir,'2601','Wang21__rhoLimNone__clean_noMTCells.RDS'))
DimPlot(wang21)

##------------------------------------------##
##    Combine with SETBP1 and annotate    ####
##------------------------------------------##
combSrat = readRDS('Results/03_SETBP1_annotation/2505/SETBP1_SGS_clean_noMTCells_annot_2505.RDS')

setbp1_wang21 = merge_seurat_objects(combSrat,wang21,keepAllGenes = F,genomeVersions = c('v38','v38'))
setbp1_wang21 = standard_clustering(setbp1_wang21)
setbp1_wang21$dataset = NA
setbp1_wang21$dataset[setbp1_wang21$donorID=='L061'] = 'MDS'
setbp1_wang21$dataset[setbp1_wang21$donorID %in% c('GOSH084','BJ9','BJ111','BJ112','BJ113','BJ114')] = 'SETBP1'
setbp1_wang21$dataset[setbp1_wang21$cellID %in% wang21$cellID] = 'wang21'
DimPlot(setbp1_wang21,group.by = 'annot',label = T,repel = T,label.box = T,label.size = 3.3,cols = c(col25,pal34H))+NoAxes() + NoLegend()+
  theme(panel.border = element_rect(fill=F,color = 'black',linewidth = 1))
setbp1_wang21$annot[setbp1_wang21$dataset == 'wang21'] = '?'
table(setbp1_wang21$annot[setbp1_wang21$seurat_clusters %in% c(44)])
setbp1_wang21$annot[setbp1_wang21$seurat_clusters %in% c(1,9,47,37,14) & setbp1_wang21$annot == '?'] = 'B.cell'
setbp1_wang21$annot[setbp1_wang21$seurat_clusters %in% c(51) & setbp1_wang21$annot == '?'] = 'Plasma_cell'
setbp1_wang21$annot[setbp1_wang21$seurat_clusters %in% c(40) & setbp1_wang21$annot == '?'] = 'DC2'
setbp1_wang21$annot[setbp1_wang21$seurat_clusters %in% c(31) & setbp1_wang21$annot == '?'] = 'Mono_CD16'
setbp1_wang21$annot[setbp1_wang21$seurat_clusters %in% c(8,6) & setbp1_wang21$annot == '?'] = 'Mono_Mac'
setbp1_wang21$annot[setbp1_wang21$seurat_clusters %in% c(2,35,24,41,16,11,0,21,29,7,38,5,12,15,3,46) & setbp1_wang21$annot == '?'] = 'T_cell'
setbp1_wang21$annot[setbp1_wang21$seurat_clusters %in% c(39) & setbp1_wang21$annot == '?'] = 'T_MAIT'
setbp1_wang21$annot[setbp1_wang21$seurat_clusters %in% c(10,20) & setbp1_wang21$annot == '?'] = 'NK_T'
setbp1_wang21$annot[setbp1_wang21$seurat_clusters %in% c(56,4,34) & setbp1_wang21$annot == '?'] = 'NK'

setbp1_wang21$donorID[is.na(setbp1_wang21$donorID)] = as.character(setbp1_wang21$orig.ident[is.na(setbp1_wang21$donorID)])
setbp1_wang21$donorID[is.na(setbp1_wang21$donorID)] = as.character(setbp1_wang21$orig.ident[is.na(setbp1_wang21$donorID)])

setbp1_wang21$assay[setbp1_wang21$dataset == 'wang21'] = "scRNA 10x5'v2(DUAL)"
setbp1_wang21$sex[setbp1_wang21$donorID %in% c('GSM5160432')] = 'male'
setbp1_wang21$sex[setbp1_wang21$donorID %in% c('GSM5160434','GSM5160435')] = 'female'



## sub-clustering B-cells (from all donor)
b_cells = subset(setbp1_wang21,subset = cellID %in% setbp1_wang21$cellID[setbp1_wang21$annot=='B.cell'])
b_cells = standard_clustering(b_cells)
b_cells = FindClusters(b_cells,resolution = 0.2)
b_cells$donorID[is.na(b_cells$donorID)] = b_cells$orig.ident[is.na(b_cells$donorID)]
DimPlot(b_cells,group.by = 'donorID',label = T,repel = T,label.box = T,cols = col25)
DimPlot(setbp1_wang21,cells.highlight = b_cells$cellID[b_cells$donorID == 24])
library(SoupX)
qm  = quickMarkers(b_cells@assays$RNA@counts,b_cells$seurat_clusters)


## sub-clustering MonoMac (from all donor)
monomac = subset(setbp1_wang21,subset= cellID %in% setbp1_wang21$cellID[setbp1_wang21$annot %in% c('DC1','DC2','?','Mono_Mac','Mono_CD16','pDC','41','50') & setbp1_wang21$assay != "snRNA 10x5'v2" & setbp1_wang21$donorID != 'L061'])
monomac = standard_clustering(monomac)
DimPlot(monomac,group.by = 'annot',label = T,repel = T,label.box = T,cols = col25)
FeaturePlot(monomac,'CLEC10A')
monomac$annot[monomac$seurat_clusters==22] = 'DC1'
monomac$annot[monomac$seurat_clusters==8] = 'DC2'
monomac$annot[monomac$seurat_clusters==11] = 'pDC'
monomac$annot[monomac$seurat_clusters==17] = 'Mono_CD16'
monomac$annot[monomac$seurat_clusters==21] = 'DC2'
DimPlot(setbp1_wang21,cells.highlight = monomac$cellID[monomac$seurat_clusters == 10])
qm  = quickMarkers(monomac@assays$RNA@counts,monomac$seurat_clusters)


setbp1_wang21$annot[setbp1_wang21$cellID %in% monomac$cellID] = monomac$annot[match(setbp1_wang21$cellID[setbp1_wang21$cellID %in% monomac$cellID],monomac$cellID)]


## Adjust annotation
setbp1_wang21$annot = setbp1_wang21$annot_oct24
setbp1_wang21$broadAnno = as.character(setbp1_wang21$annot)
setbp1_wang21$broadAnno[setbp1_wang21$broadAnno %in% c('T_CD8','T_cell','T_gd','T_MAIT')] = 'T_cells'
setbp1_wang21$broadAnno[setbp1_wang21$broadAnno %in% c('NK','NK_T')] = 'NK_T'
setbp1_wang21$broadAnno[setbp1_wang21$broadAnno %in% c('DC1','DC2')] = 'DCs'
setbp1_wang21$broadAnno[setbp1_wang21$broadAnno %in% c('Mono_Mac','Mono_CD16')] = 'MonoMac'
setbp1_wang21$broadAnno[setbp1_wang21$broadAnno %in% c('41','50','doublets','?')] = 'others'
#setbp1_wang21$broadAnno[setbp1_wang21$broadAnno %in% c('MDS?','teratoma')] = 'Tumour'


saveRDS(setbp1_wang21,'SETBP1_Wang21_combined_2410.RDS')








## General assessment of cell type distribution / proportion of cells
# Exclude BJ9 teratoma samples

df = setbp1_wang21@meta.data[!(setbp1_wang21$tissue == 'Tumour' & setbp1_wang21$donorID == 'BJ9') ,]
table(df$broadAnno)
df = df %>% filter(broadAnno !='others') %>% 
  group_by(donorID,orig.ident,broadAnno) %>% summarise(nCell = n()) %>% 
  group_by(donorID,orig.ident) %>% mutate(totalCell = sum(nCell),
                                          frac = nCell/totalCell)
df$condition = ifelse(df$donorID %in% c('BJ114','GSM5160432','GSM5160434','GSM5160435'),'normal',
                      ifelse(df$donorID == 'L061','MDS',
                             ifelse(df$donorID %in% c('GOSH084','BJ113'),'adj','hotspot')))
df$condition = factor(df$condition,c('normal','adj','hotspot','MDS'))
ggplot(df,aes(orig.ident,frac,fill=broadAnno))+
  geom_col()+
  scale_fill_manual(values = col25)+
  scale_y_continuous(breaks = c(0,0.5,1),labels = c(0,0.5,1))+
  facet_grid(condition~donorID,scales = 'free_x',space = 'free_x')+
  theme_classic(base_size = 14)+
  theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1),
        axis.text = element_text(colour='black'),
        axis.line = element_blank(),strip.background = element_blank(),
        panel.border = element_rect(fill=F,colour='black'))+xlab('')+ylab('Fraction of cells')


## Cycling cells
df = setbp1_wang21@meta.data[!(setbp1_wang21$tissue == 'Tumour' & setbp1_wang21$donorID == 'BJ9') & (setbp1_wang21$annot != 'MDS?'),]
df = df %>% filter(!broadAnno %in% c('others','strange_MonoMac','Tumour')  & donorID != 'L061' & assay == "scRNA 10x5'v2(DUAL)") %>% 
  group_by(donorID,Phase) %>% summarise(nCell = n()) %>% 
  group_by(donorID) %>% mutate(totalCell = sum(nCell),
                               frac = nCell/totalCell)
df$condition = ifelse(df$donorID %in% c('BJ114','GSM5160432','GSM5160434','GSM5160435'),'normal',
                      ifelse(df$donorID == 'L061','MDS',
                             ifelse(df$donorID %in% c('GOSH084','BJ113'),'adj','hotspot')))
df$condition = factor(df$condition,c('normal','adj','hotspot','MDS'))
df$Phase = factor(df$Phase,c('G1','S','G2M'))
ggplot(df,aes(condition,frac))+
  geom_boxplot(aes(fill=condition))+
  geom_point()+
  scale_fill_manual(values = col25[c(1,3,2)])+
  scale_y_continuous(limits = c(0,1),breaks = c(0,0.5,1),labels = c(0,0.5,1))+
  facet_grid(~Phase,scales = 'free_x',space = 'free_x')+
  theme_classic(base_size = 14)+
  theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1),
        axis.text = element_text(colour='black'),
        axis.line = element_blank(),strip.background = element_blank(),
        panel.border = element_rect(fill=F,colour='black'))+xlab('')+ylab('Fraction of cells')







##================================================================================================================================================================================================##
##    Define SETBP1 gene signatures   ####

outDir = '~/lustre_mt22/SETBP1/Results/3_pseudobulk_DESeq2/oct24_v2'
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}
setwd(outDir)


##-------------------##
##   Libraries     ####
##-------------------##

library(Seurat)
library(GenomicFeatures)
library(DESeq2)
library(ComplexHeatmap)
library(reshape2)
library(zoo)
library(tidyverse)
library(RColorBrewer)
library(readxl)
source("~/lustre_mt22/generalScripts/utils/misc.R")
source("~/lustre_mt22/generalScripts/utils/sc_utils.R")
source("~/lustre_mt22/generalScripts/utils/pseudobulk.R")

##-----------------------##
##        Params          #
##-----------------------##
#Define genomic coordinates
gtf = '/nfs/srpipe_references/downloaded_from_10X/refdata-gex-GRCh38-2020-A/genes/genes.gtf'
txdb = makeTxDbFromGFF(gtf)
gns = genes(txdb)

geneMap = read.delim('/lustre/scratch125/casm/team274sb/mt22/SETBP1/Data/SETBP1/cellranger612_count_42668_CG_SB_NB11528275_GRCh38-2020-A/filtered_feature_bc_matrix/features.tsv.gz',header = F)
colnames(geneMap) = c('ensID','geneSym','GEX')
geneMap$geneSym[duplicated(geneMap$geneSym)] = paste0(geneMap$geneSym[duplicated(geneMap$geneSym)],'.1')
geneMap$geneSym = gsub('_','-',geneMap$geneSym)
geneMap$chr = as.character(seqnames(gns[match(geneMap$ensID,gns$gene_id)]))


kMerge=11
keepCylcingCells=T
remove_cyclingCells=F
ageMatch = T



## Import srat object

setbp1_wang21 = readRDS('~/lustre_mt22/SETBP1/Results/x_published_scRNAseq_preprocessing/SETBP1_Wang21_combined_2410.RDS')

## Define cell groups for DEG
setbp1_wang21$group = ifelse(setbp1_wang21$annot %in% c('B.cell','DC1','DC2','MK','Mono_CD16','Mono_Mac','NK','NK_T','pDC','Plasma_cell'),setbp1_wang21$annot,
                        ifelse(setbp1_wang21$annot %in% c('T_CD8','T_cell','T_gd','T_MAIT'),'T.cell','others'))
table(setbp1_wang21$broadAnno)

setbp1_wang21$mutationGroup = ifelse(setbp1_wang21$condition == 'unaffected' | setbp1_wang21$donorID %in% c(24,25,26),'normal',
                                ifelse(setbp1_wang21$donorID %in% c('GOSH084','BJ113'),'adjacent','hotspot'))




##------- Using pseudobulk DESeq2 ----------####

formula = '~ %s + dataset'

outDir = file.path(outDir,'mutationGroup_dataset')
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}
setwd(outDir)
plotDir = outDir

###### compute pseudobulk for each cell type #####
# Loop over cell-type
out = list()
for(tgtIdx in c(11:(length(unique(setbp1_wang21$group))))){
  tgtCell = unique(setbp1_wang21$group)[tgtIdx] 
  message(sprintf("Consider cell type %s",tgtCell))
  if(tgtCell %in% c('others','Plasma_cell','MK')){
    next
  }
  
  # Check we have a reasonable number of counts from both settings. ie Get number of cells for the relevant cell type and genenotype
  nCellsGroup = table(factor(setbp1_wang21$mutationGroup[setbp1_wang21$group == tgtCell & setbp1_wang21$donorID != 'L061' & setbp1_wang21$assay != "snRNA 10x5'v2"],
                             levels=c('adjacent', 'hotspot','normal')))
  nCellsGroup = nCellsGroup[nCellsGroup >30]
  if(!all(nCellsGroup>30)){
    message(sprintf('There are too few cells of type %s.  Skipping...',tgtCell))
    next
  }
  #Check how many from individual donors
  nCells = table(setbp1_wang21$donorID[setbp1_wang21$group == tgtCell & setbp1_wang21$donorID != 'L061' & setbp1_wang21$assay != "snRNA 10x5'v2"])
  if(sum(nCells>30)<3){
    message(sprintf("Too few effective replicates.  Skipping..."))
    print(nCells)
    next
  }
  message("Found the following number of high confidence cells")
  print(nCells)
  
  # remove individuals with < 30 cells
  donorID_toKeep = names(nCells[nCells >= 30])
  
  
  #OK, we're going ahead, create the needed objects
  
  cellID_toKeep = setbp1_wang21$cellID[setbp1_wang21$group == tgtCell & setbp1_wang21$donorID %in% donorID_toKeep & setbp1_wang21$donorID != 'L061' & setbp1_wang21$assay != "snRNA 10x5'v2"]
  cnt_mtx = setbp1_wang21@assays$RNA@counts[,colnames(setbp1_wang21@assays$RNA@counts) %in% cellID_toKeep]  
  rownames(cnt_mtx) = geneMap$ensID[match(rownames(cnt_mtx),geneMap$geneSym)]
  
  mDat = setbp1_wang21@meta.data[colnames(cnt_mtx),c('donorID','mutationGroup','sex','dataset','cellID')]
  mDat$group = mDat$mutationGroup
  mDat$group = factor(mDat$group,c('normal','adjacent','hotspot'))
  
  pdf(file.path(plotDir,paste0(tgtCell,'.pdf')))
  out[[paste0(tgtCell)]] = compareCell_simplified(toc = cnt_mtx,
                                                  mDat = mDat,
                                                  coords = gns[rownames(cnt_mtx)],
                                                  cellTypeName=tgtCell,
                                                  tgtChrs=c(1:22),
                                                  formula='~ %s + dataset',
                                                  donorID='donorID',groupID='group',geneMap = geneMap[match(rownames(cnt_mtx),geneMap$ensID),])
  dev.off()
}


saveRDS(out,file.path(outDir,'DESeq2_byCT_byMutationGroup.dataset.RDS'))

out = readRDS(file.path(outDir,'mutationGroup/DESeq2_byCT_byMutationGroup.RDS'))

dim(out[['Mono_Mac']][['mainDE']])
deg = do.call(rbind,lapply(out,function(i){return(i[['mainDE']])}))









##------- Using DEG_FindMarkers() v2 ----------####

# hs-SGS (BJ111 + BJ112) vs normal (BJ114 + 3 donors from Wang21)

hs_markers_allCT = data.frame()
for(celltype in c('B.cell','DC2','Mono_CD16','Mono_Mac','NK','pDC','Plasma_cell')){
  print(celltype)
  #setbp1_wang21$group = ifelse(setbp1_wang21$annot == celltype & setbp1_wang21$donorID %in% c('BJ112','BJ114','BJ111'),combSrat$donorID,'others')
  setbp1_wang21$group = ifelse(setbp1_wang21$annot == celltype & setbp1_wang21$donorID %in% c('BJ112','BJ111'),'hs',
                               ifelse(setbp1_wang21$annot == celltype & setbp1_wang21$donorID %in% c('BJ114','25','24','26'),'normal','others'))
  Idents(setbp1_wang21) = setbp1_wang21$group
  hs_markers = FindMarkers(setbp1_wang21,ident.1 = 'hs',ident.2 = 'normal',
                           features = geneMap$geneSym[geneMap$geneSym %in% rownames(setbp1_wang21) & geneMap$chr %in% paste0('chr',c(1:22))],
                           logfc.threshold = 0.1)
  hs_markers$celltype = celltype
  hs_markers = annotateMarkers(hs_markers)
  hs_markers$direction = ifelse(hs_markers$avg_log2FC > 0 ,'hs_up','hs_down')
  
  hs_markers$pct.bj112 = apply(setbp1_wang21@assays$RNA@counts[match(hs_markers$geneSym,rownames(setbp1_wang21)),setbp1_wang21$cellID[setbp1_wang21$group == 'hs' & setbp1_wang21$donorID == 'BJ112']],1,function(x){sum(x>0)/length(x)})
  hs_markers$pct.bj111 = apply(setbp1_wang21@assays$RNA@counts[match(hs_markers$geneSym,rownames(setbp1_wang21)),setbp1_wang21$cellID[setbp1_wang21$group == 'hs' & setbp1_wang21$donorID == 'BJ111']],1,function(x){sum(x>0)/length(x)})
  hs_markers$pct.bj114 = apply(setbp1_wang21@assays$RNA@counts[match(hs_markers$geneSym,rownames(setbp1_wang21)),setbp1_wang21$cellID[setbp1_wang21$group == 'normal' & setbp1_wang21$donorID == 'BJ114']],1,function(x){sum(x>0)/length(x)})
  hs_markers$pct.24 = apply(setbp1_wang21@assays$RNA@counts[match(hs_markers$geneSym,rownames(setbp1_wang21)),setbp1_wang21$cellID[setbp1_wang21$group == 'normal' & setbp1_wang21$donorID == '24']],1,function(x){sum(x>0)/length(x)})
  hs_markers$pct.25 = apply(setbp1_wang21@assays$RNA@counts[match(hs_markers$geneSym,rownames(setbp1_wang21)),setbp1_wang21$cellID[setbp1_wang21$group == 'normal' & setbp1_wang21$donorID == '25']],1,function(x){sum(x>0)/length(x)})
  hs_markers$pct.26 = apply(setbp1_wang21@assays$RNA@counts[match(hs_markers$geneSym,rownames(setbp1_wang21)),setbp1_wang21$cellID[setbp1_wang21$group == 'normal' & setbp1_wang21$donorID == '26']],1,function(x){sum(x>0)/length(x)})
  hs_markers_allCT = rbind(hs_markers_allCT,hs_markers)
}

hs_markers_allCT$pct.normMax = apply(hs_markers_allCT[,c('pct.bj114','pct.24','pct.25','pct.26')],1,max)
hs_markers_allCT$pct.normMin = apply(hs_markers_allCT[,c('pct.bj114','pct.24','pct.25','pct.26')],1,min)

hs_markers_allCT$pct.hsMax = apply(hs_markers_allCT[,c('pct.bj112','pct.bj111')],1,max)
hs_markers_allCT$pct.hsMin = apply(hs_markers_allCT[,c('pct.bj112','pct.bj111')],1,min)


hs_markers_allCT$pct_high = ifelse(hs_markers_allCT$direction == 'hs_up',(hs_markers_allCT$pct.hsMin >= hs_markers_allCT$pct.normMax),
                                   (hs_markers_allCT$pct.hsMax <= hs_markers_allCT$pct.normMin))
hs_markers_allCT$pct_high_diff = ifelse(hs_markers_allCT$direction == 'hs_up',hs_markers_allCT$pct.hsMin - hs_markers_allCT$pct.normMax,
                                        hs_markers_allCT$pct.hsMax - hs_markers_allCT$pct.normMin)

write.csv(hs_markers_allCT,'~/lustre_mt22/SETBP1/Results/3_DEG_FindMarkers/hs.vs.normal_markers_allCT.csv')


##---   Define gene signatures from FM_markers  ---####

hs_markers_allCT = read.csv('~/lustre_mt22/SETBP1/Results/3_DEG_FindMarkers/hs.vs.normal_markers_allCT.csv')
table(hs_markers_allCT$pct_high)
table(hs_markers_allCT$celltype[hs_markers_allCT$pct_high],hs_markers_allCT$direction[hs_markers_allCT$pct_high])
View(hs_markers_allCT[hs_markers_allCT$pct_high == T & abs(hs_markers_allCT$pct_high_diff) >= 0.2,])

genes = hs_markers_allCT$geneSym[hs_markers_allCT$celltype == 'Mono_CD16' & hs_markers_allCT$pct_high == T & abs(hs_markers_allCT$pct_high_diff) >= 0.2]
genes = unique(hs_markers_allCT$geneSym[hs_markers_allCT$direction == 'hs_up' & hs_markers_allCT$pct_high == T & abs(hs_markers_allCT$pct_high_diff) >= 0.2])
length(genes)

Idents(setbp1_wang21) = setbp1_wang21$Phase
DotPlot(setbp1_wang21,idents = c('Mono_CD16'),group.by = 'donorID',
        features = genes)

genes = unique(topGenes$geneSym[topGenes$direction == 'hs_up'])
DotPlot(setbp1_wang21,idents = 'G1',group.by = 'annot',
        features = genes)+
  RotatedAxis()+
  xlab('') + ylab('')+
  theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1),
        axis.text = element_text(size=8))



## Celltype-signature
min_log2FC = 0.25

hs_markers_allCT_allCT = hs_markers_allCT %>% 
  filter(pct_high == T & abs(pct_high_diff) >= 0.1 & abs(avg_log2FC) >= min_log2FC) %>% 
  group_by(geneSym,direction) %>% mutate(nComp = n_distinct(celltype))
table(hs_markers_allCT_allCT$celltype,hs_markers_allCT_allCT$direction)

hs_markers_allCT_setbp1PosOnly = hs_markers_allCT %>% 
  filter(pct_high == T & abs(pct_high_diff) >= 0.1 & abs(avg_log2FC) >= min_log2FC) %>% 
  filter(celltype %in% c('B.cell','Mono_CD16','Mono_Mac','NK','Plasma_cell','pDC')) %>% 
  group_by(geneSym,direction) %>% mutate(nComp = n_distinct(celltype))

topGenes = hs_markers_allCT_allCT[hs_markers_allCT_allCT$nComp >=4,]
topGenes$pct.teratoma = apply(setbp1_wang21@assays$RNA@counts[match(topGenes$geneSym,rownames(setbp1_wang21)),setbp1_wang21$cellID[setbp1_wang21$annot == 'teratoma']],1,function(i){sum(i>0)/length(i)})
topGenes2 = hs_markers_allCT_setbp1PosOnly[hs_markers_allCT_setbp1PosOnly$nComp >= 2,]
table(topGenes2$direction)
n_distinct(topGenes2$geneSym[topGenes2$direction == 'hs_up'])
topGenes = topGenes2





nDEG_perCT = hs_markers_allCT_allCT %>% group_by(celltype) %>% summarise(nGene = n_distinct(geneSym))
nCell_perCT = table(setbp1_wang21$annot)
nDEG_perCT$nCell = nCell_perCT[match(nDEG_perCT$celltype,names(nCell_perCT))]
ggplot(nDEG_perCT,aes(log10(nCell),nGene))+
  geom_point(aes(col=celltype),size=3)+
  scale_color_manual(values = col25)+
  theme_classic()+
  theme(panel.border = element_rect(fill=F),axis.line = element_blank(),
        strip.background=element_rect(linewidth=0),
        strip.text.x = element_text(size=13),
        axis.text.x = element_text(size = 10))



library(UpSetR)
geneList = split(hs_markers_allCT_allCT$geneSym,hs_markers_allCT_allCT$celltype)
geneList = lapply(geneList, function(i){unique(i)})

upset(fromList(geneList),nsets = 20,text.scale = 1.5)




##--------------------------------------------------##
##    Scoring gene signatures in scRNAseq data    ####
##--------------------------------------------------##


hs_markers_allCT_allCT$group = paste0(hs_markers_allCT_allCT$celltype,'.',hs_markers_allCT_allCT$direction)
geneList = split(hs_markers_allCT_allCT$geneSym,hs_markers_allCT_allCT$group)
geneList = lapply(geneList, function(i){unique(i)})
geneList[['topGenes_up']] = unique(topGenes$geneSym[topGenes$direction == 'hs_up'])

setbp1_wang21 <- UCell::AddModuleScore_UCell(setbp1_wang21, features = geneList,ncores = 10)
module_type = 'degFM_v2'
write.csv(setbp1_wang21@meta.data,file.path(outDir,paste0('SETBP1merged_',module_type,'_moduleScore.csv')))
df = read.csv(file.path(outDir,paste0('SETBP1merged_',module_type,'_moduleScore.csv')))
df$mutationGroup = setbp1_wang21$mutationGroup[match(df$cellID,setbp1_wang21$cellID)]
df$mutationGroup[df$dataset == 'wang21'] = 'normal'
df$finalAnn_broad = df$annot
df$Genotype = df$mutationGroup
df$Genotype[df$Genotype =='normal'] = 'Diploid'
df$dataset_ID = ifelse(df$donorID %in% c('L061'),'MDS_AMKL_AEL',
                       ifelse(df$donorID %in% c('BJ111','BJ112','BJ113','BJ9','GOSH084'),'SGS','normalPBMC'))
df$disease = ifelse(df$donorID %in% c('L061'),'MDS',
                    ifelse(df$donorID %in% c('BJ111','BJ112','BJ113','BJ9','GOSH084'),'SGS','normalPBMC'))
df$timePoint='?'
df$module = df$TAM.MLDS.up_UCell
df$dataset_ID_og = df$dataset_ID

data = rbind(data,df[,colnames(data)])

data$finalAnn_broad[data$finalAnn_broad == 'teratoma'] = 'Cancer'
data$finalAnn_broad[data$finalAnn_broad == 'MDS?'] = 'Cancer'
df = data[data$finalAnn_broad %in% c('Cancer','MEMP_MEP','B.cell'),]
#df = df[!(df$dataset_ID %in% c('fLiver') & df$Genotype %in% c('diploid','T21')) & !(df$donorID == 'L041' & df$finalAnn_broad == 'Cancer'),]
df = df[!(df$donorID %in% c('L041','CC3') & df$finalAnn_broad == 'Cancer'),]
df = df[!(df$finalAnn_broad == 'MEMP_MEP' & df$dataset_ID %in% c('pBALL','pAML','MLDS','TAM','MDS','AMKL','AEL','infantALL','LPD')),]

# remove blasts from TP1/TP4, except for L076 and L038
df$toKeep = T
df$toKeep[!grepl('L076|L038|L156',df$donorID) & df$timePoint != 'D' & df$dataset_ID %in% c('MLDS','TAM') & df$finalAnn_broad == 'Cancer'] = F
#df$toKeep[!grepl('L076|L038|L156',df$donorID) & df$timePoint != 'D0' & df$dataset_ID %in% c('pBALL','pAML','MDS','AMKL','AEL','infantALL','LPD') & df$finalAnn_broad == 'Cancer'] = F
df = df[df$toKeep==T,]
#df$toKeep[!df$donorID %in% c('L076','L038') & !is.na(df$timePoint) & df$timePoint != 'Diagnostic' & df$dataset_ID %in% c('MLDS','TAM') & df$finalAnn_broad == 'Tumour'] = F


df$group = ifelse(df$finalAnn_broad != 'Cancer',as.character(df$finalAnn_broad),
                  ifelse(df$dataset_ID %in% c('MLDS','TAM','Recurrent TAM') & df$finalAnn_broad == 'Cancer','TAM / MLDS',
                         ifelse(df$finalAnn_broad == 'Cancer','Other leukaemia','others')))
df$group = factor(df$group,c('B.cell','MEMP_MEP','TAM / MLDS','Other leukaemia'))
df$Genotype[df$Genotype == 'complete_trisomy'] = 'Triploid'
df$Genotype[df$Genotype == 'diploid'] = 'Diploid'
df$Genotype = factor(df$Genotype,c('Diploid','T21','T18','T22','MX','Triploid','adjacent','hotspot'))

## Add AlleleIntegrator call for L038 ##
# aiRes = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x_alleleIntegrator/L038/L038_TumourEry_subClustering_AIRes_240131.csv')
# df$AIres = '?'
# df$AIres[df$cellID %in% aiRes$cellID] = aiRes$AI_output[match(df$cellID[df$cellID %in% aiRes$cellID],aiRes$cellID)]
#df$timePoint = big.srat$timePoint[match(df$cellID,big.srat$cellID)]

# Remove infantALL dataset
df = df[df$dataset_ID != 'infantALL',]

geno_cols = c('Diploid' = grey(0.2),
              'T21' = '#93221E',
              'T18' = '#3d5dad',
              'T22' = '#679551',
              'T13' = '#526691',
              'MX' = '#b18db8',
              'Triploid' = '#e07d26',
              'adjacent'=col25[3],'hotspot'=col25[1])

groupID_order = c(groupID_order,groupID_order_byMedian[grepl('_TAM:',groupID_order_byMedian)],
                  groupID_order_byMedian[grepl('Recurrent TAM:T21:L156$',groupID_order_byMedian)],
                  groupID_order_byMedian[grepl('Recurrent TAM:T21:L156_RD',groupID_order_byMedian)],
                  groupID_order_byMedian[grepl('MLDS_MLDS:',groupID_order_byMedian)],
                  groupID_order_byMedian[grepl('Other.*L061$',groupID_order_byMedian)],
                  groupID_order_byMedian[grepl('Other.*L067$',groupID_order_byMedian)],
                  groupID_order_byMedian[grepl('Other.*L010$',groupID_order_byMedian)],
                  groupID_order_byMedian[grepl('Other.*P9$',groupID_order_byMedian)],
                  groupID_order_byMedian[grepl('Other.*L069$',groupID_order_byMedian)],
                  groupID_order_byMedian[grepl('Other.*_pAML:',groupID_order_byMedian)],
                  groupID_order_byMedian[grepl('Other.*_pBALL:Diploid',groupID_order_byMedian)],
                  groupID_order_byMedian[grepl('Other.*_pBALL:T21',groupID_order_byMedian)],
                  groupID_order_byMedian[grepl('Other.*_LPD:T21',groupID_order_byMedian)],
                  groupID_order_byMedian[grepl('Other.*_SGS',groupID_order_byMedian)],
                  'Other leukaemia_SGS:adjacent:GOSH084','Other leukaemia_SGS:hotspot:BJ111',
                  groupID_order_byMedian[grepl('B.cell.*',groupID_order_byMedian)])

df=df[df$groupID %in% groupID_order,] 








df = setbp1_wang21@meta.data
ggplot(df,aes(annot,topGenes_up_UCell,fill=mutationGroup))+
  geom_boxplot()+
  #facet_wrap(vars(annot))+
  geom_hline(yintercept = 0.2)+theme_classic()+
  theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1),
        axis.text = element_text(size=8))

ggplot(df,aes(annot,TAM.MLDS.up_UCell,fill=mutationGroup))+
  geom_boxplot()+
  #facet_wrap(vars(annot))+
  #geom_hline(yintercept = 0.2)+theme_classic()+
  theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1),
        axis.text = element_text(size=8))


## Score in paed_PBMC
a = readRDS('~/lustre_mt22/SETBP1/Data/scRNAseq_ped_PBMC_ARS2022/pediatricPBMC_SETBP1_integrated.RDS')
set.seed(12234)
paedPBMC <- UCell::AddModuleScore_UCell(paedPBMC, features = geneList,ncores = 10)
write.csv(paedPBMC@meta.data,file.path(outDir,paste0('paedPBMC.ARS22_',module_type,'_moduleScore.csv')))














##----- 2. boxplot of log10CPM count for gene of interests ----##
topGenes
gene_toPlot = topGenes[topGenes$direction == 'hs_up',]
bulkSamples = import_BulkSamples(filter_lowExpr_gene=T)
tpmCnt = bulkSamples[['cpm_count']]
bulk_mdat = bulkSamples[['mdat']]
bulk_mdat$group = ifelse(bulk_mdat$setbp1_mut_inHotspot == T,'hs',
                           ifelse(bulk_mdat$setbp1_mut_inHotspot == F | bulk_mdat$setbp1_mut %in% c('D814N','G816S','P1070L'),'non-hs',bulk_mdat$setbp1_mut))
bulk_mdat$group[bulk_mdat$group == 'FALSE'] = 'no_info'
bulk_mdat$group[grepl('Healthy',bulk_mdat$cancerType)] = 'Healthy'



mtx = tpmCnt[rownames(tpmCnt) %in% gene_toPlot$ensID[gene_toPlot$direction == 'hs_up'],]
mtx = mtx[,!colnames(mtx) %in% c('ensID','geneLength')]
rownames(mtx) = gene_toPlot$geneSym[match(rownames(mtx),gene_toPlot$ensID)]
mtx = mtx[match(unique(gene_toPlot$geneSym[gene_toPlot$geneSym %in% rownames(mtx)]),rownames(mtx)),]
mtx = as.data.frame(t(as.matrix(mtx)))

mtx$sampleID = rownames(mtx)

mtx$group = paste0(bulk_mdat$source[match(mtx$sampleID,bulk_mdat$sampleID)],':',bulk_mdat$group[match(mtx$sampleID,bulk_mdat$sampleID)])


mtx = pivot_longer(mtx,cols = c(1:(which(colnames(mtx) == 'sampleID')-1)),names_to = 'gene',values_to = 'CPM')
mtx = mtx[mtx$group !='TAM',]
mtx = mtx[!is.na(mtx$CPM),]

mtx2 = mtx[mtx$group %in% c('TAM:0','TAM:1'),]
mtx2 = mtx2 %>% group_by(group,gene) %>% summarise(med = median(CPM))
mtx2 = pivot_wider(mtx2,id_cols = 'gene',names_from = 'group',values_from = 'med')
mtx2$med_diff = mtx2$`TAM:1` - mtx2$`TAM:0`



library(ggbeeswarm)
mtx$gene = factor(mtx$gene,unique(gene_toPlot$geneSym))
ggplot(mtx[mtx$gene %in% unique(gene_toPlot$geneSym)[81:100],],aes(group,log(CPM)))+
  geom_boxplot(aes(fill=group),outlier.colour = 'white',position = 'dodge', alpha = 1,width=0.5,linewidth=0.3)+
  geom_point(size=0.5,alpha=0.7,col=grey(0))+
  scale_fill_manual(values = c(col25[2],col25[3],col25[1],col25[2],col25[3],grey(0.7),col25[3],grey(0.7),col25[3])) +
  facet_wrap(vars(gene),nrow=3,scales = 'free_y')+
  theme_classic()+
  #ggtitle(title)+xlab('')+ylab('Module score')+
  theme(panel.border = element_rect(fill=F),axis.line = element_blank(),
        strip.background=element_rect(linewidth=0),
        axis.text.x = element_text(size = 10,angle = 90, vjust = 0.5,hjust = 1))


mtx = mtx[mtx$gene %in% c('CD83','CCNL1','ARIH1','STK17B','MBD2','PTEN','ARID4B','DDX6','BEST1','TACC1','SMCHD1','MAP3K2','MAP4K4','ROCK1','AGFG1','CHD1','CLK1','USP3','SOS2','CDV3','KMT2C','ITCH','MEX3C'),]
mtx$gene = factor(mtx$gene,c('CDKN2C', 'CA1','CD81','PTK2','TMX4','MSRB3','ZHX1','ZNF141','PRAME','CD44','CD74','CPA3','ARHGEF12','MAGI1','PGM1','RRAS2','DSTN','LNCAROD','SLC5A3','BST2','ZNF345','LINS1'))
ggplot(mtx,aes(group,log(CPM)))+
  geom_boxplot(aes(fill=group),outlier.colour = 'white',position = 'dodge', alpha = 1,width=0.5,linewidth=0.3)+
  geom_point(size=0.5,alpha=0.7,col=grey(0))+
  scale_fill_manual(values = c(col25[2],col25[3],col25[1],col25[2],col25[3],grey(0.7),col25[3],grey(0.7),col25[3])) +
  #scale_fill_manual(values = c('white','#c7065a',col25[4],rep(colAlpha(col25[1],0.4),2),rep(grey(0.7),5))) +
  facet_wrap(vars(gene),nrow=3,scales = 'free_y')+
  theme_classic()+
  xlab('')+
  #ggtitle(title)+xlab('')+ylab('Module score')+
  theme(panel.border = element_rect(fill=F),axis.line = element_blank(),
        strip.background=element_rect(linewidth=0),
        strip.text.x = element_text(size=13),
        axis.text.x = element_text(size = 10,angle = 90, vjust = 0.5,hjust = 1))+
  ggtitle('Bulk RNA-seq dataset')


setbp1_wang21$group_4dotplot = ifelse(setbp1_wang21$annot %in% c('MDS?','teratoma'),setbp1_wang21$annot,
                                      setbp1_wang21$donorID)
DotPlot(setbp1_wang21,group.by = 'group_4dotplot',
        features = c('CD83','CCNL1','ARIH1','STK17B','MBD2','PTEN','ARID4B','DDX6','BEST1','TACC1','SMCHD1','MAP3K2','MAP4K4','ROCK1','AGFG1','CHD1','CLK1','USP3','SOS2','CDV3','KMT2C','ITCH','MEX3C'))+
  RotatedAxis()+
  xlab('') + ylab('')+
  theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1),
        axis.text = element_text(size=8))







