## Generates the Figures for paper ##

#################
# Materials map #
#################

## Main Figures ##






## Supplementary Figures ##



##----    Set working directory  -----##
setwd('~/lustre_mt22/SETBP1/')


#------------------------#
##      Libraries     ####
#------------------------#
library(tidyverse)
library(Seurat)
library(readxl)
source('/lustre/scratch126/casm/team274sb/mt22/CN_methods/scripts/finalScripts/R/misc.R')
source('~/lustre_mt22/generalScripts/utils/sc_utils.R')
source('~/lustre_mt22/generalScripts/utils/misc.R')

plotDir='~/lustre_mt22/SETBP1/manuscriptDraft_0325/Plots'
if(!dir.exists(plotDir)){
  dir.create(plotDir,recursive = T)
}


##----------------------------------##
##    Set general color schemes   ####
##----------------------------------##



##-----------------------------------##
##      EGA requests for SETBP1    ####
##-----------------------------------##
projectMani = read_excel('~/lustre_mt22/projectManifest.xlsx','SETBP1')

# sc/sn-RNAseq
setbp1_dataset = projectMani[!is.na(projectMani$assay) & !grepl('^L\\d+|HCA_',projectMani$DonorID) &
                               grepl('scRNA|snRNA|TCR',projectMani$assay),]
setbp1_dataset.sub = setbp1_dataset[,c('DonorID','externalID','PDID','Tissue','Sex','mutation','assay','sangerSampleID','SequenceScape_study','Supplier_ID')]


write.csv(setbp1_dataset.sub,'~/lustre_mt22/SETBP1/scRNAseq.data_info_for_EGA.csv')


# WGS
setbp1_wgs_dataset = projectMani[!is.na(projectMani$DonorID) & !is.na(projectMani$assay) & 
                                   projectMani$assay == 'WGS'  & projectMani$DonorID != 'L046',]

setbp1_wgs_dataset$Condition = ifelse(setbp1_wgs_dataset$DonorID %in% c('BJ111','BJ112','BJ9'),'SGS',
                                      ifelse(setbp1_wgs_dataset$DonorID %in% c('BJ113','GOSH084'),'aSGS',
                                             ifelse(setbp1_wgs_dataset$DonorID == 'L061','diploid, MDS with SETBP1 mutation',
                                                    ifelse(setbp1_wgs_dataset$DonorID == 'L067','diploid, MDS with GATA1 mutation','Healthy'))))
setbp1_wgs_dataset$mutation[setbp1_wgs_dataset$DonorID == 'L067'] = 'GATA1 mutation'

setbp1_dataset.sub = setbp1_wgs_dataset[,c('DonorID','externalID','Condition','Tissue','Sex','mutation','assay','PDID','canapp_project','SequenceScape_study','Supplier_ID','sangerSampleID')]
colnames(setbp1_dataset.sub) = c('DonorID','SampleID','Condition','Tissue','Sex','Mutation','Assay','PDID','Canapp_project','SequenceScape_study','Supplier_ID','SangerSampleID')
table(setbp1_dataset.sub$DonorID,setbp1_dataset.sub$Condition)
table(setbp1_dataset.sub$Mutation,setbp1_dataset.sub$Condition)

setbp1_dataset.sub = setbp1_dataset.sub[setbp1_dataset.sub$DonorID != 'L067',]
write.csv(setbp1_dataset.sub,'~/lustre_mt22/SETBP1/WGS.data_info_for_EGA.csv')



##---------------------------##
##    Dataset overview     ####
##---------------------------##
setbp1_srat_fp = '~/lustre_mt22/SETBP1/Results/x_published_scRNAseq_preprocessing/SETBP1_Wang21_combined_2410.RDS'
setbp1_mdat_fp = '~/lustre_mt22/SETBP1/Results/2_annotation/SETBP1_merged_processed_annot_tmp_oct24.csv'




fig1a_datasetSummary = function(){
  ## Import mdat dataset ##
  mdat = read.csv(setbp1_mdat_fp)
  # Remove doublets
  mdat = mdat[mdat$annot != 'doublets',]
  mdat$condition = ifelse(!mdat$condition %in% c('MDS_SETBP1','unaffected','MDS_GATA1'),paste0('SGS (',mdat$condition,')'),mdat$condition)
  mdat$condition[mdat$condition == 'MDS_SETBP1'] = 'MDS_G870S'
  
  mdat$disease2 = factor(mdat$disease2,c('SGS','MDS','healthy'))
  mdat$source = ifelse(mdat$dataset == 'wang21','GSE168732','this_study')
  
  ## Write table of number of cells per donorID
  dataset = as.data.frame(table(mdat$source,mdat$donorID,mdat$setbp1_mutation_status,mdat$condition,mdat$tissue,mdat$sex,mdat$assay,mdat$orig.ident))
  colnames(dataset) = c('Source','Donor_ID','SETBP1_mutation_status','Condition','Tissue','Sex','Assay','Sample_ID','nCell')
  dataset$Donor_ID = factor(dataset$Donor_ID,c('GOSH084','BJ113',
                                               'BJ111','BJ112','BJ9','L061','L067',
                                               'BJ114','GSM5160432','GSM5160434','GSM5160435'))
  dataset = dataset[order(dataset$Donor_ID),]
  
  # Add PDID
  projectMani = read_excel('~/lustre_mt22/projectManifest.xlsx',sheet = 'SETBP1')
  projectMani = projectMani[projectMani$assay == 'WGS' & !is.na(projectMani$DonorID) & !projectMani$PDID %in% c('PD54850b','PD54851b','PD62642a') & !grepl('Saliva|Cells|fibroblast|^BM$',projectMani$Tissue),]
  projectMani$DonorID[projectMani$DonorID == 'USA_case'] = 'BJ9'
  
  dataset$PDID = projectMani$PDID[match(dataset$Donor_ID,projectMani$DonorID)]
  # dd$cellRanger_version = ifelse(dd$assay == "scRNA 3' v3.1",'cellranger_302','cellranger_700')
  # dd$refGenome_version = ifelse(dd$assay == "scRNA 3' v3.1",'GRCh38 1.2.0','GRCh38 2020-A')
  
  write.csv(dataset[dataset$nCell >0,],file.path(plotDir,'..',paste0('TableS1_dataset_nCellBy10XSample.csv')),row.names = F)
  
  
  

  ## Make heatmap of samples
  mdat$donorID2 = paste0(mdat$donorID, ' ',mdat$condition)
  dataset = as.data.frame(table(mdat$donorID2,mdat$tissue))
  colnames(dataset) = c('donorID','Tissue','nCell')
  
  dataset$sampleAvailable = ifelse(dataset$nCell > 0,1,0)
  dataset = pivot_wider(dataset,id_cols = 'Tissue',names_from = 'donorID',values_from = 'sampleAvailable') 
  dataset = column_to_rownames(dataset,'Tissue')
  dataset = (as.matrix(dataset))
  
  library(ComplexHeatmap)
  
  plotFun_SETBP1_dataset = function(noFrame=FALSE,noPlot=FALSE){
    par(mar=c(0.1,0.1,1,0.1))
    botAnno = HeatmapAnnotation(df = data.frame('setbp1_mut' = mdat$setbp1_mutation_status[match(colnames(dataset),mdat$donorID2)],
                                           'germline_vs_somatic' = ifelse(grepl('L061',colnames(dataset)),'somatic',
                                                                          ifelse(grepl('SGS',colnames(dataset)),'germline','unaffected'))),
                                col = list(setbp1_mut = c('wild_type' = 'white','hotspot' = '#cc5d0f','non_hotspot' = '#1b4891'),
                                           germline_vs_somatic = c('somatic' = pal34H[14],'germline'=pal34H[2],'unaffected'=grey(1))))
                                #annotation_legend_param = list(Genotype = list(direction = "horizontal")))
    
    botAnno = HeatmapAnnotation(df = data.frame('setbp1_mut' = mdat$setbp1_mutation_status[match(colnames(dataset),mdat$donorID2)]),
                                col = list(setbp1_mut = c('wild_type' = 'white','hotspot' = '#cc5d0f','non_hotspot' = '#1b4891')))
    
    
    
    hm = Heatmap(as.matrix(dataset),name = 'sample',col = c("1"=grey(0.2),"0"=grey(1)),
                 cluster_rows = F,cluster_columns = F,row_names_side = "left",show_column_dend = F,
                 column_split = mdat$disease2[match(colnames(dataset),mdat$donorID2)],
                 border = T,rect_gp = gpar(col = "black", lwd = 1.2),bottom_annotation = botAnno)
    draw(hm)
  }
  
  saveFig(file.path(plotDir,paste0('Fig1_SETBP1_dataset')),plotFun_SETBP1_dataset,rawData=mdat,width = 5.53,height = 3.15,res = 500,useDingbats = F)
  
}





srat = readRDS(setbp1_srat_fp)
DimPlot(srat,group.by = 'condition',cols = col25)
sc_srat = subset(srat,subset = cellID %in% srat2$cellID[srat2$assay == "scRNA 10x5'v2(DUAL)"])

## Import L067
otherLeuk_srat_fp = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/otherLeuk/otherLeuk_clean_annotated_noUnknown_2408.RDS'
otherLeuk = readRDS(otherLeuk_srat_fp)
L067 = subset(otherLeuk,subset = donorID == 'L067')
rm(otherLeuk)

sc_srat = merge_seurat_objects(sc_srat,L067,keepAllGenes = F,genomeVersions = c('v38','v38'))
sc_srat = standard_clustering(sc_srat)
saveRDS(sc_srat,'~/lustre_mt22/SETBP1/Results/x_published_scRNAseq_preprocessing/SETBP1_Wang21_scRNAseq.only_2502.RDS')

## Harmony to regress out patient specific effect
sc_srat2 = standard_clustering(sc_srat,runHarmony = T,harmonyVar = 'donorID')

## Refine annotation
sc_srat$finalAnn_broad = as.character(sc_srat$annot)
sc_srat$finalAnn_broad[sc_srat$finalAnn_broad %in% c('T_CD4','T_CD8','T_cell','T_cells','T_gd','T_MAIT')] = 'T_cells'
sc_srat$finalAnn_broad[sc_srat$finalAnn_broad %in% c('B.cell','naive.B')] = 'B_cells'
sc_srat$finalAnn_broad[sc_srat$finalAnn_broad %in% c('EE','LE')] = 'Ery'
sc_srat$finalAnn_broad[sc_srat$finalAnn_broad %in% c('Plasma_cell','Plasma.cell')] = 'Plasma.cell'
sc_srat$finalAnn_broad[sc_srat$finalAnn_broad == 'MDS?' & sc_srat$donorID != 'L061'] = 'doublets'

DimPlot(sc_srat,group.by = 'finalAnn_broad',cols = c(col25,pal34H),label = T,repel = T,label.box = T)
DimPlot(sc_srat2,group.by = 'seurat_clusters',cols = c(col25,pal34H),label = T,label.box = T,repel = T) + NoLegend()
DimPlot(sc_srat,cells.highlight = sc_srat$cellID[sc_srat$finalAnn_broad == '41'])
DimPlot(sc_srat,cells.highlight = sc_srat2$cellID[sc_srat2$seurat_clusters == 39])
DimPlot(sc_srat,cells.highlight = sc_srat$cellID[sc_srat$finalAnn_broad == 'MDS?' & sc_srat$donorID != 'L061'])

sc_srat$group = ifelse(sc_srat$seurat_clusters %in% c(10,12,3,40,29),'g1',ifelse(sc_srat$seurat_clusters == 11,'g2','others'))

library(SoupX)
qm = quickMarkers(sc_srat@assays$RNA@counts[,sc_srat$cellID[sc_srat$group %in% c('g1','g2')]],sc_srat$group[sc_srat$group %in% c('g1','g2')])


sc_srat2$group = ifelse(sc_srat2$seurat_clusters %in% c(0,11,10,19,12,3,5,15),'g1',ifelse(sc_srat2$seurat_clusters %in% c(7,21,6),'g2','others'))
sc_srat2$group = ifelse(sc_srat2$seurat_clusters %in% c(24) & sc_srat2$finalAnn_broad == 'MDS?','g1',
                        ifelse(sc_srat2$seurat_clusters %in% c(24) & sc_srat2$finalAnn_broad != 'MDS?','g3',
                               ifelse(sc_srat2$finalAnn_broad2 == 'T_cells','g4',
                                      ifelse(sc_srat2$seurat_clusters %in% c(33),'g2','others'))))
sc_srat2$group = ifelse(sc_srat2$seurat_clusters == 17 & sc_srat2$finalAnn_broad %in% c('41','MDS?','Tumour','T_cells','NK'),sc_srat2$finalAnn_broad,
                        'others')
qm = quickMarkers(sc_srat2@assays$RNA@counts[,sc_srat2$cellID[sc_srat2$group %in% c('g1','g3')]],sc_srat2$group[sc_srat2$group %in% c('g1','g3')])
qm = quickMarkers(sc_srat2@assays$RNA@counts[,sc_srat2$cellID[!sc_srat2$group %in% c('others')]],sc_srat2$group[!sc_srat2$group %in% c('others')])
qm = quickMarkers(sc_srat2@assays$RNA@counts[,sc_srat2$cellID[!sc_srat2$finalAnn_broad2 %in% c('doublets','50','?')]],sc_srat2$finalAnn_broad2[!sc_srat2$finalAnn_broad2 %in% c('doublets','50','?')])
DimPlot(sc_srat2,cells.highlight = sc_srat$cellID[sc_srat$cellID %in% sc_srat2$cellID[sc_srat2$seurat_clusters %in% c(21,44,6,38,37,47)] & sc_srat$finalAnn_broad == '?'])
table(sc_srat2$donorID[sc_srat2$seurat_clusters == 24],sc_srat2$finalAnn_broad[sc_srat2$seurat_clusters == 24])
DimPlot(sc_srat2,cells.highlight = sc_srat2$cellID[sc_srat2$seurat_clusters == 24 & sc_srat2$finalAnn_broad == 'MDS?'])
DimPlot(sc_srat,cells.highlight = sc_srat2$cellID[sc_srat2$finalAnn_broad == 'B_cells' & !sc_srat2$seurat_clusters %in% c(30,9,34,2,18)])


## Fix annotation
sc_srat$finalAnn_broad[sc_srat$cellID %in% sc_srat2$cellID[sc_srat2$seurat_clusters %in% c(21,44,6,38,37,47)] & sc_srat$finalAnn_broad %in% c('T_cells','?')]  = 'NK_T'

## Sub-cluster NK / T cells
t.cells = subset(sc_srat,subset = finalAnn_broad %in% c('NK','NK_T','T_cells'))
t.cells = standard_clustering(t.cells)
sc_srat$finalAnn_broad[sc_srat$cellID %in% t.cells$cellID[t.cells$seurat_clusters %in% c(25,11,33,26)]] = 'NK_T'
t.cells$finalAnn_broad = sc_srat$finalAnn_broad[match(t.cells$cellID,sc_srat$cellID)]
sc_srat2$finalAnn_broad = sc_srat$finalAnn_broad[match(sc_srat2$cellID,sc_srat$cellID)]
DimPlot(t.cells,group.by = 'finalAnn_broad',cols = col25,label = T,repel = T,label.box = T)
DimPlot(t.cells,group.by = 'seurat_clusters',label = T,repel = T,label.box = T)+NoLegend()
DimPlot(t.cells,cells.highlight = t.cells$cellID[t.cells$seurat_clusters == 26])
DimPlot(t.cells,cells.highlight = sc_srat2$cellID[sc_srat2$seurat_clusters == 21 & sc_srat2$finalAnn_broad=='T_cells'])
FeaturePlot(t.cells,c('NCAM1','CD3D','KLRD1','NKG7','CD8A','CD4'))

## Sub-cluster MonoMac cells
mono.mac = subset(sc_srat,subset = finalAnn_broad %in% c('DC1','DC2','pDC','Mono_CD16','Mono_Mac','strange_MonoMac'))
mono.mac = standard_clustering(mono.mac)
DimPlot(mono.mac,group.by = 'finalAnn_broad',cols = col25,label = T,repel = T,label.box = T)
DimPlot(mono.mac,group.by = 'seurat_clusters',cols = c(col25,pal34H),label = T,repel = T,label.box = T)
DimPlot(sc_srat2,cells.highlight = mono.mac$cellID[mono.mac$seurat_clusters == 13])
FeaturePlot(mono.mac,c('CD14','FCGR3A','CLEC4C','JCHAIN'))
qm_mono = quickMarkers(mono.mac@assays$RNA@counts,mono.mac$finalAnn_broad)


DotPlot(sc_srat,#idents = unique(combSrat$seurat_clusters[is.na(combSrat$cluster_ann)]),
        group.by = 'finalAnn_broad',features = unique(c('IGFBP1','IGFBP2','HMGA2','MEG3','RBPMS',
                                                        'CD34','CD38','SPINK2','MLLT3','PRSS57', # HSC_MPP
                                                         'SERPINB1', 'GATA1',	'GATA2', 'TESPA1',	'CTNNBL1',#MEMP
                                                         'FCER1A', 'ITGA2B', 'HBD','KLF1','GYPA','PLEK', # MEP
                                                         'ZBTB16','LTB', 'CD52',# ILC precursor
                                                         'IL7R','CD3D','CD3E','GZMA','CD4','CD8A', #Early.lymphoid_T.lymphocyte
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
                                                         'CLEC9A',#DC1
                                                         'CLEC10A','CD1C', # DC2 
                                                         'IL3RA', # DC2
                                                         'CLEC4C', #pDC
                                                         
                                                         'CD27' # Plasma cells
                                                         
                                                         
        )))+RotatedAxis() + ggtitle('Celltype marker genes') + xlab('')+
  theme(axis.text.x = element_text(size=8,angle=90,vjust = 0.5,hjust = 1))



l061 = subset(sc_srat,subset = donorID == 'L061')
l061 = standard_clustering(l061)
DimPlot(l061,group.by = 'Phase',cols = col25,label = T,repel = T)
DimPlot(l061,cells.highlight = l061$cellID[l061$finalAnn_broad == 'MDS?'])
table(l061$finalAnn_broad)


## Finalise the annotation
sc_srat$finalAnn_broad[sc_srat$finalAnn_broad %in% c('MDS?','Tumour')] = 'MDS'
sc_srat$broadLineage = as.character(sc_srat$finalAnn_broad)
sc_srat$broadLineage[sc_srat$broadLineage %in% c('T_cells','NK_T','NK')] = 'NK.T'
sc_srat$broadLineage[sc_srat$broadLineage %in% c('DC1','DC2','Mono_CD16','Mono_Mac','pDC','strange_MonoMac')] = 'Myeloid'
sc_srat$broadLineage[sc_srat$broadLineage %in% c('Ery','MK')] = 'Ery.MK'
sc_srat$broadLineage[sc_srat$broadLineage %in% c('?','50','41')] = 'unknown'
sc_srat$broadLineage[sc_srat$broadLineage %in% c('Plasma.cell')] = 'B_cells'

table(sc_srat$broadLineage,sc_srat$finalAnn_broad)
DimPlot(sc_srat,group.by = 'disease',cols = c(col25,pal34H),label = T,repel = T,label.box = T)

mdat = cbind(sc_srat@meta.data,sc_srat@reductions$umap@cell.embeddings)
write.csv(mdat,'~/lustre_mt22/SETBP1/Results/x_published_scRNAseq_preprocessing/SETBP1_Wang21_combined_annot_2410.csv')



sc_srat = subset(sc_srat,subset = finalAnn_broad %in% unique(sc_srat$finalAnn_broad[!sc_srat$finalAnn_broad %in% c('?','41','50','doublets')]))
sc_srat = standard_clustering(sc_srat)

sc_srat$condition[sc_srat$donorID == 'L067'] = 'MDS_GATA1'
sc_srat$condition[is.na(sc_srat$condition)] = 'unaffected'
sc_srat$disease = ifelse(sc_srat$donorID %in% c("GOSH084",'BJ113'),'aSGS',
                      ifelse(sc_srat$donorID %in% c("BJ9",'BJ111','BJ112'),'SGS',
                             ifelse(sc_srat$donorID %in% c('L061','L067'),'MDS',
                                    ifelse(sc_srat$donorID == 'BJ114','adult','children'))))

sc_srat$assay[sc_srat$assay == "scRNA 10x5'v2(DUAL)"] = "scRNA 10x5'v2"
sc_srat$setbp1_mutation = sc_srat$condition
sc_srat$setbp1_mutation[sc_srat$setbp1_mutation == 'MDS_SETBP1'] = 'G870S'
sc_srat$setbp1_mutation_status = ifelse(sc_srat$donorID %in% c('BJ9','BJ111','BJ112','L061'),'hotspot',
                                     ifelse(sc_srat$donorID %in% c('GOSH084','BJ113'),'non_hotspot','wild_type'))
sc_srat$disease2 = ifelse(sc_srat$disease %in% c('adult','children'),'healthy',
                       ifelse(sc_srat$disease %in% c('aSGS'),'SGS',sc_srat$disease))
sc_srat$disease2 = factor(sc_srat$disease2,c('SGS','MDS','healthy'))

sc_srat$tissue[sc_srat$dataset == 'wang21'] = 'PBMC'

saveRDS(sc_srat,'~/lustre_mt22/SETBP1/Results/x_published_scRNAseq_preprocessing/SETBP1_Wang21_scRNAseq.only_noUnknowns_2502.RDS')


mdat = cbind(srat@meta.data,srat@reductions$umap@cell.embeddings)
write.csv(mdat,'~/lustre_mt22/SETBP1/Results/x_published_scRNAseq_preprocessing/SETBP1_Wang21_scRNAseq.only_noUnknowns_mdat_2502.csv')


setbp1_mdat_fp = '~/lustre_mt22/SETBP1/Results/x_published_scRNAseq_preprocessing/SETBP1_Wang21_scRNAseq.only_noUnknowns_mdat_2502.csv'
fig1b_SETBP1.MDS_UMAP = function(){
  
  if(file.exists(file.path(plotDir,'Fig1b_SETBP1.MDS_UMAP_rawData.tsv'))){
    dd = read.delim(file.path(plotDir,'Fig1b_SETBP1.MDS_UMAP_rawData.tsv'),header = T,sep = '\t')
    
  }else{
    mdat = read.csv(setbp1_mdat_fp)
    mdat$disease[mdat$donorID == 'BJ113'] = 'aSGS'
    
    dd = mdat[,c("cellID","donorID","broadLineage",'finalAnn_broad','disease','UMAP_1','UMAP_2')]
    dd = dd[dd$finalAnn_broad != 'doublets',]
    dd = dd[sample(1:nrow(dd),nrow(dd)),]
  }
  
  disease_col = c('adult' = colAlpha(grey(0.6),0.3),
                  'children' = colAlpha(grey(0.75),0.3),
                  'SGS' = colAlpha('#cc5d0f',0.55),
                  'aSGS' = colAlpha('#1b4891',0.55),
                  'MDS' = colAlpha('#a62121',0.3))
  
  SETBP1_mutation_col = c('wild_type' = colAlpha(grey(0.6),0.3),
                          'hotspot' = colAlpha('#cc5d0f',0.55),
                          'non_hotspot' = colAlpha('#1b4891',0.55))
  dd$setbp1_mutation_status = ifelse(dd$donorID %in% c('BJ9','BJ111','BJ112','L061'),'hotspot',
                                     ifelse(dd$donorID %in% c('GOSH084','BJ113'),'non_hotspot','wild_type'))
  
  plotFun_setbp1_mutation = function(noFrame=FALSE,noPlot=FALSE){
    par(mar=c(0,0,0.1,0))
    
    plot(dd$UMAP_1,dd$UMAP_2,
         las=1,
         type='n',
         cex.main = 0.85,xaxt='n',yaxt='n',
         xlab='',ylab='',
         frame.plot=F)
    
    if(!noPlot){
      points(dd$UMAP_1,dd$UMAP_2,
             col = SETBP1_mutation_col[as.character(dd$setbp1_mutation_status)],
             pch = 19,
             cex=0.07)
    }
  }
  
  saveFig(file.path(plotDir,'Fig1b_SETBP1.MDS_setbp1.mutation_UMAP'),plotFun_setbp1_mutation,rawData=mdat,width = 4.1,height = 3.7,res = 500,useDingbats = F)
  
  
  broadLin_col=c('B_cells' = pal34H[34],
                 'Ery.MK' = col25[8],
                 'MDS' = 'darkred',
                 'Myeloid' = pal34H[9],
                 'NK.T' = pal34H[14]
                 )
  plotFun_broadLin = function(noFrame=FALSE,noPlot=FALSE){
    par(mar=c(0,0,0.1,0))
    
    plot(dd$UMAP_1,dd$UMAP_2,
         las=1,
         type='n',
         cex.main = 0.85,xaxt='n',yaxt='n',
         xlab='',ylab='',
         frame.plot=F)
    
    if(!noPlot){
      points(dd$UMAP_1,dd$UMAP_2,
             col = broadLin_col[as.character(dd$broadLineage)],
             pch = 19,
             cex=0.07)
    }
  }
  
  saveFig(file.path(plotDir,'Fig1b_SETBP1.MDS_broadLin_UMAP'),plotFun_broadLin,rawData=mdat,width = 4.1,height = 3.7,res = 500,useDingbats = F)
  
  
  
  ## Lineage proportion of cells
  dd = mdat %>% group_by(donorID,condition,disease,setbp1_mutation_status,broadLineage) %>% summarise(nCell = n()) %>% 
    group_by(donorID,condition,disease,setbp1_mutation_status) %>% mutate(total_nCell = sum(nCell))
  dd$frac = dd$nCell/dd$total_nCell
  
  SETBP1_mutation_col = c('wild_type' = colAlpha(grey(0.6),1),
                          'hotspot' = colAlpha('#cc5d0f',1),
                          'non_hotspot' = colAlpha('#1b4891',1))
  
  dd$disease[dd$disease == 'aSGS'] = 'atypical SGS'
  dd$disease[dd$disease == 'adult'] = 'normal (PD66164)'
  dd$disease[dd$disease == 'children'] = 'normal (Wang et al., 2021)'
  dd$disease = factor(dd$disease,c('MDS','SGS','atypical SGS','normal (PD66164)','normal (Wang et al., 2021)'))
  dd = dd[dd$disease != 'MDS',]
  
  dd$broadLineage = factor(dd$broadLineage,c('Ery.MK','Myeloid','B_cells','NK.T'))
  
  plotFun_broadLin_proportion = function(noFrame=FALSE,noPlot=FALSE){
    p = ggplot(dd,aes(disease,frac,fill=setbp1_mutation_status))+
      geom_boxplot(outlier.shape = NA,alpha=0.9)+
      geom_point()+
      scale_fill_manual(values = SETBP1_mutation_col)+
      facet_grid(~broadLineage,scales = 'free_x',space = 'free_x')+
      theme_classic(base_size = 11)+xlab('')+
      ylab('Fraction of cells')+
      theme(panel.border = element_rect(fill=F,colour = 'black'),axis.line = element_blank(),
            strip.background=element_blank(),
            strip.text.x = element_text(size=10,colour = 'black'),
            axis.ticks = element_line(colour = 'black',linewidth = 0.2),
            axis.text.x = element_text(size = 9,angle = 90, vjust = 0.5,hjust = 1,colour = 'black'),
            axis.text.y = element_text(size=10),
            axis.text = element_text(colour = 'black'),
            legend.text = element_text(colour = 'black',size = 10),
            legend.title = element_text(colour = 'black',size = 12))
    
    print(p)
  }
  
  saveFig(file.path(plotDir,'SupFig1d_lineageFraction'),plotFun_broadLin_proportion,rawData=dd,width = 5.8,height = 4.2,res = 500,useDingbats = F)
  
  
  wilcox.test(dd$frac[dd$broadLineage == 'NK.T' & dd$disease %in% c('SGS','atypical SGS')],
              dd$frac[dd$broadLineage == 'NK.T' & dd$disease %in% c('normal (PD66164)','normal (Wang et al., 2021)')],alternative = 'less')
  wilcox.test(dd$frac[dd$broadLineage == 'Myeloid' & dd$disease %in% c('SGS','atypical SGS')],
              dd$frac[dd$broadLineage == 'Myeloid' & dd$disease %in% c('normal (PD66164)','normal (Wang et al., 2021)')],alternative = 'greater')
  
  prop.test(x=c(dd$nCell[dd$broadLineage == 'Myeloid' & dd$disease %in% c('SGS','atypical SGS')],
                dd$nCell[dd$broadLineage == 'Myeloid' & dd$disease %in% c('normal (PD66164)','normal (Wang et al., 2021)')]),
            n=c(dd$total_nCell[dd$broadLineage == 'Myeloid' & dd$disease %in% c('SGS','atypical SGS')],
                dd$total_nCell[dd$broadLineage == 'Myeloid' & dd$disease %in% c('normal (PD66164)','normal (Wang et al., 2021)')]), 
            p = NULL, alternative = "greater",
            correct = FALSE)
  
  prop.test(x=c(dd$nCell[dd$broadLineage == 'Myeloid' & dd$disease %in% c('SGS','atypical SGS')],
                dd$nCell[dd$broadLineage == 'Myeloid' & dd$disease %in% c('normal (PD66164)','normal (Wang et al., 2021)')]),
            n=c(dd$total_nCell[dd$broadLineage == 'Myeloid' & dd$disease %in% c('SGS','atypical SGS')],
                dd$total_nCell[dd$broadLineage == 'Myeloid' & dd$disease %in% c('normal (PD66164)','normal (Wang et al., 2021)')]), 
            p = NULL, alternative = "greater",
            correct = FALSE)
  
  
  
  
  
  
  
  
  donorID_col=c('PD60299' = pal34H[9],
                 'PD66174' = pal34H[27],
                # SGS
                 'PD60412' = pal34H[17],
                 'PD66162' = pal34H[18],
                 'PD66163' = brewer.pal(n=10,'Paired')[7],
                # Normal
                'PD66164' = grey(0.4),
                'GSM5160432' = grey(0.85),
                'GSM5160434' = grey(0.7),
                'GSM5160435' = grey(0.55),
                # MDS
                'PD61858' = pal34H[3],
                'PD61857' = pal37H[14]
  )
  
  #show_col(pal37H)
  
  dd$donorID[dd$donorID == 'GOSH084'] = 'PD60299'
  dd$donorID[dd$donorID == 'BJ113'] = 'PD66174'
  dd$donorID[dd$donorID == 'BJ9'] = 'PD60412'
  dd$donorID[dd$donorID == 'BJ111'] = 'PD66162'
  dd$donorID[dd$donorID == 'BJ112'] = 'PD66163'
  dd$donorID[dd$donorID == 'BJ114'] = 'PD66164'
  dd$donorID[dd$donorID == 'L061'] = 'PD61858'
  dd$donorID[dd$donorID == 'L067'] = 'PD61857'
  
  plotFun_donorID = function(noFrame=FALSE,noPlot=FALSE){
    par(mar=c(0,0,0.1,0))
    
    plot(dd$UMAP_1,dd$UMAP_2,
         las=1,
         type='n',
         cex.main = 0.85,xaxt='n',yaxt='n',
         xlab='',ylab='',
         frame.plot=F)
    
    if(!noPlot){
      points(dd$UMAP_1,dd$UMAP_2,
             col = donorID_col[as.character(dd$donorID)],
             pch = 19,
             cex=0.07)
    }
  }
  
  saveFig(file.path(plotDir,'SupFig1b_SETBP1.MDS_donorID_UMAP'),plotFun_donorID,rawData=mdat,width = 4.1,height = 3.7,res = 500,useDingbats = F)
  
  
  
  
  
}


mds = subset(srat,subset = donorID %in% c('L061','L067'))
mds = standard_clustering(mds)
saveRDS(mds,'~/lustre_mt22/SETBP1/Results/MDS_L061.L067_HARM_2502.RDS')
saveRDS(mds,'~/lustre_mt22/SETBP1/Results/MDS_L061.L067_2502.RDS')

DimPlot(mds,group.by = 'donorID',cols = col25)
DimPlot(mds,cells.highlight = mds$cellID[mds$finalAnn_broad == 'MDS' & mds$donorID == 'L067'])

mdat = cbind(mds@meta.data,mds@reductions$umap@cell.embeddings)
write.csv(mdat,'~/lustre_mt22/SETBP1/Results/MDS_L061.L067_mdat_2502.csv')

mds_mdat_fp = '~/lustre_mt22/SETBP1/Results/MDS_L061.L067_mdat_2502.csv'

fig1_MDS_UMAP = function(){
  
  if(file.exists(file.path(plotDir,'Fig1b_SETBP1.MDS_UMAP_rawData.tsv'))){
    dd = read.delim(file.path(plotDir,'Fig1b_SETBP1.MDS_UMAP_rawData.tsv'),header = T,sep = '\t')
    
  }else{
    mdat = read.csv(mds_mdat_fp)
    
    dd = mdat[,c("cellID","donorID","broadLineage",'finalAnn_broad','disease','UMAP_1','UMAP_2')]
    dd = dd[dd$finalAnn_broad != 'doublets',]
    dd = dd[sample(1:nrow(dd),nrow(dd)),]
  }
  
  dd$group = ifelse(dd$finalAnn_broad == 'MDS',paste0(dd$donorID,'_MDS'),dd$donorID)
  plotFun_MDS_UMAP = function(noFrame=FALSE,noPlot=FALSE){
    p = ggplot(dd,aes(UMAP_1,UMAP_2,col=group))+
      geom_point(data=dd[dd$group != 'L067_MDS',],size=0.1)+
      geom_point(data=dd[dd$group == 'L067_MDS',],size=0.3)+
      scale_color_manual(values = c('L061'=grey(0.5),
                                    'L067'=grey(0.8),
                                    'L061_MDS'='#cc5d0f',
                                    'L067_MDS'=col25[4]))+
      theme_classic(base_size = 11)+xlab('')+ylab('')+
      theme(panel.border = element_blank(),
            axis.line = element_blank(),legend.position = 'none',
            axis.ticks = element_blank(),
            axis.text = element_blank())
    
    print(p)
  }
  
  saveFig(file.path(plotDir,'Fig1_MDS_celltype_UMAP'),plotFun_MDS_UMAP,rawData=mdat,width = 3,height = 2.9,res = 500,useDingbats = F)
  
  
  broadLin_col=c('B_cells' = pal34H[34],
                 'Ery.MK' = col25[8],
                 'MDS' = 'darkred',
                 'Myeloid' = pal34H[9],
                 'NK.T' = pal34H[14]
  )
  plotFun_MDS_broadLin_UMAP = function(noFrame=FALSE,noPlot=FALSE){
    p = ggplot(dd,aes(UMAP_1,UMAP_2,col=broadLineage))+
      geom_point(size=0.1)+
      scale_color_manual(values = col25)+
      theme_classic(base_size = 11)+xlab('')+ylab('')+
      theme(panel.border = element_blank(),
            axis.line = element_blank(),legend.position = 'right',
            axis.ticks = element_blank(),
            axis.text = element_blank())
    
    print(p)
  }
  
  saveFig(file.path(plotDir,'Fig1_MDS_celltype_UMAP'),plotFun_MDS_UMAP,rawData=mdat,width = 3,height = 2.9,res = 500,useDingbats = F)
  
  
  
}





sgs = subset(srat,subset = donorID %in% c('BJ111','BJ112','BJ113','BJ114','BJ9','GOSH084','GSM5160432','GSM5160434','GSM5160435'))
sgs = subset(srat,subset = donorID %in% c('BJ111','BJ112','BJ113','BJ114','BJ9','GOSH084'))

sgs = standard_clustering(sgs)
#saveRDS(sgs,'~/lustre_mt22/SETBP1/Results/MDS_L061.L067_HARM_2502.RDS')
#saveRDS(sgs,'~/lustre_mt22/SETBP1/Results/SGS_scRNAseq_2502.RDS')
saveRDS(sgs,'~/lustre_mt22/SETBP1/Results/SGS.inhouseDataOnly_scRNAseq_2502.RDS')

sgs$disease = ifelse(sgs$donorID %in% c("GOSH084",'BJ113'),'aSGS',
                         ifelse(sgs$donorID %in% c("BJ9",'BJ111','BJ112'),'SGS',
                                ifelse(sgs$donorID %in% c('L061','L067'),'MDS',
                                       ifelse(sgs$donorID == 'BJ114','adult','children'))))

DimPlot(sgs,group.by = 'disease',cols = col25)

mdat = cbind(sgs@meta.data,sgs@reductions$umap@cell.embeddings)
#write.csv(mdat,'~/lustre_mt22/SETBP1/Results/SGS_scRNAseq_mdat_2502.csv')
write.csv(mdat,'~/lustre_mt22/SETBP1/Results/SGS.inhouseDataOnly_scRNAseq_mdat_2502.csv')

sgs_mdat_fp = '~/lustre_mt22/SETBP1/Results/SGS.inhouseDataOnly_scRNAseq_mdat_2502.csv'

fig1_SGS_UMAP = function(){
  
  if(file.exists(file.path(plotDir,'Fig1b_SETBP1.MDS_UMAP_rawData.tsv'))){
    dd = read.delim(file.path(plotDir,'Fig1b_SETBP1.MDS_UMAP_rawData.tsv'),header = T,sep = '\t')
    
  }else{
    mdat = read.csv(sgs_mdat_fp)
    
    dd = mdat[,c("cellID","donorID","broadLineage",'finalAnn_broad','disease','UMAP_1','UMAP_2')]
    dd = dd[dd$finalAnn_broad != 'doublets',]
    dd = dd[sample(1:nrow(dd),nrow(dd)),]
  }
  
  SETBP1_mutation_col = c('wild_type' = colAlpha(grey(0.6),0.3),
                          'hotspot' = colAlpha('#cc5d0f',0.55),
                          'non_hotspot' = colAlpha('#1b4891',0.55))
  dd$setbp1_mutation_status = ifelse(dd$donorID %in% c('BJ9','BJ111','BJ112','L061'),'hotspot',
                                     ifelse(dd$donorID %in% c('GOSH084','BJ113'),'non_hotspot','wild_type'))
  
  plotFun_setbp1_mutation = function(noFrame=FALSE,noPlot=FALSE){
    par(mar=c(0,0,0.1,0))
    
    plot(dd$UMAP_1,dd$UMAP_2,
         las=1,
         type='n',
         cex.main = 0.85,xaxt='n',yaxt='n',
         xlab='',ylab='',
         frame.plot=F)
    
    if(!noPlot){
      points(dd$UMAP_1,dd$UMAP_2,
             col = SETBP1_mutation_col[as.character(dd$setbp1_mutation_status)],
             pch = 19,
             cex=0.07)
    }
  }
  
  saveFig(file.path(plotDir,'Fig1_SGS.inhouseDataOnly_setbp1.mutation_UMAP'),plotFun_setbp1_mutation,rawData=mdat,width = 4.1,height = 3.7,res = 500,useDingbats = F)
  
  
  
  
  
  broadLin_col=c('B_cells' = pal34H[34],
                 'Ery.MK' = col25[8],
                 'MDS' = 'darkred',
                 'Myeloid' = pal34H[9],
                 'NK.T' = pal34H[14]
  )
  plotFun_SGS_broadLin_UMAP = function(noFrame=FALSE,noPlot=FALSE){
    p = ggplot(dd,aes(UMAP_1,UMAP_2,col=broadLineage))+
      geom_point(size=0.1)+
      scale_color_manual(values = col25)+
      theme_classic(base_size = 11)+xlab('')+ylab('')+
      theme(panel.border = element_blank(),
            axis.line = element_blank(),legend.position = 'right',
            axis.ticks = element_blank(),
            axis.text = element_blank())
    
    print(p)
  }
  
  saveFig(file.path(plotDir,'Fig2_SGS_celltype_UMAP'),plotFun_SGS_broadLin_UMAP,rawData=mdat,width = 3,height = 2.9,res = 500,useDingbats = F)
  
  
  
  
  
  
  donorID_col=c('PD60299' = '#6da0de',
                'PD66174' = pal34H[27],
                # SGS
                'PD60412' = brewer.pal(n=10,'Paired')[7],
                'PD66162' = pal34H[18],
                'PD66163' = pal34H[17],
                # Normal
                'PD66164' = grey(0.4),
                'GSM5160432' = grey(0.85),
                'GSM5160434' = grey(0.7),
                'GSM5160435' = grey(0.55),
                # MDS
                'PD61858' = pal34H[3],
                'PD61857' = pal37H[12]
  )
  
  #show_col(pal37H)
  
  dd$donorID[dd$donorID == 'GOSH084'] = 'PD60299'
  dd$donorID[dd$donorID == 'BJ113'] = 'PD66174'
  dd$donorID[dd$donorID == 'BJ9'] = 'PD60412'
  dd$donorID[dd$donorID == 'BJ111'] = 'PD66162'
  dd$donorID[dd$donorID == 'BJ112'] = 'PD66163'
  dd$donorID[dd$donorID == 'BJ114'] = 'PD66164'
  dd$donorID[dd$donorID == 'L061'] = 'PD61858'
  dd$donorID[dd$donorID == 'L067'] = 'PD61857'
  
  plotFun_SGS_donorID_UMAP = function(noFrame=FALSE,noPlot=FALSE){
    par(mar=c(0,0,0.1,0))
    
    plot(dd$UMAP_1,dd$UMAP_2,
         las=1,
         type='n',
         cex.main = 0.85,xaxt='n',yaxt='n',
         xlab='',ylab='',
         frame.plot=F)
    
    if(!noPlot){
      points(dd$UMAP_1,dd$UMAP_2,
             col = donorID_col[as.character(dd$donorID)],
             pch = 19,
             cex=0.07)
    }
  }
  
  saveFig(file.path(plotDir,'SupFig1_SGS_donorID_UMAP'),plotFun_SGS_donorID_UMAP,rawData=dd,width = 3.1,height = 3.1,res = 500,useDingbats = F)
  
  
}





markers = (c('PTPRC','SERPINB1', # Immune cells
             'CD34','HMGA2','MEG3','RBPMS', # MDS
             'TPSAB1','KIT', # MDS
            'MLLT3','PRSS57', # HSC_MPP
            'GATA2','GATA1',#MEMP
            'KLF1','HBD','GYPA',# early Erythroid
            'ALAS2', # mid.erythroid
            'HBA1','BPGM', # late.erythroid
            
            'ITGA2B', #Megakaryocyte
            'NRGN','PPBP','PF4','PLEK',  # Platelets
            
            'CD68','HMOX1','ITGAX',#Macrophages
            'ITGAM','CD14','FCGR3A',# Monocytes
            
            
            'IRF8',	#DC.precursor
            'FLT3', #DCs
            'CLEC9A',#DC1
            'CLEC10A','CD1C', # DC2 
            'IL3RA', # DC2
            'CLEC4C', #pDC
            
            
            
            'CD79B','CD19',# pro-b-cells
            # naive_b_cell = CD19+ CD27-
            # memory_b_cell = CD19+ and CD27 + 
            'CD79A',# pre-b-cells
            'TCL1A','MS4A1',  # B-cells
            
            'CD27', # Plasma cells
            
            #'LTB', 'CD52',# ILC precursor
            'IL7R','CD3D','CD8A', #Early.lymphoid_T.lymphocyte
            #'TRDV2','TRGV9', # gamma delta T-cell
            'PRF1','GZMA', # effector T cells
            'ZBTB16','NKG7','KLRD1' #NK
            
))

setbp1_srat_fp = '~/lustre_mt22/SETBP1/Results/x_published_scRNAseq_preprocessing/SETBP1_Wang21_scRNAseq.only_noUnknowns_2502.RDS'

supFig1_SETBP1.MDS_dotPlot = function(){
  library(viridis)
  
  # Import the seurat object
  srat = readRDS(setbp1_srat_fp)
  mdat = read.csv(setbp1_mdat_fp)
  
  srat$finalAnn_broad[srat$finalAnn_broad == 'strange_MonoMac'] = 'Mono_Mac'
  srat$finalAnn_broad = factor(srat$finalAnn_broad,c('MDS','Ery','MK',
                                                     'Mono_Mac','Mono_CD16','DC1','DC2','pDC',
                                                     'B_cells','Plasma.cell','T_cells','NK_T','NK'))
  direction = 'vertical'
  if(direction == 'vertical'){
    genes = rev(markers)
  }else{
    genes = markers
  }
  
  plotFun = function(noFrame=FALSE,noPlot=FALSE){
    Idents(srat) = srat$finalAnn_broad
    p = DotPlot(srat,#idents = unique(srat$finalAnn_broad),
                #cols = c("#EBFFE5", "#244b05"),
                #cols = c(colAlpha('#F1F5FA',1),'#425580'),
                cols = c(colAlpha(grey(0.95),0.8),'black'),
                #cols = c(grey(0.99), grey(0.2)),
                #group.by = 'seurat_clusters',
                #idents = unique(srat$finalAnn[srat$finalAnn != 'others']),
                features = genes)+RotatedAxis() 
    
    
    if(direction == 'vertical'){
      p = p + coord_flip() + 
        scale_y_discrete(position = "right")+
        theme(axis.text.x = element_text(size=9,angle = 90,vjust = 0.5,hjust = 0),
              axis.text.y = element_text(size=9),
              legend.title = element_text(size=8),
              legend.text = element_text(size=8),
              legend.position = 'top') + xlab('') + ylab('')
    }else if(direction == 'horizontal'){
      p = p +
        theme(axis.text.y = element_text(size=15),
              axis.text.x = element_text(size=8,angle = 90,vjust = 0.5,hjust = 1),
              legend.title = element_text(size=8),
              legend.text = element_text(size=8),
              legend.position = 'top') + xlab('') + ylab('')
    }
    print(p)
  }
  
  if(direction == 'vertical'){
    saveFig(file.path(plotDir,'SupFig1_SETBP1.MDS_Celltype_DotPlot_vertical'),plotFun,width = 4,height =10.5,res = 500)  
  }else{
    saveFig(file.path(plotDir,'SupFig1_SETBP1.MDS_Celltype_DotPlot_horizontal'),plotFun,width = 9.5,height = 4.6,res = 500)  
  }
  
}


supFig1_SETBP1_expression = function(){
  srat = readRDS(setbp1_srat_fp)
  srat$finalAnn_broad[srat$finalAnn_broad == 'strange_MonoMac'] = 'Mono_Mac'
  srat$finalAnn_broad = factor(srat$finalAnn_broad,c('MDS','Ery','MK',
                                                     'Mono_Mac','Mono_CD16','DC1','DC2','pDC',
                                                     'B_cells','Plasma.cell','T_cells','NK_T','NK'))
  srat$broadLineage = factor(srat$broadLineage,c('MDS','Ery.MK','Myeloid','B_cells','NK.T'))
  
  srat$group = paste0(srat$finalAnn_broad,'.',srat$donorID)
  avgExpr = AverageExpression(srat,features = 'SETBP1',group.by = 'group')
  avgExpr = as.data.frame(t(avgExpr$RNA))
  avgExpr$group = rownames(avgExpr)
  avgExpr$donorID = srat$donorID[match(avgExpr$group,srat$group)]
  avgExpr$disease = ifelse(avgExpr$donorID %in% c("GOSH084",'BJ113'),'aSGS',
                           ifelse(avgExpr$donorID %in% c("BJ9",'BJ111','BJ112'),'SGS',
                                  ifelse(avgExpr$donorID %in% c('L061','L067'),'MDS',
                                         ifelse(avgExpr$donorID == 'BJ114','adult','children'))))
  
  avgExpr$condition = srat$condition[match(avgExpr$group,srat$group)]
  avgExpr$celltype = srat$finalAnn_broad[match(avgExpr$group,srat$group)]
  avgExpr$celltype = factor(avgExpr$celltype,c('MDS','Ery','MK',
                                               'Mono_Mac','Mono_CD16','DC1','DC2','pDC',
                                               'B_cells','Plasma.cell','T_cells','NK_T','NK'))
  colnames(avgExpr)[1] = 'expr'
  
  
  
  setbp1_expr = as.data.frame(srat@assays$RNA@data['SETBP1',])
  colnames(setbp1_expr) = 'norm_expr'
  setbp1_expr$group = srat$group[match(rownames(setbp1_expr),srat$cellID)]
  setbp1_expr$celltype = srat$finalAnn_broad[match(rownames(setbp1_expr),srat$cellID)]
  setbp1_expr$celltype = factor(setbp1_expr$celltype,c('MDS','Ery','MK',
                                                       'Mono_Mac','Mono_CD16','DC1','DC2','pDC',
                                                       'B_cells','Plasma.cell','T_cells','NK_T','NK'))
  setbp1_expr$group_avg = avgExpr$expr[match(setbp1_expr$group,avgExpr$group)]
  setbp1_expr = setbp1_expr %>% group_by(group) %>% mutate(nCell = sum(norm_expr > 0),
                                                              frac = nCell / length(norm_expr))
  setbp1_expr$broadLineage = srat$broadLineage[match(setbp1_expr$celltype,srat$finalAnn_broad)]
  setbp1_expr$donorID = srat$donorID[match(setbp1_expr$group,srat$group)]
  
  avgExpr$frac_express = setbp1_expr$frac[match(avgExpr$group,setbp1_expr$group)]
  avgExpr$broadLineage = setbp1_expr$broadLineage[match(avgExpr$celltype,setbp1_expr$celltype)]
  
  setbp1_expr.sub = do.call(rbind,lapply(split(setbp1_expr,setbp1_expr$group),function(x){
    if(nrow(x) == 1){
      return(x)
    }else{
      return(x[sample(1:nrow(x),max(2,nrow(x)*0.05)),])
    }
  }))
  
  broadLin_col=c('B_cells' = pal34H[34],
                 'Ery.MK' = col25[8],
                 'MDS' = 'darkred',
                 'Myeloid' = pal34H[9],
                 'NK.T' = pal34H[14]
  )
  
  avgExpr$donorID[avgExpr$donorID == 'GOSH084'] = 'PD60299'
  avgExpr$donorID[avgExpr$donorID == 'BJ113'] = 'PD66174'
  avgExpr$donorID[avgExpr$donorID == 'BJ9'] = 'PD60412'
  avgExpr$donorID[avgExpr$donorID == 'BJ111'] = 'PD66162'
  avgExpr$donorID[avgExpr$donorID == 'BJ112'] = 'PD66163'
  avgExpr$donorID[avgExpr$donorID == 'BJ114'] = 'PD66164'
  avgExpr$donorID[avgExpr$donorID == 'L061'] = 'PD61858'
  avgExpr$donorID[avgExpr$donorID == 'L067'] = 'PD61857'
  
  avgExpr$disease[avgExpr$disease %in% c('adult','children')] = 'healthy'
  avgExpr$disease = factor(avgExpr$disease,c('MDS','SGS','aSGS','healthy'))
  
  donorID_col=c('PD60299' = '#8cc9ed',
                'PD66174' = pal34H[27],
                # SGS
                'PD60412' = brewer.pal(n=10,'Paired')[7],
                'PD66162' = pal34H[18],
                'PD66163' = pal34H[17],
                # Normal
                'PD66164' = grey(0),
                'GSM5160432' = grey(0.8),
                'GSM5160434' = grey(0.65),
                'GSM5160435' = grey(0.4),
                # MDS
                'PD61858' = pal34H[3],
                'PD61857' = pal37H[12]
  )
  plotFun_setbp1_expressionLevel = function(noFrame=FALSE,noPlot=FALSE){
    p1 = ggplot(setbp1_expr.sub,aes(celltype))  +
      geom_quasirandom(aes(y=norm_expr),size=0.5,width = 0.3,alpha=0.3,col=grey(0.8))+
      geom_point(data=avgExpr,aes(y=expr,size = frac_express,col=donorID))+
      #scale_size(breaks = seq(0,0.45,0.1), limits = c(0,0.45)) +
      scale_color_manual(values = donorID_col,name='Donor ID')+
      facet_grid(.~broadLineage,scales = 'free_x',space = 'free_x')+
      theme_classic(base_size = 13)+
      ggtitle('')+xlab('')+ylab('Normalised SETBP1 expression level')+
      theme(panel.border = element_rect(fill=F,colour = 'black'),axis.line = element_blank(),
            strip.background=element_blank(),
            strip.text.x = element_text(size=10,colour = 'black'),
            strip.text.y = element_text(size=10,colour = 'black'),
            axis.ticks = element_line(colour = 'black',linewidth = 0.2),
            axis.text.x = element_text(size = 10,angle = 90, vjust = 0.5,hjust = 1,colour = 'black'),
            axis.text.y = element_text(size=10),
            axis.text = element_text(colour = 'black'),
            legend.text = element_text(colour = 'black',size = 10),
            legend.title = element_text(colour = 'black',size = 12))
    
    
    
    avgExpr = avgExpr[!(avgExpr$disease == 'MDS' & avgExpr$celltype != 'MDS'),]
    p1 = ggplot(avgExpr,aes(disease,expr))  +
      geom_boxplot(outlier.shape = NA)+
      geom_quasirandom(aes(size = frac_express,col=donorID),alpha=0.9,width = 0.1)+
      #scale_size(breaks = seq(0,0.45,0.1), limits = c(0,0.45)) +
      scale_color_manual(values = donorID_col,name='Donor ID')+
      facet_grid(.~broadLineage+celltype,scales = 'free_x',space = 'free_x')+
      theme_classic(base_size = 13)+
      ggtitle('')+xlab('')+ylab('Average normalised SETBP1 expression')+
      theme(panel.border = element_rect(fill=F,colour = 'black'),axis.line = element_blank(),
            strip.background=element_blank(),
            strip.text.x = element_text(size=10,colour = 'black'),
            strip.text.y = element_text(size=10,colour = 'black'),
            axis.ticks = element_line(colour = 'black',linewidth = 0.2),
            axis.text.x = element_text(size = 10,angle = 90, vjust = 0.5,hjust = 1,colour = 'black'),
            axis.text.y = element_text(size=10),
            axis.text = element_text(colour = 'black'),
            legend.text = element_text(colour = 'black',size = 10),
            legend.title = element_text(colour = 'black',size = 12))
    
    print(p1)
  }
  
  saveFig(file.path(plotDir,'SupFig1_SETBP1_expressionLevel'),plotFun_setbp1_expressionLevel,rawData = setbp1_expr,width = 6,height =4,res = 500)  
  saveFig(file.path(plotDir,'SupFig1_SETBP1_expressionLevel_v2'),plotFun_setbp1_expressionLevel,rawData = setbp1_expr,width = 11,height =5,res = 500)  

}



fig1_SETBP1.hs_Leuk.signature = function(){
  
  ## 1. Import list of aCML - DEGs
  bulk_DEGs_aCML = readRDS('~/lustre_mt22/SETBP1/Results/4_SETBP1_impact_inLeuk/feb25/aCML_hs.vs.wt_edgeR.RDS')
  tt_aCML = bulk_DEGs_aCML[['tt']]
  tt_aCML = tt_aCML[order(abs(tt_aCML$logFC),decreasing = T),]
  tt_aCML$direction = ifelse(tt_aCML$logFC > 0,'SETBP1hs_up','SETBP1hs_down')
  table(tt_aCML$direction)
  dd = as.data.frame(table(tt_aCML$direction))
  dd$Freq[dd$Var1 == 'SETBP1hs_down'] = -dd$Freq[dd$Var1 == 'SETBP1hs_down']
  
  ## Plot bar plot of number of DEGs
  plotFun_nDEG = function(noFrame=FALSE,noPlot=FALSE){
    
    p = ggplot(dd,aes(Var1,Freq,fill=Var1))+
      geom_col(width = 0.58)+
      #scale_fill_manual(values = c(col25[1],col25[2],colAlpha(col25[1:2],alphas = 0.6)))+
      scale_fill_manual(values = c('#1a4a87','#a4282c'),name='')+
      geom_hline(yintercept = 0,col='black',lwd=0.5)+
      scale_y_continuous(breaks = c(-100,0,200,400,600),limits = c(-110,620)) +
      # scale_y_break(c(-40,-1870),scales = 4) +
      theme_classic(base_size = 11)+xlab('')+ylab('# DEGs')+
      theme(panel.border = element_rect(fill=F,colour = 'black',linewidth = 1),
            axis.line = element_blank(),#legend.position = 'bottom',
            axis.ticks = element_line(colour = 'black'),
            legend.title = element_text(size=10,colour = 'black'),
            legend.text = element_text(size=8,colour = 'black'),legend.key.size = unit(0.5,'cm'),
            axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,colour = 'black'),
            axis.text = element_text(color='black'))
    print(p)
    
  }
  
  saveFig(file.path(plotDir,paste0('Fig1_aCML_edgeR.DEGs')),plotFun_nDEG,rawData=dd,width = 3.15,height = 2.9,res = 500,useDingbats = F)
  
}
