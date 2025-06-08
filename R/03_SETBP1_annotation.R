# Annotate SETBP1 scRNAseq dataset, using:
# 1. CellTypist
# 2. Manual assessment of marker genes expression

setwd('~/SETBP1_GoF_mutations')
library(Seurat)
library(Matrix)
library(tidyverse)
source("R/utils/logisticRegressionCellTypist.R")

outDir = 'Results/03_SETBP1_annotation/2505'
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}

## Import SETBP1 srat object ------------------------------------------------------------------
srat = readRDS('Results/02_SETBP1_scQC/2505/SETBP1_SGS_clean_noMTCells.RDS')


## CellTypist Logistic Regression ------------------------------------------------------------- 
## using CellTypist's default lowLevel Immume model
mtx = srat@assays$RNA@counts
message('Predicting cell labels....')
lr_output_fp = file.path(outDir,'lowLevelImmumeREF_LRwCT_output.RDS')
if(!file.exists(lr_output_fp)){
  lr_output = runCelltypist(cnts=mtx)
  saveRDS(lr_output,lr_output_fp)
  
  message('Doing plots....')
  pdf(gsub('output.RDS','logit.pdf',lr_output_fp),width = 30,height = 30)
  library(ComplexHeatmap)
  p = plotByGroup(fit=lr_output[['logitMat']],cellGroupings=lr_output[['labMat']][,2],colSpaceTruncate=TRUE)
  p
  p = plotByGroup(fit=lr_output[['logitMat']],cellGroupings=srat$seurat_clusters[match(rownames(lr_output[[1]]),rownames(srat@meta.data))],colSpaceTruncate=TRUE)
  p
  dev.off()
  
}else{
  lr_output = readRDS(lr_output_fp)
}


message('Calculating softmax')
## Softmax to assign low_conf cells
softmax_p = lr_output[['logitMat']]
softmax_p = t(apply(softmax_p,1,function(x){sapply(x,function(xi){xi = exp(xi)/sum(exp(x))})}))


## Assign labels to each cell
#  for each cell, find the name with the best match - ie. highest probability
#  assign the celltype if prob >=0.95
#  assign as low_conf if prob >=0.8
#  else - unknown
softmax_label = apply(softmax_p, 1, function(x){
  max_score = max(x)
  second_max_score = max(x[x!=max_score])
  diff = max_score - second_max_score
  
  
  lab = ifelse(sum(x == max(x)) > 1,
               paste(c('ambiguous',names(x[x==max(x)])),sep = ':'),
               ifelse(max(x) >= 0.95,
                      ifelse(diff >= 0.1, names(x[x==max_score]), 
                             paste(c('ambiguous',names(x[x==max_score]),names(x[x==second_max_score])),sep = ':')),
                      ifelse(max(x) >= 0.85,
                             ifelse(diff >= 0.1, paste0('lowConf_',names(x[x==max_score])), 
                                    paste(c('lowConf_ambiguous',names(x[x==max_score]),names(x[x==second_max_score])),sep = ':')),
                             'unknown')))
  return(lab)
})



#### Add annotation to the new object
srat[[paste0('LR_predicted_label_lowImm')]] = lr_output[[2]]$predicted_labels[match(rownames(srat@meta.data),lr_output[[2]]$X)]
srat[[paste0('LR_softmax_predicted_label_lowImm')]] = softmax_label[match(rownames(srat@meta.data),names(softmax_label))]

write.csv(srat@meta.data,file.path(outDir,'SETBP1_SGS_clean_noMTCells_annot_wCellTypist_2505.csv'))

## Finalise the annotation
srat$LR_predicted_label_lowImm2 = gsub(' |-|\\+|\\/','.',srat$LR_predicted_label_lowImm)
a = as.data.frame(table(srat$LR_predicted_label_lowImm,srat$LR_softmax_predicted_label_lowImm))
a$Var1 = gsub(' |-|\\+|\\/','.',a$Var1)
a = a[a$Freq > 0,]
a = a[a$Var1 != a$Var2,]
colnames(a) = c('predicted','soft_max','Freq')


DimPlot(srat,group.by = 'seurat_clusters',repel = T,label = T,label.box = T)
DimPlot(srat,group.by = 'annot_tmp',repel = T,label = T,label.box = T) + NoLegend()

DimPlot(srat,cells.highlight = srat$cellID[srat$LR_predicted_label_lowImm2=='Mast.cells' & 
                                             srat$LR_softmax_predicted_label_lowImm == 'unknown'])

DimPlot(srat,cells.highlight = srat$cellID[srat$annot_tmp == 'Mast.cells'] )

srat$annot_tmp = as.character(srat$seurat_clusters)
srat$annot_tmp[srat$annot_tmp == '28'] = '28_doublets?'
srat$annot_tmp[srat$LR_predicted_label_lowImm2 == srat$LR_softmax_predicted_label_lowImm] = srat$LR_predicted_label_lowImm2[srat$LR_predicted_label_lowImm2 == srat$LR_softmax_predicted_label_lowImm]
srat$annot_tmp[srat$annot_tmp == 'Pro.B.cells'] = 'Pro.B.cells_doublets?'

srat$broadAnno = as.character(srat$annot_tmp)
srat$broadAnno[srat$broadAnno == 'Megakaryocytes.platelets'] = 'MegK'
srat$broadAnno[srat$broadAnno == 'Tcm.Naive.cytotoxic.T.cells'] = 'Tcell_cytotoxic'
srat$broadAnno[srat$broadAnno == 'CD16..NK.cells'] = 'NK_CD16+'
srat$broadAnno[srat$broadAnno == 'Tem.Effector.helper.T.cells'] = 'Tcell_helper'
srat$broadAnno[srat$broadAnno == 'Classical.monocytes'] = 'Monocyte_CD14'
srat$broadAnno[srat$broadAnno == 'Non.classical.monocytes'] = 'Monocyte_CD16'

srat$broadAnno[srat$broadAnno == 'Tcm.Naive.helper.T.cells'] = 'Tcell_helper_naive'
srat$broadAnno[srat$broadAnno == 'Regulatory.T.cells'] = 'Tcell_regulatory'
srat$broadAnno[srat$broadAnno == 'MAIT.cells'] = 'MAIT.cell'
srat$broadAnno[srat$broadAnno == 'Plasma.cells'] = 'Plasma.cell'
srat$broadAnno[srat$broadAnno == 'Naive.B.cells'] = 'B.cell_naive'
srat$broadAnno[srat$broadAnno == 'Tem.Temra.cytotoxic.T.cells'] = 'Tcell_cytotoxic_temra?'
srat$broadAnno[srat$broadAnno == 'Tem.Trm.cytotoxic.T.cells'] = 'Tcell_cytotoxic_trm?'



srat = readRDS(file.path(outDir,'SETBP1_SGS_clean_noMTCells_annot_2505.RDS'))






markers = c('PTPRC','SERPINB1', # Immune cells
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
)

DotPlot(srat,group.by = 'annot_tmp',features = markers)+
  RotatedAxis() + 
  ggtitle('Celltype marker genes') + xlab('')+
  theme(axis.text.x = element_text(size=8,angle=90,vjust = 0.5,hjust = 1))
