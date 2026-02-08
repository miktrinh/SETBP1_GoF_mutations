

##---------------##
#### Libraries ####
##---------------##


library(Seurat)
library(GenomicFeatures)
library(DESeq2)
library(ComplexHeatmap)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(readxl)
library(ggsci)
source("~/lustre_mt22/generalScripts/utils/misc.R")
source("~/lustre_mt22/generalScripts/utils/sc_utils.R")
source("~/lustre_mt22/generalScripts/utils/pseudobulk.R")



##-----------------------------------------------------------------------##
#### 1a. Import relevant datasets: TF + membrane protein + Gene map    ####
##-----------------------------------------------------------------------##

#Define genomic coordinates
gtf = '/nfs/srpipe_references/downloaded_from_10X/refdata-b37-2.1.0/genes/gene_annotations.gtf.gz'
txdb = makeTxDbFromGFF(gtf)
gns = genes(txdb)
geneMap$chr = as.character(seqnames(gns[match(geneMap$ensID,gns$gene_id)]))



outDir = '~/lustre_mt22/SETBP1/Results/8_bulkRNA_SETBP1'
##--------------------------------------------##
##        DEG on BEAT_AML_v2 dataset        ####
##--------------------------------------------##
beat_sce_path = file.path(outDir,'BEAT_AML_sub_sce.RDS')
beat_sce = readRDS(beat_sce_path)
beat_mdat = as.data.frame(colData(beat_sce))
beat_rowData = as.data.frame(rowData(beat_sce))
beat_rowData$chr = as.character(seqnames(gns[match(beat_rowData$ensembl_gene_id,gns$gene_id)]))
colnames(beat_rowData)[1:2] = c('ensID','geneSym')

## 1. Prepare raw count matrix for DESeq2
# Raw count matrix
beat_rawCnt = assays(beat_sce)[['counts_raw']]

# Drop bullshit and Uninteresting genes
# Drop genes NOT on the main chromosomes
toc = beat_rawCnt[rownames(beat_rawCnt) %in% beat_rowData$ensID[beat_rowData$chr %in% c(1:23,'X')],!colnames(beat_rawCnt) %in% c('BA3212R','BA3300R')]
# Drop Non-expressed genes
w = which(rowSums(toc)>0)
toc = toc[w,]

## 2. Prepare colData for DESeq2
colDat = beat_mdat[match(colnames(toc),beat_mdat$sampleID),]
rownames(colDat) = colDat$sampleID
colDat$setbp1_mut_inHotspot[is.na(colDat$setbp1_mut_inHotspot)] = 'none'

colDat$group = ifelse(grepl('Healthy',colDat$cancerType),'healthy',
                      ifelse(colDat$setbp1_mut_inHotspot == 'none','diseased_none',colDat$setbp1_mut_inHotspot))
                      
colDat$group = factor(colDat$group,levels = c('healthy','diseased_none','outside_hotspot','inside_hotspot'))                            
colDat$clin_age = ifelse(is.na(colDat$clin_age),'unknown',colDat$clin_age)


## 3. DEGs

beat_res = bulkDESeq2(toc = toc, 
                      mDat = colDat[match(colnames(toc),colDat$sampleID),c('sampleID','age','cancerType','cancerSubClass','group')],
                      cellTypeName = 'BEAT_AML',
                      formula = '~ %s',donorID = 'sampleID',groupID = 'group',
                      geneMap = beat_rowData)


outDir = '~/lustre_mt22/SETBP1/Results/13_bulkDEGs_BEATaml/aug23'
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}
saveRDS(beat_res,file.path(outDir,'bulk_DESEQ2_BEATaml_output.RDS'))
beat_res = readRDS(file.path(outDir,'bulk_DESEQ2_BEATaml_output.RDS'))

beat_results = as.data.frame(beat_res[['res']])
beat_results = annotateGenes(beat_results,geneMap = beat_rowData)
table(beat_results$contrast)


beat_deg = beat_res[['all_de']]
table(beat_deg$contrast)


write.csv(beat_deg,file.path(outDir,'bulk_DESEQ2_BEATaml_deg.csv'))
all_de = read.csv(file.path(outDir,'bulk_DESEQ2_BEATaml_deg.csv'))
## add direction to DEG
all_de$direction = ifelse(all_de$log2FoldChange >0,'down_HCA','up_HCA')
all_de$list_id = paste0(all_de$contrast,':',all_de$direction)








#----------------------------------
# DEGs all AML vs Healthy
colDat$group = ifelse(grepl('Healthy',colDat$cancerType),'healthy','diseased')
colDat$group = factor(colDat$group,levels = c('healthy','diseased'))                            

## 3. DEGs

beat_res = bulkDESeq2(toc = toc, 
                      mDat = colDat[match(colnames(toc),colDat$sampleID),c('sampleID','age','cancerType','cancerSubClass','group')],
                      cellTypeName = 'BEAT_AML',
                      formula = '~ %s',donorID = 'sampleID',groupID = 'group',
                      geneMap = beat_rowData)
saveRDS(beat_res,file.path(outDir,'bulk_DESEQ2_BEATaml_allAMLvsAllHealthyoutput.RDS'))
beat_res = readRDS(file.path(outDir,'bulk_DESEQ2_BEATaml_allAMLvsAllHealthyoutput.RDS'))

beat_results = as.data.frame(beat_res[['res']])
beat_results = annotateGenes(beat_results,geneMap = beat_rowData)
table(beat_results$contrast)


beat_deg = beat_res[['all_de']]
table(beat_deg$contrast)

write.csv(beat_deg,file.path(outDir,'bulk_DESEQ2_BEATaml_allAMLvsAllHealthy_deg.csv'))
all_de = read.csv(file.path(outDir,'bulk_DESEQ2_BEATaml_allAMLvsAllHealthy_deg.csv'))
## add direction to DEG
all_de$direction = ifelse(all_de$log2FoldChange >0,'down_HCA','up_HCA')
all_de$list_id = paste0(all_de$contrast,':',all_de$direction)


#---------------------------------------
# DEGs BEAT_AML SETBP1.hs vs SETBP1.none
colDat$group = ifelse(grepl('Healthy',colDat$cancerType),'healthy',as.character(colDat$setbp1_mut_inHotspot))
colDat$group = factor(colDat$group,levels = c('none','inside_hotspot','outside_hotspot','healthy'))

## 3. DEGs

beat_res = bulkDESeq2(toc = toc,
                      mDat = colDat[match(colnames(toc),colDat$sampleID),c('sampleID','age','cancerType','cancerSubClass','group')],
                      cellTypeName = 'BEAT_AML',
                      formula = '~ %s',donorID = 'sampleID',groupID = 'group',
                      geneMap = beat_rowData)
saveRDS(beat_res,file.path(outDir,'bulk_DESEQ2_BEATaml_hs.vs.noneDiseased_output.RDS'))
beat_res = readRDS(file.path(outDir,'bulk_DESEQ2_BEATaml_hs.vs.noneDiseased_output.RDS'))

beat_results = as.data.frame(beat_res[['res']])
beat_results = annotateGenes(beat_results,geneMap = beat_rowData)
table(beat_results$contrast)


beat_deg = beat_res[['all_de']]
table(beat_deg$contrast)


write.csv(beat_deg,file.path(outDir,'bulk_DESEQ2_BEATaml_hs.vs.noneDiseased_deg.csv'))
all_de = read.csv(file.path(outDir,'bulk_DESEQ2_BEATaml_hs.vs.noneDiseased_deg.csv'))
## add direction to DEG
all_de$direction = ifelse(all_de$log2FoldChange >0,'down_HCA','up_HCA')
all_de$list_id = paste0(all_de$contrast,':',all_de$direction)



#----------------------------------
## 3b. DEGs against CD34+ healthy cells only (above, the reference is Healthy_CD34 + healthy_BM_MNC)
colDat$group = ifelse(grepl('Healthy.*CD34',colDat$cancerType),'healthy_CD34pos',
                      ifelse(colDat$cancerType == 'Healthy Individual BM MNC','healthy_bmmnc',
                             ifelse(colDat$setbp1_mut_inHotspot == 'none','diseased_none',colDat$setbp1_mut_inHotspot)))

colDat$group = factor(colDat$group,levels = c('healthy_CD34pos','healthy_bmmnc','diseased_none','outside_hotspot','inside_hotspot'))                            

beat_res = bulkDESeq2(toc = toc, 
                      mDat = colDat[match(colnames(toc),colDat$sampleID),c('sampleID','age','cancerType','cancerSubClass','group')],
                      cellTypeName = 'BEAT_AML',
                      formula = '~ %s',donorID = 'sampleID',groupID = 'group',
                      geneMap = beat_rowData)

outDir = '~/lustre_mt22/SETBP1/Results/13_bulkDEGs_BEATaml'
saveRDS(beat_res,file.path(outDir,'bulk_DESEQ2_BEATaml_healthyCD34ref_output.RDS'))
beat_res = readRDS(file.path(outDir,'bulk_DESEQ2_BEATaml_healthyCD34ref_output.RDS'))

beat_results = as.data.frame(beat_res[['res']])
beat_results = annotateGenes(beat_results,geneMap = beat_rowData)
table(beat_results$contrast)


beat_deg = beat_res[['all_de']]
table(beat_deg$contrast)


write.csv(beat_deg,file.path(outDir,'bulk_DESEQ2_BEATaml_healthyCD34ref_deg.csv'))
all_de = read.csv(file.path(outDir,'bulk_DESEQ2_BEATaml_healthyCD34ref_deg.csv'))
## add direction to DEG
all_de$direction = ifelse(all_de$log2FoldChange >0,'down_HCA','up_HCA')
all_de$list_id = paste0(all_de$contrast,':',all_de$direction)

table(res$padj < 0.05)
res = all_de[all_de$contrast == 'group_inside_hotspot_vs_healthy_CD34pos',]
plot(res$log2FoldChange,-log10(res$padj),
     col = ifelse(res$padj< 0.01 & abs(res$log2FoldChange)>1.5,'red','black'),
     cex=0.5,
     xlab='logFC',
     ylab='-log10(adjusted p-value)',
     main='SETBP1.hs AML vs CD34+ Healthy BM'
)

table(res$log2FoldChange[res$padj< 0.01 & abs(res$log2FoldChange)>1.5]<0)

#----------------------------------










# From here on, copy the code from 8_curate_SETBP1_bulkRNA
de_list = split(all_de$geneSym,all_de$list_id)

library(UpSetR)
upset(data=fromList(de_list),nsets = 10)


all_de = all_de %>% group_by(geneSym,direction) %>% mutate(n_contrast = n_distinct(contrast),
                                                           nDiseasedNone = ifelse(sum(grepl('diseased_none',unique(contrast)))>0,T,F),
                                                           nInHotspot = ifelse(sum(grepl('inside_hotspot',unique(contrast)))>0,T,F),
                                                           nOutHotspot = ifelse(sum(grepl('outside_hotspot',unique(contrast)))>0,T,F))
# a = all_de %>% group_by(geneSym) %>% summarise(n = n_distinct(list_id))

#### SETBP1 mutation signature DE genes ####
# these are DEGs (up and down) which are common between hotspot/non_hotspot vs healthy, but NOT in diseased_none

setbp1_sig_up = intersect(de_list[['group_inside_hotspot_vs_healthy:down_HCA']],de_list[['group_outside_hotspot_vs_healthy:down_HCA']])
setbp1_sig_up = setbp1_sig_up[!setbp1_sig_up %in% de_list[['group_diseased_none_vs_healthy:down_HCA']]]


setbp1_sig = all_de[all_de$n_contrast == 2 & !all_de$nDiseasedNone,]
table(setbp1_sig$list_id)

devtools::install_github("gaospecial/ggVennDiagram")
library("ggVennDiagram")

# Default plot
ggVennDiagram(de_list[grepl('down_HCA$',names(de_list))],category.names = c('nonSETBP1','Hotspot','nonHotspot'),edge_size = 0)+
  ggplot2::scale_fill_gradient(low="lightblue",high = "darkblue")


#########-------------- <- ########
## DEGs from sc dataset
all_deg.sc = read.csv('~/lustre_mt22/SETBP1/Results/3_pseudobulk_DESeq2/jun23/DESeq2_byCT_byDataset_allDEGs.csv')
#all_deg.sc$batch = gsub('dataset_broad_HCA_vs_','',all_deg.sc$contrast)
all_deg.sc$batch = gsub('_vs_healthyControl|group_','',all_deg.sc$contrast)
all_deg.sc$list_id = paste0(all_deg.sc$celltype,':',all_deg.sc$batch,':',all_deg.sc$direction)
all_deg.sc = all_deg.sc[all_deg.sc$DE ==T,]


##------------------------------------------------------------------##
##  Extract common and unique DEGs per CT between GOSH84 and BJ9  ####
##------------------------------------------------------------------##
all_deg.sc = all_deg.sc %>% group_by(geneSym,direction) %>% mutate(nCT = n(),
                                                                   nCT_bj9 = n_distinct(list_id[grepl('BJ9',list_id)]),
                                                                   nCT_gosh84 = n_distinct(list_id[grepl('GOSH84',list_id)]),
                                                                   nDataset = n_distinct(batch))

all_deg.sc = all_deg.sc %>% group_by(geneSym,direction,ct) %>% mutate(nDataset_byCT = n_distinct(batch))





# Default plot
all_deg.sc$list_id = paste0(all_deg.sc$batch,':',all_deg.sc$direction)
de_list.sc = lapply(split(all_deg.sc$geneSym,all_deg.sc$list_id),function(s){unique(s)})

ggVennDiagram(de_list.sc[grepl('down_HCA$',names(de_list.sc))],category.names = c('BJ9','GOSH84'),edge_size = 0)+
  ggplot2::scale_fill_gradient(low="lightblue",high = "darkblue")


x = append(de_list,de_list.sc)
ggVennDiagram(x[grepl('down_HCA$',names(x))],category.names = c('nonSETBP1','Hotspot','nonHotspot', 'BJ9','GOSH84'),edge_size = 0,
              #fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
              stroke_size = 0.5, set_name_size = 4)+
  ggplot2::scale_fill_gradient(low="lightblue",high = "darkblue")




## Hotspot mutation effect
x = list(hotspot = unique(de_list[['group_inside_hotspot_vs_healthy:down_HCA']]),
         GOSH84 = unique(de_list.sc[['dataset2_GOSH84_vs_HCA:down_HCA']]),
         BJ9 = unique(de_list.sc[['dataset2_BJ9_vs_HCA:down_HCA']]))

ggVennDiagram(x,category.names = c('hotspot', 'GOSH84','BJ9'),edge_size = 0,
              #fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
              stroke_size = 0.5, set_name_size = 4)+
  ggplot2::scale_fill_gradient(low="white",high = "darkred")

hotspot_DEGs_bulk = unique(de_list[['group_inside_hotspot_vs_healthy:down_HCA']])
hotspot_DEGs_bulk = hotspot_DEGs_bulk[!hotspot_DEGs_bulk %in% c(unique(de_list[['group_diseased_none_vs_healthy:down_HCA']]),
                                                               unique(de_list[['group_outside_hotspot_vs_healthy:down_HCA']]))]
sc_bj9 = x[['BJ9']]
sc_bj9 = sc_bj9[!sc_bj9 %in% x[['GOSH84']]]
hotspot_DEGs = intersect(hotspot_DEGs_bulk,sc_bj9)
hotspot_DEGs = all_deg.sc[all_deg.sc$geneSym %in% hotspot_DEGs & all_deg.sc$log2FoldChange>0,]
table(hotspot_DEGs$ct)
table(hotspot_DEGs$ct)


## Which of SETBP1 target genes were also DEGs
library(readxl)
piazza_18_genes = read_excel('~/lustre_mt22/SETBP1/Piazza_etal_2018/41467_2018_4462_MOESM9_ESM.xlsx',skip = 1)
View(hotspot_DEGs[hotspot_DEGs$geneSym %in% piazza_18_genes$`Approved Symbol`,])
DotPlot(combSrat_harm,group.by = 'ann',scale = T,
        features = piazza_18_genes$`Approved Symbol`[piazza_18_genes$log2FoldChange>0])+theme(axis.text = element_text(size=8)) + RotatedAxis()





## ENRICHR on this set of genes
# EnrichR
## Enrichment analysis ###
library(enrichR)
listEnrichrSites()
setEnrichrSite("Enrichr") # Human genes
websiteLive <- TRUE
dbs <- listEnrichrDbs()
if (is.null(dbs)) websiteLive <- FALSE
if (websiteLive) head(dbs)

dbs <- c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021", "GO_Biological_Process_2021",
         'HDSigDB_Human_2021','KEGG_2021_Human','MSigDB_Hallmark_2020','MSigDB_Oncogenic_Signatures',
         'TF_Perturbations_Followed_by_Expression','WikiPathways_2019_Human',
         'Cancer_Cell_Line_Encyclopedia','CCLE_Proteomics_2020','CellMarker_Augmented_2021',
         'Disease_Perturbations_from_GEO_down','Disease_Perturbations_from_GEO_up',
         'Drug_Perturbations_from_GEO_2014','Drug_Perturbations_from_GEO_down','Drug_Perturbations_from_GEO_up','DrugMatrix','IDG_Drug_Targets_2022',
         'Elsevier_Pathway_Collection','Enrichr_Submissions_TF-Gene_Coocurrence','Gene_Perturbations_from_GEO_down','Gene_Perturbations_from_GEO_up')

if (websiteLive) {
  up <- enrichr(hotspot_DEGs$geneSym[hotspot_DEGs$log2FoldChange >0 & hotspot_DEGs$celltype %in% c('Mono_CD14','Mono_CD16','DC2') & grepl('BJ9',hotspot_DEGs$contrast)], dbs)
  
  #down <- enrichr(up_MEP_notProB$geneSym[up_MEP_notProB$log2FoldChange <0], dbs)
}

up[[9]][c(1,3,4,5),]
plotEnrich(up[[9]], showTerms = 50, numChar = 40, y = "Count", orderBy = "P.value")
up.proB[[5]][c(1,2,6,12,14,16),]
up.proB[[9]][c(1,2,5,6,7,8,10,11,15),]
plotEnrich(up.proB[[9]], showTerms = 50, numChar = 40, y = "Count", orderBy = "P.value")

# Which of these genes were also DEGs in GOSH84
a = all_deg.sc[all_deg.sc]



























#####==================== aCML ====================≠≠≠###
##      DEGs of aCML dataset      ###
##--------------------------------------------##
##        DEG on aCML dataset        ####
##--------------------------------------------##
aCML_sce_path = file.path('~/lustre_mt22/SETBP1/Data/aCML_sce.RDS')
aCML_sce = readRDS(aCML_sce_path)
aCML_mdat = as.data.frame(colData(aCML_sce))
aCML_rowData = as.data.frame(rowData(aCML_sce))
#aCML_rowData$chr = as.character(seqnames(gns[match(aCML_rowData$Gene_ID,gns$gene_id)]))
colnames(aCML_rowData) = c('chr','geneSym','ensID','label','exonLength')

## 1. Prepare raw count matrix for DESeq2
# Raw count matrix
aCML_rawCnt = assays(aCML_sce)[['counts_raw']]

# Drop bullshit and Uninteresting genes
# Drop genes NOT on the main chromosomes
toc = aCML_rawCnt[rownames(aCML_rawCnt) %in% aCML_rowData$ensID[aCML_rowData$chr %in% paste0('chr',c(1:23,'X'))],]
dim(toc)
# Drop Non-expressed genes
w = which(rowSums(toc)>0)
toc = toc[w,]

## 2. Prepare colData for DESeq2
colDat = aCML_mdat[match(colnames(toc),aCML_mdat$sampleID),]
rownames(colDat) = colDat$sampleID
#colDat$setbp1_mut_inHotspot[is.na(colDat$setbp1_mut_inHotspot)] = 'none'

colDat$group = ifelse(grepl('wild_type',colDat$setbp1_mut),'wt',
                      ifelse(colDat$setbp1_mut == 'E858K','outside_hotspot','inside_hotspot'))

colDat$group = factor(colDat$group,levels = c('wt','inside_hotspot','outside_hotspot'))                            
#colDat$clin_age = ifelse(is.na(colDat$clin_age),'unknown',colDat$clin_age)


## 3. DEGs

aCML_res = bulkDESeq2(toc = toc, 
                      mDat = colDat[match(colnames(toc),colDat$sampleID),c('sampleID','group')],
                      cellTypeName = 'aCML',
                      formula = '~ %s',donorID = 'sampleID',groupID = 'group',
                      geneMap = aCML_rowData)


outDir = '~/lustre_mt22/SETBP1/Results/13_bulkDEGs_BEATaml/aug23'
saveRDS(aCML_res,file.path(outDir,'bulk_DESEQ2_aCML_Mut.vs.WT_output.RDS'))
aCML_res = readRDS(file.path(outDir,'bulk_DESEQ2_aCML_Mut.vs.WT_output.RDS'))

aCML_results = as.data.frame(aCML_res[['res']])
aCML_results = annotateGenes(aCML_results,geneMap = aCML_rowData)
table(aCML_results$contrast)
View(aCML_results[aCML_results$contrast == 'group_inside_hotspot_vs_wt' & 
                    aCML_results$geneSym %in% mod$geneSym[mod$direction == 'up_HCA'],])

aCML_deg = aCML_res[['all_de']]
table(aCML_deg$contrast)


write.csv(aCML_deg,file.path(outDir,'bulk_DESEQ2_aCML_Mut.vs.WT_deg.csv'))
aCML_deg = read.csv(file.path(outDir,'bulk_DESEQ2_aCML_Mut.vs.WT_deg.csv'))
## add direction to DEG
aCML_deg$direction = ifelse(aCML_deg$log2FoldChange >0,'down_WT','up_WT')
aCML_deg$list_id = paste0(aCML_deg$contrast,':',aCML_deg$direction)



res = aCML_deg[aCML_deg$contrast == 'group_inside_hotspot_vs_wt',]
table(res$padj < 0.05)
plot(res$log2FoldChange,-log10(res$padj),
     col = ifelse(res$padj< 0.05 & abs(res$log2FoldChange)>1,'red','black'),
     cex=0.5,
     xlab='logFC',
     ylab='-log10(adjusted p-value)',
     main='SETBP1.hs vs SETBP1.wt aCML'
)

table(res$log2FoldChange[res$padj< 0.05 & abs(res$log2FoldChange)>1] < 0)






## Plot heatmap of expression of these genes
aCML_deg.sub = aCML_deg[aCML_deg$contrast == 'group_inside_hotspot_vs_wt' & aCML_deg$padj < 0.01,]
genesSameDirection_aCML = c(mod$geneSym[mod$direction == 'up_HCA' & mod$geneSym %in% aCML_deg.sub$geneSym[aCML_deg.sub$log2FoldChange < 0 & aCML_deg.sub$padj < 0.01 ]],
                       mod$geneSym[mod$direction == 'down_HCA' & mod$geneSym %in% aCML_deg.sub$geneSym[aCML_deg.sub$log2FoldChange > 0 & aCML_deg.sub$padj < 0.01 ]])


tpm = aCML_exprData[rownames(aCML_exprData) %in% c(mod$ensID[mod$geneSym %in% genesSameDirection_aCML]),]
dim(tpm)
rownames(tpm) = mod$geneSym[match(rownames(tpm),mod$ensID)]
# remove genes with 0 everywhere
tpm = tpm[rowSums(tpm) > 0,]


# z-scaled by gene across samples
#tpm = as.data.frame(scale(t(tpm)))
mtx = tpm
group = aCML_mdat$setbp1_mut
group = group[match(colnames(mtx),aCML_mdat$sampleID)]
groupCols = col25[-6][1:n_distinct(group)]
names(groupCols) = unique(group)

hm = Heatmap(t((t(mtx))),show_row_dend = F,show_column_dend = F,
             show_row_names = T,show_column_names = F,
             cluster_rows = T,row_names_gp = gpar(fontsize=8,
                                                  col=ifelse(mod$isTF[match(rownames(mtx),mod$geneSym)],'red',
                                                             ifelse(mod$isCSM[match(rownames(mtx),mod$geneSym)],'blue','black'))),
             bottom_annotation = HeatmapAnnotation(group = group,
                                                   col = list(group = groupCols)),
             split=mod$direction[match(rownames(mtx),mod$geneSym)],
             column_split = group)

draw(hm)

rowOrder = row_order(draw(hm))



DotPlot(combSrat,group.by = 'dataset2',scale = T,
        #idents = unique(combSrat$ann[!grepl('unknown|Mast|MK|DC\\.precursor|doublets',combSrat$ann)]),
        features = rownames(mtx)[rev(c(rowOrder[[1]],rowOrder[[2]]))])+
  coord_flip() + theme(axis.text.x = element_text(size=12),
                       axis.text.y = element_text(size=10))