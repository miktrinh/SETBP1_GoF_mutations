### Prepare for deconvolution with the 1 bulk data from AML?? with SETBP1 mutation

### This script was first created in June 2022
### revisited in March 2023
#setwd('~/lustre_mt22/SETBP1/Results/June/CellSignalAnalysis/')


outDir = '~/lustre_mt22/SETBP1/Results/10_CSA/mar23_v1'
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}

setwd(outDir)

#############
# Libraries #
#############
# Load libraries
library(tidyverse)
library(Matrix)
library(Seurat)
library(biomaRt)
library(readxl)
library(SingleCellExperiment)
library(ComplexHeatmap)
library(patchwork)
source("/lustre/scratch117/casm/team274/mt22/generalScripts/utils/sc_utils.R")
source("/lustre/scratch117/casm/team274/mt22/generalScripts/utils/misc.R")

#### Useful functions ####
run_CSA = function(ref_srat,ref_dataset,scREF_prefix,geneWeight_fp = '~/lustre_mt22/Wilms/CellSignalAnalysis/geneWeights.tsv',geneMap=NULL){
  if(is.null(geneMap)){
    geneMap = read.csv('~/lustre_mt22/Aneuploidy/scGeneMap.csv')  
  }
  
  # Run CSA on the whole reference dataset
  #create matrix of raw counts
  data = ref_srat@assays$RNA@counts
  #convert geneSym to ensembl
  m = match(rownames(data),geneMap$geneSym)
  if(sum(is.na(m))>0){
    data = data[-which(is.na(m)),]
  }
  
  
  rownames(data) = geneMap$ensID[m]
  
  #write matrix
  #scREF_prefix = ref_dataset
  #scREF_prefix_sub = strsplit(scREF_prefix,split='_')
  #scREF_prefix_sub = scREF_prefix_sub[1]
  if(!file.exists(paste0("scData/",ref_dataset,".mtx"))){
    writeMM(obj = data, file = paste0("scData/",ref_dataset,".mtx"))
    #write rowNames file- keep prefix to all files the same- this is what you will put into command later
    write.table(rownames(data), file=paste0('scData/',ref_dataset,'_rowNames.tsv'), sep='\n', col.names = F,row.names = F, quote = FALSE)
    #writing corresponding column names (cell annotations)
    ref_srat@meta.data$cellID = rownames(ref_srat@meta.data)
    columns = ref_srat@meta.data[,c('cellID','annot')]
    
    columns$colNames = paste0(columns$annot,':',columns$cellID)
    write.table(columns$colNames, file=paste0('scData/',ref_dataset,'_columnNames.tsv'), sep='\n',col.names = F, row.names = F, quote = FALSE)
    
  }
  
  system(sprintf('mkdir -p %s',paste0(ref_dataset,'/',scREF_prefix)))
  cmd = paste0('python ~/lustre_mt22/Wilms/scripts/CellSignalAnalysis/cellSignalAnalysisV2_cellSignalAnalysis.py -b bulkData/bulk_samples.txt -s scData/',ref_dataset,
               ' -w ',geneWeight_fp,' -dg ',ref_dataset,'/',scREF_prefix,'/',scREF_prefix,'_out')
  message(sprintf('Running job on %s now',scREF_prefix))
  system(cmd)
  print(cmd)  
  
  message('Job finished!')
  
  #return(scREF_prefix)
}

#arm R with the functions - run this whole massive block
normaliseExposures = function(tgt){
  #Is this just the base?
  if(file.exists(paste0(tgt,'_fitExposures.tsv')))
    tgt = paste0(tgt,'_fitExposures.tsv')
  fit = read.table(tgt,sep='\t',header=TRUE,na.strings = T)
  
  #fit$celltype = rownames(fit$celltype)
  #rownames(fit) = fit$celltype
  #fit = fit[,-1]
  #Extract the goodness of fit rows
  gofNoms = c('pR2','fitCount','obsCount')
  gof = fit[gofNoms,]
  gof['log2(countRatio)',] = log2(unlist(gof['fitCount',]/gof['obsCount',]))
  #And the exposure table
  exposures = fit[-match(gofNoms,rownames(fit)),]
  #Normalise
  exposures = t(t(exposures)/colSums(exposures))
  #Done
  return(list(exposures=exposures,gof=gof,raw=fit[-match(gofNoms,rownames(fit)),]))
}


plot_CSA = function(fit,mdat){
  ### Plotting CSA heatmap #####
  exposureScale=c(0,0.5)
  #Colours for exposures
  hmCols = c('#ffffff','#f0f0f0','#d9d9d9','#bdbdbd','#969696','#737373','#525252','#252525','#000000')
  #Colours for pR2 metric
  pR2Cols  = c('#fff5eb','#fee6ce','#fdd0a2','#fdae6b','#fd8d3c','#f16913','#d94801','#a63603','#7f2704')
  #Colours for log library sizes
  libCols = c('#fcfbfd','#efedf5','#dadaeb','#bcbddc','#9e9ac8','#807dba','#6a51a3','#54278f','#3f007d')
  #And library size ratio
  libRatCols = c('#8c510a','#bf812d','#dfc27d','#f6e8c3','#f5f5f5','#c7eae5','#80cdc1','#35978f','#01665e')
  
  #Create the bottom annotation
  gof = fit$gof
  gof['fitCount',] = log10(gof['fitCount',])
  gof['obsCount',] = log10(gof['obsCount',])
  rownames(gof) = gsub('(.*Count)','log10(\\1)',rownames(gof))
  #gof = gof[,!grepl('^NA',colnames(gof))]
  #Decide on range for library size
  libRange = range(gof[grep('Count',rownames(gof)),])
  #Convert colours to ramps
  bCols = circlize::colorRamp2(seq(0,1,length.out=length(pR2Cols)),pR2Cols)
  lCols = circlize::colorRamp2(seq(libRange[1],libRange[2],length.out=length(libCols)),libCols)
  hmColObj = circlize::colorRamp2(seq(exposureScale[1],exposureScale[2],length.out=length(hmCols)),hmCols)
  lrCols = circlize::colorRamp2(seq(-1,1,length.out=length(libRatCols)),libRatCols)
  
  
  #####Add additional annotation
  # Reordering fitExposure sample order (columns) to match the order in bulk_samples (metadata) 
  # m = match(colnames(fit$exposures),mdat$Sample)
  # sum(is.na(m))
  # mdat = mdat[m[!is.na(m)],]
  
  group = mdat$bulk_cat
  supgroup = mdat$Subgroup
  
  
  #Colours for mutant-wildtype
  library(RColorBrewer)
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  
  n <- length(unique(mdat$bulk_cat))
  group_Cols = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))[1:n]
  names(group_Cols) = unique(mdat$bulk_cat)
  
  
  n <- length(unique(mdat$Subgroup))
  subgroup_Cols = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))[1:n]
  names(subgroup_Cols) = unique(mdat$Subgroup)
  
  #t = mdat %>% group_by(bulk_cat,Subgroup) %>% summarise(n=n())
  #t$col = rep(c('black','brown'),nrow(t)/2)[1:nrow(t)]
  #subgroup_Cols = rep(t$col,t$n)
  
  botAnno = HeatmapAnnotation(df = as.data.frame(t(gof)),
                              group = mdat$bulk_cat,
                              supgroup = mdat$Subgroup,
                              annotation_name_side = 'left',
                              col = list(pR2 =bCols,
                                         `log10(fitCount)` = lCols,
                                         `log10(obsCount)` = lCols,
                                         `log2(countRatio)` = lrCols,
                                         group = group_Cols))
  #supgroup = subgroup_Cols))
  
  # order fit$exposure rows by alphabetical order
  #fit$exposures = fit$exposures[,!grepl('^NA',colnames(fit$exposures))]
  row_names = rownames(fit$exposures)[rownames(fit$exposures) != 'Intercept']
  row_names = c('Intercept',row_names[order(row_names)])
  fit$exposures = fit$exposures[row_names,]
  
  
  ## Remove rows with no signal
  #row_sum = rowSums(fit$exposures)
  #fit$exposures = fit$exposures[which(row_sum >= 1),]
  #fit$exposures = fit$exposures[,order(colnames(fit$exposures),decreasing = F)]
  
  scREF_tag = sapply(strsplit(rownames(fit$exposures),split='_'),'[',1)
  levels_order = c('Intercept','SCPs','GOSH84','BJ9','HCA','Others','diploid')
  if('aBM' %in% unique(scREF_tag)){
    levels_order = c(levels_order,'aDiploid')
  }
  if('fBM' %in% unique(scREF_tag)){
    levels_order = c(levels_order,'fDiploid')
  }
  
  
  split_order = factor(scREF_tag,levels = c(levels_order,unique(scREF_tag[!scREF_tag %in% levels_order])))   
  
  
  hm = Heatmap(fit$exposures,
               col=hmColObj,
               name='Exposures',
               bottom_annotation=botAnno,
               split = split_order,
               #column_split = paste0(mdat$Subgroup[match(colnames(fit$exposures),mdat$Sample)],'_',mdat$Source[match(colnames(fit$exposures),mdat$Sample)]),
               column_split = mdat$bulk_cat,
               column_names_gp = gpar(fontsize=5,col=subgroup_Cols),
               column_title_gp = gpar(fontsize=5),column_title_rot = 90,
               row_names_gp = gpar(fontsize=5),
               show_row_names = T,
               cluster_rows = F, cluster_columns=F,
               column_names_max_height = unit(8, "cm"),
               row_gap = unit(6, "mm"),border_gp = gpar(col = "black", lty = 1),
               show_column_names = T,column_gap = unit(5, "mm"))
  
  
  return(hm)
}


# Global parameters ####
skipIfExists=T
geneMap = read.table('~/lustre_mt22/SETBP1/Data/SETBP1/cellranger612_count_42668_CG_SB_NB11528275_GRCh38-2020-A/filtered_feature_bc_matrix/features.tsv.gz')
colnames(geneMap) = c('ensID','geneSym','Feature','expr')
geneMap$geneSym[duplicated(geneMap$geneSym)] = paste0(geneMap$geneSym[duplicated(geneMap$geneSym)],'.1')




#####################
# Prepare bulk data #
#####################

for(dir in c('bulkData','scData','results')){
  if(!dir.exists(dir)){
    dir.create(dir,recursive = T)
  }
}

#bulk_sources = c('aCML','stJudes','TCGA_LAML','TARGET_AML','BEAT_AML')
bulk_sources = c('stJudes','BEAT_AML','JMML')


### 1. Generate one counts.txt file per bulkRNA samples #####
# Read in the main bulk TPM file - NEED RAW COUNT THOUGH!!!

bulk_samples = tibble()

if('aCML' %in% bulk_sources){
  aCML_sce_path = '~/lustre_mt22/SETBP1/Data/aCML_sce.RDS'
  if(!file.exists(aCML_sce_path)){
    bulk_aCML_dir = '/lustre/scratch125/casm/team274sb/mt22/SETBP1/Data/aCML/'
    if(length(list.files(bulk_aCML_dir)) == 1 & sum(grepl('.tar$',list.files(bulk_aCML_dir))) == 1){
      untar(list.files(bulk_aCML_dir), exdir = bulk_aCML_dir)
    }
    
    ### Write raw count
    sample_list = c() # to be added to the list of samples for CSA
    rawCnt = tibble() # Generate table of raw counts for all samples
    for(sample in list.files(bulk_aCML_dir,full.names = T)){
      if(grepl('GSE42146',sample)){
        next()
      }
      sampleID = gsub('.txt.gz$','',basename(sample))
      sample_list = c(sample_list,paste0(sampleID,'_counts.txt'))
      
      count = read.delim(sample,sep = '\t')
      w = grepl('Gene_ID|Sum_Exon_Length_Gene|Reads_Counts_accepted_hits_CMLPhn',colnames(count))
      tmp = count[,w]
      colnames(tmp) = c('ensID','geneLength',sampleID)
      
      if(ncol(rawCnt) == 0){
        rawCnt = tmp
      }else{
        rawCnt = merge(rawCnt,tmp,all=T,by=c('ensID','geneLength')) #cbind(rawCnt,data.frame(sampleID = raw_count[,3]))  #
      }
      
      if(!file.exists(paste0('bulkData/',sampleID,'_counts.txt'))){
        write.table(tmp,paste0('bulkData/',sampleID,'_counts.txt'),sep = '\t',quote = F,row.names = F)  
      }
    }
    
    #colnames(rawCnt) = gsub('-','.',gsub('_counts.txt$','',sample_list))
    # Prepare raw count and TPM count
    rawCnt = column_to_rownames(rawCnt,var = 'ensID')
    rawCnt = rawCnt[,-c(1)]
    
    ### Write TPM count
    if(!file.exists('~/lustre_mt22/SETBP1/Data/aCML_tpm.csv')){
      geneLen = as.vector(tmp[,2])
      
      rpk = apply(rawCnt, 2, function(x) x/(as.vector(geneLen)/1000))
      #normalize by the sample size using rpk values
      tpm = apply(rpk, 2, function(x) x / (sum(as.numeric(x)) / 1e6))  %>% as.data.frame()
      colNames = colnames(tpm)
      tpm$ensID = tmp$ensID
      tpm$geneLength = tmp$geneLength
      tpm = cbind(tpm[,c(ncol(tpm)-1,ncol(tpm))],tpm[,colNames])
      write.csv(tpm,'~/lustre_mt22/SETBP1/Data/aCML_tpm.csv')
    }else{
      tpm = read.delim('~/lustre_mt22/SETBP1/Data/aCML_tpm.csv',sep = ',')
      if('X' %in% colnames(tpm)){
        tpm=tpm[,-c(colnames(tpm) == 'X')]
      }
    }
    ## Extract SETBP1 expression level in TPM
    setbp1_aCML = tpm[tpm$ensID == 'ENSG00000152217',-c(1,2)] %>% t()
    
    
    ### Generate metadata
    # Import metadata file
    mdat = read.delim('~/lustre_mt22/SETBP1/Data/aCML/GSE42146_series_matrix.txt',skip = 33)
    mdat = mdat %>% as.matrix() %>% t() %>% as.data.frame()
    colnames(mdat) = gsub('!','',mdat[1,])
    mdat = mdat[-1,]
    rownames(mdat) = colnames(rawCnt)
    mdat$Sample_characteristics_ch1 = ifelse(grepl('wild type',mdat$Sample_characteristics_ch1),'wild_type',sapply(strsplit(mdat$Sample_characteristics_ch1,split = ' '),'[',3))
    
    m = match(colnames(rawCnt),rownames(mdat))
    
    aCML_mdat = data.frame(source = rep('aCML',nrow(mdat)),
                           sample_file = sample_list,
                           sampleID = colnames(rawCnt),
                           caseID = gsub('^.*_Gene_Expression_','',colnames(rawCnt)),
                           age = 'adult',
                           sex = '??',
                           cancerType = 'aCML',
                           cancerSubClass = 'aCML',
                           setbp1_mut = mdat$Sample_characteristics_ch1,
                           setbp1_mut_inHotspot = (mdat$Sample_characteristics_ch1 != 'wild_type'),
                           setbp1_tpm = setbp1_aCML[,1])
    
    
    # Generate a singleCellExperiment instance
    aCML_sce <- SingleCellExperiment(list(counts_raw=rawCnt,counts_tpm=tpm[,-c(1,2)]),
                                     colData=aCML_mdat,
                                     rowData=count[,c(1:5)],
                                     metadata=list(genome='hg18'))
    saveRDS(aCML_sce,aCML_sce_path)  
  }else{
    aCML_sce = readRDS(aCML_sce_path)
    aCML_mdat = as.data.frame(colData(aCML_sce))
  }
  
  bulk_samples = rbind(bulk_samples,aCML_mdat)
}






if('stJudes' %in% bulk_sources){
  stJudes_sce_path = '~/lustre_mt22/SETBP1/Data/StJudes_subset/stJudes_sce.RDS'
  if(!file.exists(stJudes_sce_path)){
    # the 1 relevant sample from StJude is ~/lustre_mt22/SETBP1/Data/Bulk_leukemia_StJudes/salmon_out/SJHM030201_D1.RNA-Seq.quant.genes.sf
    #stjude_aml = read.delim('~/lustre_mt22/SETBP1/Data/Bulk_leukemia_StJudes/salmon_out/SJHM030201_D1.RNA-Seq.quant.genes.sf')
    
    ## Read in pre-downloaded raw count (from Matt / Ellie)
    # This is for Gene ID quantification instead of transcript ID
    stjude_bulkrna = read.table('/lustre/scratch125/casm/team274sb/mt22/SETBP1/Data/Bulk_leukemia_StJudes/salmon_out/stJudes_bulk_txiConvGlen.csv',header = T)
    #stjude_bulkrna2 = read.delim('/lustre/scratch117/casm/team274/my4/oldScratch/Temp/stemnessData/StJudesTxCounts.tsv',sep = '\t',header = T)
    
    
    ### Generate metadata
    # List of samples with SETBP1 mutations, downloaded from stJude cloud website
    setbp1_mut_samples = read_excel('~/lustre_mt22/SETBP1/Data/StJudes_SETBP1_mutation_manifest.xlsx')
    # Import metadata
    mdat = read.delim('/lustre/scratch125/casm/team274sb/my4/oldScratch/Temp/stemnessData/StJudesBasicMetadata.tsv',sep = '\t',header = T)
    
    # randomly select 5 from each categories for negative controls: AML with/without SETBP1 mutations AND other diseases
    set.seed(1234)
    mdat = mdat %>% 
      dplyr::filter(sample_name %in% gsub('.RNA.Seq.quant.sf$','',colnames(stjude_bulkrna))) %>% 
      dplyr::filter(attr_sex != 'Not Available') %>% 
      dplyr::filter(grepl('RNA-Seq.bam$',file_path)) %>% 
      group_by(sj_diseases) %>% 
      mutate(idx = as.numeric(as.factor(sample_name)), 
             selected = idx %in% sample(1:max(idx),ifelse(unique(max(idx))>5,5,max(idx)))) %>% 
      dplyr::filter((sample_name %in% setbp1_mut_samples$Sample) | (selected & sj_diseases %in% c('AML','BALL','HM','NBL','TALL')))
    
    if(n_distinct(mdat$sample_name) != nrow(mdat)){
      message('Found duplicated sample names in St Judes data')
    }
    
    
    
    ### Generate raw count
    sample_list=c()
    rawCnt = tibble()
    # import list of StJudes SETBP1 mutation samples
    for(sampleID in c(mdat$sample_name)){
      sample_list = c(sample_list,paste0(sampleID,'_counts.txt'))
      
      if(!sum(grepl(sampleID,colnames(stjude_bulkrna)))){
        print('No sample data!')
      }
      tmp = stjude_bulkrna[,c("GeneName","geneLength",colnames(stjude_bulkrna)[grepl(sampleID,colnames(stjude_bulkrna))])]
      
      colnames(tmp) = c('ensID','geneLength',sampleID)
      
      if(ncol(rawCnt)==0){
        rawCnt = tmp
      }else{
        rawCnt = merge(rawCnt,tmp,all=T,by=c('ensID','geneLength'))
      }
      
      
      if(!file.exists(paste0('bulkData/',sampleID,'_counts.txt'))){
        write.table(tmp,paste0('bulkData/',sampleID,'_counts.txt'),sep = '\t',quote = F,row.names = F)  
      }
    }
    
    #colnames(rawCnt) = gsub('_counts.txt$','',sample_list)
    # Prepare raw count and TPM count
    rawCnt = column_to_rownames(rawCnt,var = 'ensID')
    rawCnt = rawCnt[,-c(1)]
    
    ### Write TPM count
    if(!file.exists('~/lustre_mt22/SETBP1/Data/StJudes_subset/stJudessub_tpm.csv')){
      geneLen = as.vector(tmp[,2])
      
      rpk = apply(rawCnt, 2, function(x) x/(as.vector(geneLen)/1000))
      #normalize by the sample size using rpk values
      tpm = apply(rpk, 2, function(x) x / (sum(as.numeric(x)) / 1e6))  %>% as.data.frame()
      colNames = colnames(tpm)
      tpm$ensID = tmp$ensID
      tpm$geneLength = tmp$geneLength
      tpm = cbind(tpm[,c(ncol(tpm)-1,ncol(tpm))],tpm[,colNames])
      write.csv(tpm,'~/lustre_mt22/SETBP1/Data/StJudes_subset/stJudessub_tpm.csv')
    }else{
      tpm = read.delim('~/lustre_mt22/SETBP1/Data/StJudes_subset/stJudessub_tpm.csv',sep = ',')
      if('X' %in% colnames(tpm)){
        tpm=tpm[,-c(colnames(tpm) == 'X')]
      }
    }
    
    
    ## Extract SETBP1 expression level in TPM
    setbp1_stJudes = tpm[tpm$ensID == 'ENSG00000152217',-c(1,2)] %>% t()
    
    
    ### Generate metadata
    m = match(colnames(rawCnt),mdat$sample_name)
    mdat = mdat[m,]
    
    stJude_mdat = data.frame(source = rep('StJude',nrow(mdat)),
                             sample_file = sample_list,
                             sampleID = mdat$sample_name,
                             caseID = mdat$sample_name,
                             age = ifelse(is.na(mdat$attr_age_at_diagnosis),'pediatric',mdat$attr_age_at_diagnosis),
                             sex = mdat$attr_sex,
                             cancerType = mdat$sj_diseases,
                             cancerSubClass = mdat$attr_diagnosis,
                             setbp1_mut = ifelse(is.na(setbp1_mut_samples$Mutation[match(mdat$sample_name,setbp1_mut_samples$Sample)]),'none',
                                                 setbp1_mut_samples$Mutation[match(mdat$sample_name,setbp1_mut_samples$Sample)]),
                             setbp1_mut_inHotspot = ifelse(mdat$sample_name %in% setbp1_mut_samples$Sample,'outside_hotspot','none'),
                             setbp1_tpm = setbp1_stJudes[,1])
    # I've manually checked the SETBP1 mutations, all of which are outside of the hotspot region. E858K is a common somatic mutation for SETBP1
    
    # Generate a singleCellExperiment instance
    stJude_sce <- SingleCellExperiment(list(counts_raw=rawCnt, counts_tpm=tpm[,-c(1,2)]),
                                       colData=stJude_mdat,
                                       rowData=stjude_bulkrna[,c("GeneName","geneLength")],
                                       metadata=list(sample_metadata=mdat,
                                                     setbp1_mut_sample_metadata=setbp1_mut_samples[setbp1_mut_samples$Sample %in% colnames(rawCnt),]))
    saveRDS(stJude_sce,stJudes_sce_path)  
  }else{
    stJude_sce = readRDS(stJudes_sce_path)
    stJude_mdat = as.data.frame(colData(stJude_sce))
  }
  
  
  bulk_samples = rbind(bulk_samples,stJude_mdat)
}





if('JMML' %in% bulk_sources){
  jmml_sce_path = '~/lustre_mt22/SETBP1/Data/JMML/JMML_sce.RDS'
  if(!file.exists(jmml_sce_path)){
    # Read in TPM bulk count
    tpmCnt = read.table('~/lustre_mt22/SETBP1/Data/JMML/GSE111709_JMML_BULK_tpm.txt')
    # Read in clinical data info
    jmml_clinical_mdat = read_excel('~/lustre_mt22/SETBP1/Data/JMML/JEM_20180853_TableS1.xlsx')
    # Read in sample metadata
    jmml_sample_mdat = read_delim('~/lustre_mt22/SETBP1/Data/JMML/GSE111895-GPL18573_series_matrix.txt',
                                  delim = "\t", escape_double = FALSE, col_names = FALSE, 
                                  trim_ws = TRUE, skip = 27)
    # re-orientate the dataset (transposing and add column names)
    jmml_sample_mdat = as.data.frame(t(as.matrix(jmml_sample_mdat)))
    colnames(jmml_sample_mdat) = gsub('^!','',jmml_sample_mdat[1,])
    jmml_sample_mdat = jmml_sample_mdat[-1,]
    jmml_sample_mdat$patientID = sapply(strsplit(jmml_sample_mdat$`Sample_title`,split='_'),'[',2)
    jmml_sample_mdat2 = merge(jmml_sample_mdat,jmml_clinical_mdat,by.x=c('patientID'),by.y=c('Patient ID'),all=T)
    
    ### Generate raw count
    sample_list=c()
    rawCnt = tibble()
    
    # JMML dataset was mapped to hg19_grch37 reference genome
    # This is the same as BEAT-AML
    geneAnn_jmml_fp = '~/lustre_mt22/SETBP1/Data/BEAT_AML/geneAnn_ens75_hg19_feb2014.csv'
    geneAnn = read.table(geneAnn_jmml_fp)
    
    
    ### Add cell type name to the sampleID
    for(i in 1:ncol(tpmCnt)){
      sampleID = colnames(tpmCnt)[i]
      sample_list = c(sample_list,paste0(sampleID,'_counts.txt'))
      
      data = data.frame(ensID = rownames(tpmCnt),
                        geneLength = geneAnn$length[match(rownames(tpmCnt),geneAnn$ensembl_gene_id)],
                        sample = tpmCnt[,i])
      colnames(data)[3] = sampleID
      
      ## Write count file
      if(!dir.exists('~/lustre_mt22/SETBP1/Results/June/CellSignalAnalysis/bulkData/JMML_Bulk')){
        dir.create('~/lustre_mt22/SETBP1/Results/June/CellSignalAnalysis/bulkData/JMML_Bulk')
      }
      if(!file.exists(paste0('~/lustre_mt22/SETBP1/Results/June/CellSignalAnalysis/bulkData/JMML_Bulk/',sampleID,'_TPMcounts.txt'))){
        write.table(data,paste0('~/lustre_mt22/SETBP1/Results/June/CellSignalAnalysis/bulkData/JMML_Bulk/',sampleID,'_TPMcounts.txt'),sep = '\t',quote = F,row.names = F)
      }
    }
    
    ### Write TPM count
    if(!file.exists('~/lustre_mt22/SETBP1/Data/JMML/JMML_tpm.csv')){
      tpmCnt$ensID = rownames(tpmCnt)
      tpmCnt$GeneName = geneAnn$external_gene_id[match(rownames(tpmCnt),geneAnn$ensembl_gene_id)]
      tpmCnt$geneLength = geneAnn$length[match(rownames(tpmCnt),geneAnn$ensembl_gene_id)]
      tpmCnt = tpmCnt[,c('ensID','GeneName','geneLength',colnames(tpmCnt)[!colnames(tpmCnt) %in% c('ensID','geneLength','GeneName')])]
      
      write.csv(tpmCnt,'~/lustre_mt22/SETBP1/Data/JMML/JMML_tpm.csv')
    }else{
      tpmCnt = read.delim('~/lustre_mt22/SETBP1/Data/JMML/JMML_tpm.csv',sep = ',')
      if('X' %in% colnames(tpmCnt)){
        tpmCnt=tpmCnt[,-c(colnames(tpmCnt) == 'X')]
      }
    }
    
    ## Extract SETBP1 expression level in TPM
    setbp1_jmml = tpmCnt[rownames(tpmCnt) == 'ENSG00000152217',-c(1:3)] %>% t()
    
    ### Generate metadata
    m = match(colnames(tpmCnt)[-c(1:3)],jmml_sample_mdat2$Sample_title)
    mdat = jmml_sample_mdat2[m,]
    
    jmml_mdat = data.frame(source = rep('JMML',nrow(mdat)),
                           sample_file = sample_list,
                           sampleID = mdat$Sample_title,
                           caseID = mdat$Sample_title,
                           age = ifelse(is.na(mdat$`Age at presentation`),'cord_blood',mdat$`Age at presentation`),
                           sex = ifelse(is.na(mdat$Gender),'unknown',mdat$Gender),
                           cancerType = 'JMML',
                           cancerSubClass = 'JMML',
                           setbp1_mut = ifelse(is.na(mdat$Mutation),'none',
                                               ifelse(grepl('SETBP1 C.2612 T>C',mdat$Mutation), 'p.I871S',
                                                      ifelse(grepl('SETBP1c.2608G>C',mdat$Mutation),'p.G870R','none'))),
                           setbp1_mut_inHotspot = ifelse(is.na(mdat$Mutation),'none',
                                                         ifelse(grepl('SETBP1',mdat$Mutation),'inside_hotspot','none')),
                           setbp1_tpm = setbp1_jmml[,1])
    
    # Generate a singleCellExperiment instance
    jmml_sce <- SingleCellExperiment(list(counts_tpm=tpmCnt[,-c(1:3)]),
                                     colData=jmml_mdat,
                                     rowData=tpmCnt[,c("ensID","GeneName","geneLength")],
                                     metadata=list(genome_version = 'GRCh37',
                                                   ensembl_version = '75_feb2014release',
                                                   sample_metadata=jmml_sample_mdat,
                                                   sample_clinical_metadata=jmml_clinical_mdat))
    saveRDS(jmml_sce,jmml_sce_path)  
  }else{
    jmml_sce = readRDS(jmml_sce_path)
    jmml_mdat = as.data.frame(colData(jmml_sce))
  }
  
  
  bulk_samples = rbind(bulk_samples,jmml_mdat)
}






if('BEAT_AML' %in% bulk_sources){
  #beat_sce_path = '~/lustre_mt22/SETBP1/Data/BEAT_AML/BEAT_AML_sce.RDS'
  beat_sce_path = 'BEAT_AML_sub_sce.RDS'
  if(!file.exists(beat_sce_path)){
    ### Generate metadata
    # Donor-level metadata
    beat_donor_mdat = read.csv('~/lustre_mt22/SETBP1/Data/BEAT_AML/All_sample_data.csv',skip = 3)
    # Keep donor/samples with RNAseq data only 
    beat_donor_mdat.sub = beat_donor_mdat[beat_donor_mdat$corresponding_rnaseq_sample != '',]
    # --> results in 648 unique donorID, of which only 134 are in the beat_sample_mdat.sub object below
    n_distinct(beat_donor_mdat.sub$PatientID)
    
    # Sample-level metadata
    beat_sample_mdat = read_excel('~/lustre_mt22/SETBP1/Data/BEAT_AML/beataml_waves1to4_sample_mapping.xlsx')
    ## check for number of overlapping patient ID
    table(unique(beat_sample_mdat$dbgap_subject_id) %in% unique(beat_donor_mdat$PatientID))
    unique(beat_sample_mdat$dbgap_subject_id)[!unique(beat_sample_mdat$dbgap_subject_id) %in% unique(beat_donor_mdat$PatientID)]
    #3 Patient dbGapID are not present in All_sample_data.csv file. Ignore these cases as they arent SETBP1 mutation samples / or no DNA info anyway
    # Keep donor/samples with RNAseq data only 
    beat_sample_mdat.sub = beat_sample_mdat[!is.na(beat_sample_mdat$dbgap_rnaseq_sample),]
    # --> results in 653 unique donorID (patientId)
    n_distinct(beat_sample_mdat.sub$dbgap_subject_id)
    
    # Merge donor and sample metadata
    beat_mdat = merge(beat_donor_mdat,beat_sample_mdat,by.x=c('PatientID','corresponding_rnaseq_sample'),by.y=c('dbgap_subject_id','dbgap_rnaseq_sample'))
    n_distinct(beat_mdat$PatientID) # 648 Patient ID
    
    ### BEAT_AML v1 includes samples from wave 1+2
    ### BEAT-AML v2 includes samples from wave 1+2, and additional samples from wave 3+4
    ### All relevant BEAT-AML v2 data (including their supplementary table S1) can be downloaded from github page: https://biodev.github.io/BeatAML2/
    ### Dataset includes Healthy control and Diseased samples. 
    ### Some patients have samples collected at multiple timepoints (eg. diagnosis, relapse, remission, residual etc.). However, not all timepoints have both RNA and DNA data
    table(beat_donor_mdat$DiseaseStageAtSpecimenCollection)
    
    ### In terms of variants analysis:
    ### there are two files where variants/mutations were reported:
    # 1. clinical summary: variants reported here are derived solely from clinical genotyping results (FAQ: https://biodev.github.io/BeatAML2/)
    # 2. WES/targeted sequencing output: variants reported after processing and filtering steps as detailed in table S6-Variant Filtering from the Tyner et al. 2018 paper
    ### This means that the variants for the same DNA sampleID may not match. Clinically, SETBP1 may have not been included in the gene panel. WES data analysis may have filtered out the SETBP1 mutations detected clinically
    ### Additionally, not all DNA samples have matching RNA samples
    
    ### Strategy for selecting samples with SETBP1 mutation:
    # 1. samples must have matching RNA and DNA data 
    # 2. the DNA data was reported as having a SETBP1 mutations in either the "clinical summary" file, or WES_mutation file, as well as evidenced from paper 1 Tyner et al. 2018 (if patients were from wave 1+2)
    # --> ie. taking the union of all the available information 
    # This strategy results in 12 RNA samples, with supporting DNA data with SETBP1 mutation, from 11 different patients.
    #   5 samples only have support from clinical data, not from WES_mutation. All of these 5 samples have SETBP1 mutations outside of the hotspot
    #   1 sample have support from BOTH clinical and WES_mutations
    #   6 samples only have support from WES_mutation data only
    
    
    
    
    
    #### WES Mutations data
    beat_mutation_mdat = read.delim('~/lustre_mt22/SETBP1/Data/BEAT_AML/beataml_wes_wv1to4_mutations_dbgap.txt')
    # Check for SETBP1 mutation info
    beat_mutation_mdat$setbp1_mut = ifelse(beat_mutation_mdat$symbol == 'SETBP1',beat_mutation_mdat$hgvsp_short,'none')
    beat_mutation_mdat$setbp1_mut_in_hotspot = ifelse((beat_mutation_mdat$symbol == 'SETBP1'),ifelse(beat_mutation_mdat$protein_position == '1070/1596','outside_hotspot','inside_hotspot'),'none')
    # Summarise mutation info by DNA sampleID
    beat_mutation_mdat.summary = beat_mutation_mdat %>% group_by(dbgap_sample_id) %>% summarise(n_mut = n_distinct(hgvsc),setbp1_mut_in_hotspot=paste0(unique(setbp1_mut_in_hotspot),collapse = ','),
                                                                                                setbp1_mut=paste0(unique(setbp1_mut),collapse = ','))
    beat_mutation_mdat.summary$setbp1_mut_in_hotspot = ifelse(grepl(',',beat_mutation_mdat.summary$setbp1_mut_in_hotspot),gsub('none,','',beat_mutation_mdat.summary$setbp1_mut_in_hotspot),as.character(beat_mutation_mdat.summary$setbp1_mut_in_hotspot))
    beat_mutation_mdat.summary$setbp1_mut = ifelse(grepl(',',beat_mutation_mdat.summary$setbp1_mut),gsub('none,','',beat_mutation_mdat.summary$setbp1_mut),as.character(beat_mutation_mdat.summary$setbp1_mut))
    
    # Keep only relevant DNA sampleID (these are samples with matching RNAsamples)
    beat_mutation_mdat.sub = beat_mutation_mdat[beat_mutation_mdat$dbgap_sample_id %in% beat_mdat$corresponding_dnaseq_sample,]
    
    
    #### Clinical Summary data
    clinical_mut = read_excel('~/lustre_mt22/SETBP1/Data/BEAT_AML/beataml_wv1to4_clinical.xlsx')
    n_distinct(clinical_mut$dbgap_subject_id) # 805 patients
    # two patientIDs are not in donor_mdat, these are the same 2 patientID which are present in sample_mdat but NOT donor_mdat
    table(clinical_mut$dbgap_subject_id[!clinical_mut$dbgap_subject_id %in% beat_donor_mdat$PatientID])
    
    # Check for SETBP1 mutation info
    clinical_mut$setbp1_mut = ifelse(grepl('SETBP1',clinical_mut$variantSummary),gsub("\\|.*$",'',str_extract(clinical_mut$variantSummary,pattern = 'SETBP1.*')),'none')
    clinical_mut$setbp1_mut_in_hotspot = ifelse((clinical_mut$setbp1_mut != 'none'),ifelse(grepl('D868N',clinical_mut$setbp1_mut),'inside_hotspot','outside_hotspot'),'none')
    clinical_mut.summary = clinical_mut %>% filter(!is.na(dbgap_dnaseq_sample)) %>% 
      group_by(dbgap_dnaseq_sample) %>% 
      summarise(setbp1_mut_in_hotspot=paste0(unique(setbp1_mut_in_hotspot),collapse = ','),
                setbp1_mut_long=paste0(unique(setbp1_mut),collapse = ','),
                setbp1_mut=ifelse(setbp1_mut_long == 'none','none',gsub(' .*$|;.*$|,.*$','',str_extract(setbp1_mut_long,'p\\..*'))))
    clinical_mut.summary$setbp1_mut[clinical_mut.summary$setbp1_mut_long == 'SETBP1 (H1100R; MAF 51%)'] = 'p.H1100R'
    
    # check that all 6 samples with SETBP1mut from clinical summary data are present in beat_mutation_mdat.summary
    table(clinical_mut.summary$dbgap_dnaseq_sample[clinical_mut.summary$setbp1_mut != 'none'] %in% beat_mutation_mdat.summary$dbgap_sample_id)
    
    # Update beat_mutation_mdat.summary
    clinical_mut.summary$wes_setbp1_mut = beat_mutation_mdat.summary$setbp1_mut[match(clinical_mut.summary$dbgap_dnaseq_sample,beat_mutation_mdat.summary$dbgap_sample_id)]
    w = clinical_mut.summary$dbgap_dnaseq_sample[(clinical_mut.summary$wes_setbp1_mut == 'none') & (clinical_mut.summary$setbp1_mut !='none')]
    
    beat_mutation_mdat.summary$setbp1_mut[match(w,beat_mutation_mdat.summary$dbgap_sample_id)] = clinical_mut.summary$setbp1_mut[match(w,clinical_mut.summary$dbgap_dnaseq_sample)]
    beat_mutation_mdat.summary$setbp1_mut_in_hotspot[match(w,beat_mutation_mdat.summary$dbgap_sample_id)] = clinical_mut.summary$setbp1_mut_in_hotspot[match(w,clinical_mut.summary$dbgap_dnaseq_sample)]
    
    
    # Add SETBP1 mutation info to main metadata
    beat_mdat$wes_n_mut = beat_mutation_mdat.summary$n_mut[match(beat_mdat$corresponding_dnaseq_sample,beat_mutation_mdat.summary$dbgap_sample_id)]
    beat_mdat$wes_clin_setbp1_mut = beat_mutation_mdat.summary$setbp1_mut[match(beat_mdat$corresponding_dnaseq_sample,beat_mutation_mdat.summary$dbgap_sample_id)]
    beat_mdat$wes_clin_setbp1_mut_in_hotspot = beat_mutation_mdat.summary$setbp1_mut_in_hotspot[match(beat_mdat$corresponding_dnaseq_sample,beat_mutation_mdat.summary$dbgap_sample_id)]
    
    ## Add additional information from clinical summary
    colnames(clinical_mut) = paste0('clin:',colnames(clinical_mut))
    beat_mdat_clin = merge(beat_mdat, clinical_mut,by.x=c('PatientID','corresponding_dnaseq_sample'),by.y=c('clin:dbgap_subject_id','clin:dbgap_dnaseq_sample'),all.x=T)
    
    
    
    ### Healthy control samples
    # no age information for healthy control samples...
    healthy_samples = beat_mdat_clin[beat_mdat_clin$rna_control != 'No',]
    
    
    
    # Select samples to keep: 
    # Keep samples with SETBP1 mutations 
    # AND samples without SETBP1 mutation but in the 'young' category (10yrs or younger)
    # I have checked that all young samples do not have SETBP1 mutations
    
    # Overall, there are 62 samples from BEAT-AML v2 selected
    set.seed(143)
    df1 = beat_mdat_clin[!is.na(beat_mdat_clin$wes_clin_setbp1_mut) & beat_mdat_clin$wes_clin_setbp1_mut != 'none',]
    df1$sampleCat = 'Disease_SETBP1'
    df2 = beat_mdat_clin %>% filter(!is.na(`clin:ageAtDiagnosis`) & beat_mdat_clin$AgeCategory == 'young' & beat_mdat_clin$`clin:ageAtDiagnosis` <=10)
    df2$sampleCat = 'Disease_noneSETBP1'
    
    healthy_samples$sampleCat = 'Healthy_control'
    
    beat_mdat.sub = rbind(df1,df2, 
                          healthy_samples)
    
    
    write.csv(beat_mdat.sub,'BEATAML.v2_cureated_selected_samples_mdat.csv')
    
    ### Generate raw count
    beat_raw_count = read.delim('~/lustre_mt22/SETBP1/Data/BEAT_AML/beataml_waves1to4_counts_dbgap.txt')
    # Subset raw count matrix to include only samples with SETBP1 mutations
    beat_raw_count = beat_raw_count[,colnames(beat_raw_count) %in% c(beat_mdat.sub$corresponding_rnaseq_sample,colnames(beat_raw_count)[1:4])]
    rownames(beat_raw_count) = beat_raw_count$stable_id
    
    
    ### Ensembl build 75 gene models on GRCh37
    geneAnn_BEAT_fp = '~/lustre_mt22/SETBP1/Data/BEAT_AML/geneAnn_ens75_hg19_feb2014.csv'
    if(file.exists(geneAnn_BEAT_fp)){
      geneAnn = read.table(geneAnn_BEAT_fp)
    }else{
      library(biomaRt)
      library(EDASeq)
      mart=useMart("ensembl",host="https://feb2014.archive.ensembl.org", dataset = "hsapiens_gene_ensembl")
      geneAnn = getBM(c('ensembl_gene_id','external_gene_id','chromosome_name'),mart = mart)
      geneAnn_bp_length = getGeneLengthAndGCContent(geneAnn$ensembl_gene_id, "hsa")
      geneAnn = merge(geneAnn,geneAnn_bp_length,by.x = 'ensembl_gene_id',by.y=0)
      geneAnn = geneAnn[!is.na(geneAnn$length),]
      write.table(geneAnn,'~/lustre_mt22/SETBP1/Data/BEAT_AML/geneAnn_ens75_hg19_feb2014.csv',col.names = T)  
    }
    
    # Write individual bulk count file for deconvolution
    beat_raw_count = merge(geneAnn,beat_raw_count,by.x = 'ensembl_gene_id',by.y=0)
    rowData = beat_raw_count[,c(1:8)]
    
    # save this table of raw count for BEAT-AML v2
    
    sample_list = c()
    rawCnt = tibble()
    for(sample in colnames(beat_raw_count)[-c(1:8)]){
      sampleID = paste0('BEAT_AML_',sample,'_counts.txt')
      sample_list = c(sample_list,sampleID)
      tmp = beat_raw_count[,c('ensembl_gene_id','length',sample)]
      colnames(tmp)[1:2]=c('ensID',	'geneLength')
      
      if(!file.exists(paste0('~/lustre_mt22/SETBP1/Results/June/CellSignalAnalysis/bulkData/','BEAT_AML_',sample,'_counts.txt'))){
        write.table(tmp,paste0('~/lustre_mt22/SETBP1/Results/June/CellSignalAnalysis/bulkData/','BEAT_AML_',sample,'_counts.txt'),sep = '\t',quote = F,row.names = F)  
      }
      
      if(ncol(rawCnt)==0){
        rawCnt = tmp
      }else{
        rawCnt = merge(rawCnt,tmp,by = c('ensID',	'geneLength'),all=T)
      }
      
    }
    
    # Prepare raw count and TPM count
    rawCnt = column_to_rownames(rawCnt,var = 'ensID')
    rawCnt = rawCnt[,-c(1)]
    
    #if(!file.exists('~/lustre_mt22/SETBP1/Data/BEAT_AML/BEAT_AML_sub_tpm.csv')){
    if(!file.exists('BEAT_AML_sub_tpm.csv')){
      geneLen = rowData$length
      
      rpk = apply(rawCnt, 2, function(x) x/(as.vector(geneLen)/1000))
      #normalize by the sample size using rpk values
      tpm = apply(rpk, 2, function(x) x / (sum(as.numeric(x)) / 1e6))  %>% as.data.frame()
      colNames = colnames(tpm)
      tpm$ensID = tmp$ensID
      tpm$geneLength = tmp$geneLength
      tpm = cbind(tpm[,c(ncol(tpm)-1,ncol(tpm))],tpm[,colNames])
      #write.csv(tpm,'~/lustre_mt22/SETBP1/Data/BEAT_AML/BEAT_AML_sub_tpm.csv')
      write.csv(tpm,'BEAT_AML_sub_tpm.csv')
    }else{
      #tpm = read.delim('~/lustre_mt22/SETBP1/Data/BEAT_AML/BEAT_AML_sub_tpm.csv',sep = ',')
      tpm = read.delim('BEAT_AML_sub_tpm.csv',sep = ',')
      if('X' %in% colnames(tpm)){
        tpm=tpm[,-c(colnames(tpm) == 'X')]
      }
    }
    
    # Extract SETBP1 expression level in TPM
    setbp1_beatAML = tpm[tpm$ensID == 'ENSG00000152217',-c(1,2)] %>% t()
    
    ### Generate metadata
    m = match(colnames(rawCnt),beat_mdat.sub$corresponding_rnaseq_sample)
    beat_mdat.sub = beat_mdat.sub[m,]
    
    beat_mdat = data.frame(source = rep('BEAT_AML',nrow(beat_mdat.sub)),
                           sample_file = sample_list,
                           sampleID = beat_mdat.sub$corresponding_rnaseq_sample,
                           caseID = beat_mdat.sub$PatientID,
                           age = beat_mdat.sub$AgeCategory,
                           clin_age = beat_mdat.sub$`clin:ageAtDiagnosis`,
                           sex = beat_mdat.sub$Gender,
                           cancerType = beat_mdat.sub$SpecificDxAtAcquisition,
                           cancerSubClass = beat_mdat.sub$FAB_BlastMorphology,
                           setbp1_mut = beat_mdat.sub$wes_clin_setbp1_mut,
                           setbp1_mut_inHotspot = beat_mdat.sub$wes_clin_setbp1_mut_in_hotspot,
                           setbp1_tpm = setbp1_beatAML[,1])
    
    # Generate a singleCellExperiment instance
    beat_sce <- SingleCellExperiment(list(counts_raw=rawCnt, counts_tpm=tpm[,-c(1,2)]),
                                     colData=beat_mdat,
                                     rowData=rowData,
                                     metadata=list(genome_version = 'ensemble_release_75_feb2014',
                                                   beat_donor_mdat = beat_donor_mdat,
                                                   beat_sample_mdat = beat_sample_mdat,
                                                   beat_mutation_mdat = beat_mutation_mdat))
    saveRDS(beat_sce,beat_sce_path)    
  }else{
    beat_sce = readRDS(beat_sce_path)
    beat_mdat = as.data.frame(colData(beat_sce))
  }
  beat_mdat2 = beat_mdat[,colnames(beat_mdat) !='age']
  colnames(beat_mdat2)[colnames(beat_mdat2) == 'clin_age'] = 'age'
  bulk_samples = rbind(bulk_samples,beat_mdat2)
}








########## Generate list of bulk samples to run CSA on #####
bulk_sources = c('fAdrenal_SCPs','stJudes','BEAT_AML','JMML')

#### Fix some spelling / Categorise age ##
bulk_samples$age[is.na(bulk_samples$age)] = 'unknown'
bulk_samples$age[bulk_samples$age == '14 months'] = 1.2
bulk_samples$age[grepl('yr',bulk_samples$age)] = as.numeric(gsub(' yr| yrs','',bulk_samples$age[grepl('yr',bulk_samples$age)]))

bulk_samples$ageCat = ifelse(!bulk_samples$age %in% c('young','older','oldest','',"adult","pediatric",'unknown','cord_blood'),
                             ifelse(as.numeric(bulk_samples$age) < 10,'<10',
                                    ifelse(as.numeric(bulk_samples$age) < 20,'<20',
                                           ifelse(as.numeric(bulk_samples$age) < 30, '<30',
                                                  ifelse(as.numeric(bulk_samples$age) < 50, '<40',
                                                         ifelse(as.numeric(bulk_samples$age) < 70,'<70','>70'))))),as.character(bulk_samples$age))


bulk_samples$sex = ifelse(bulk_samples$sex == 'male','Male',
                          ifelse(bulk_samples$sex == 'female','Female','unknown'))

bulk_samples$cancerSubClass = gsub('^m','M',bulk_samples$cancerSubClass)

bulk_samples$setbp1_mut = gsub('^p\\.','',bulk_samples$setbp1_mut)
bulk_samples$setbp1_mut[is.na(bulk_samples$setbp1_mut)] = 'none'

bulk_samples$setbp1_mut_inHotspot[bulk_samples$setbp1_mut == 'none'] = F
bulk_samples$setbp1_mut_inHotspot = ifelse(bulk_samples$setbp1_mut_inHotspot == 'none',F,
                                           ifelse(bulk_samples$setbp1_mut_inHotspot == '??','no_info',
                                                  ifelse(bulk_samples$setbp1_mut_inHotspot == 'inside_hotspot',T,
                                                         ifelse(bulk_samples$setbp1_mut_inHotspot == 'outside_hotspot',F,bulk_samples$setbp1_mut_inHotspot))))


# Add full path to bulk_samples$sample_file
bulk_samples$sample_file_fullpath = paste0('/lustre/scratch125/casm/team274sb/mt22/mt22/SETBP1/Results/June/CellSignalAnalysis/bulkData/',bulk_samples$sample_file)
bulk_samples$sample_file_fullpath[bulk_samples$source == 'JMML'] = gsub('_counts.txt','_TPMcounts.txt',bulk_samples$sample_file_fullpath[bulk_samples$source == 'JMML'])

## fAdr_SCPs ##

# read metadata
metadata = read.delim('~/lustre_mt22/Down_Leukemia/results/CellSignalAnalysis/bulkData/fAdrenal_SCPs/fAdrenal_SCPs_metadata.txt')
metadata$source = metadata$source
metadata$sampleID = metadata$SampleID
metadata$caseID = metadata$SampleID
metadata$age = 'unknown'
metadata$sex = 'unknown'
metadata$cancerType = metadata$SampleID
metadata$cancerSubClass = metadata$SampleID
metadata$setbp1_mut = 'none'
metadata$setbp1_mut_inHotspot = 'none'
metadata$setbp1_tpm = 1
metadata$ageCat = 'unknown'
metadata$sample_file_fullpath = '/lustre/scratch125/casm/team274sb/mt22/Down_Leukemia/results/CellSignalAnalysis/bulkData/fAdrenal_SCPs/fAdrenal_SCPs_counts.txt'

metadata = metadata[,colnames(bulk_samples)]
bulk_samples = rbind(bulk_samples,metadata)



write.table(bulk_samples$sample_file_fullpath,'bulkData/bulk_samples.txt',sep = '\n',quote = F,row.names = F,col.names = F)
write.table(bulk_samples,'bulkData/bulk_samples_metadata.txt',sep = '\t',quote = F,row.names = F,col.names = T)






################################
# Prepare single cell REF data #
################################
#### 2. Prepping single cell reference ####

### Generate SCPs fAdr ref
fAdr = readRDS('/lustre/scratch117/casm/team274/mt22/abnormal_karyotypes/adrenal/adrREF.rds') # 10X indexed GRCh38 1.2.0 reference
fAdr@meta.data$cellID = rownames(fAdr@meta.data)
cells_toKeep = fAdr@meta.data$cellID[fAdr@meta.data$cell_type == 'SCPs']
set.seed(1234)
cells_toKeep = cells_toKeep[-sample(1:length(cells_toKeep),200)]
fAdr = subset(fAdr,subset = cellID %in% cells_toKeep)
fAdr@meta.data$annot = fAdr@meta.data$cell_type
fAdr@meta.data$finalAnn = fAdr@meta.data$cell_type
fAdr@meta.data$Genotype = 'diploid'



#ref_datasets = c('setbp1_hca_G1','setbp1_G1')
ref_datasets = c('setbp1_G1')
for(ref_dataset in ref_datasets){
  message(sprintf('Generating sc reference from %s',ref_dataset))
  
  if(grepl('setbp1_hca|setbp1|hca',ref_dataset)){
    #combSrat = readRDS('~/lustre_mt22/SETBP1/Results/June/processed_ann_combSrat_SETBP1_ref24.RDS')
    combSrat = readRDS('~/lustre_mt22/SETBP1/Results/2_setbp1_HCA_annotation/jan23/SETBP1_REF24_mergedProcessed_annotated_noLowCnt.RDS')
    combSrat@meta.data$donorID[is.na(combSrat@meta.data$donorID)] = 'HCA'
  }
  
  if(!grepl('hca',ref_dataset)){
    # remove HCA
    combSrat = subset(combSrat,subset = batch %in% c('BJ9_1', 'GOSH84_1', 'GOSH84_2'))
  }else if(!grepl('setbp1',ref_dataset)){
    combSrat = subset(combSrat,subset = batch %in% c('HCA'))
  }
  
  ref_srat = combSrat
  
  
  
  
  
  
  if(grepl('aBM',ref_dataset)){
    aBM_fp = '~/lustre_mt22/Down_Leukemia/results/CellSignalAnalysis/scData/aBM_Laura.RDS'  
    if(!file.exists(aBM_fp)){
      aBM_mtx = readMM('/home/jovyan/lustre_mt22/REF_datasets/fetal_ref_tissues/AdultBM_NEW_RAW.mtx')
      aBM_mtx = t(as.matrix(aBM_mtx))
      aBM_mdat = read.csv('/home/jovyan/lustre_mt22/REF_datasets/fetal_ref_tissues/AdultBM_NEW_cellLabels.csv',header = F)
      aBM_cellID = read.csv('/home/jovyan/lustre_mt22/REF_datasets/fetal_ref_tissues/AdultBM_NEW_cellNames.csv',header = F)
      rownames(aBM_mdat) = aBM_cellID$V1
      colnames(aBM_mdat)[1] = 'cell.label'
      
      aBM_genes = read.csv('/home/jovyan/lustre_mt22/REF_datasets/fetal_ref_tissues/AdultBM_NEW_genes.csv',header = F)
      colnames(aBM_genes)[1] = 'geneSym'
      rownames(aBM_mtx) = aBM_genes$geneSym
      colnames(aBM_mtx) = aBM_cellID$V1
      aBM = CreateSeuratObject(counts = aBM_mtx,meta.data = aBM_mdat)
      rm(aBM_mtx)
      aBM@misc$geneMap = aBM_genes
      aBM <- NormalizeData(aBM)
      aBM <- FindVariableFeatures(aBM, selection.method = "vst")
      aBM <- ScaleData(aBM, features = rownames(aBM))
      # CellCycle Scoring
      s.genes <- cc.genes$s.genes
      g2m.genes <- cc.genes$g2m.genes
      aBM <- CellCycleScoring(aBM, s.features = s.genes, g2m.features = g2m.genes)
      
      aBM@meta.data$Genotype = 'diploid'
      #aBM@meta.data$Sex = aBM@meta.data$sex
      aBM@meta.data$annot = as.character(aBM@meta.data$cell.label)
      aBM@meta.data$finalAnn = aBM@meta.data$annot
      
      saveRDS(aBM,aBM_fp)  
    }else{
      aBM = readRDS(aBM_fp)
      aBM@meta.data$batch = 'aBM'
      aBM@meta.data$Sex = aBM@meta.data$sex
      
    }
    ref_srat = merge_seurat_objects(ref_srat,aBM,keepAllGenes = F,genomeVersions = c('v38','v38'))
  }
  
  if(grepl('norm_fBM',ref_dataset)){
    #fBM_fp = '~/lustre_mt22/Down_Leukemia/results/CellSignalAnalysis/scData/fBM_Laura.RDS'
    fBM_fp = '~/lustre_mt22/REF_datasets/fBM_Jardine21/fBM_diploid_Laura_sratObj.RDS'  
    if(!file.exists(fBM_fp)){
      # #fBM_mtx = readMM('/home/jovyan/lustre_mt22/REF_datasets/fig1b_fbm_scaled_gex_updated_dr_20210104.mtx')
      # #fBM_mtx = t(as.matrix(fBM_mtx))
      # fBM_mtx = readMM('/home/jovyan/lustre_mt22/REF_datasets/fetal_ref_tissues/FetalBM_split_SoupXD_BroadLab400wSCP_raw.mtx')
      # fBM_mtx = as.matrix(fBM_mtx)
      # 
      # #fBM_mdat = read.csv('/home/jovyan/lustre_mt22/REF_datasets/fig1b_fbm_scaled_gex_updated_dr_20210104_OBS.csv')
      # fBM_cellname = read.csv('/home/jovyan/lustre_mt22/REF_datasets/fetal_ref_tissues/FetalBM_split_SoupXD_BroadLab400wSCP_raw_cellNames.csv',header = F)
      # fBM_cellLab = read.csv('/home/jovyan/lustre_mt22/REF_datasets/fetal_ref_tissues/FetalBM_split_SoupXD_BroadLab400wSCP_raw_cellLabels.csv',header = F)
      # 
      # fBM_mdat = data.frame(cellID = fBM_cellname$V1,cell.labels = fBM_cellLab$V1)
      # rownames(fBM_mdat) = fBM_mdat$cellID
      # 
      # #fBM_genes = read.csv('/home/jovyan/lustre_mt22/REF_datasets/fig1b_fbm_scaled_gex_updated_dr_20210104_VAR.csv')
      # #colnames(fBM_genes)[1] = 'geneSym'
      # #rownames(fBM_mtx) = fBM_genes$geneSym
      # #colnames(fBM_mtx) = fBM_mdat$index
      # fBM_genes = read.csv('/home/jovyan/lustre_mt22/REF_datasets/fetal_ref_tissues/FetalBM_split_SoupXD_BroadLab400wSCP_raw_rowNames.tsv',header = F)
      # fBM_geneSym = read.csv('/home/jovyan/lustre_mt22/REF_datasets/fetal_ref_tissues/FetalBM_split_SoupXD_BroadLab400wSCP_raw_genes.csv',header = F)
      # fBM_genes = data.frame(geneSym = fBM_geneSym$V1,ensID = fBM_genes$V1)
      # 
      # rownames(fBM_mtx) = fBM_geneSym$V1
      # colnames(fBM_mtx) = fBM_mdat$cellID
      # 
      # # Remove SCP
      # fBM_mtx = fBM_mtx[,!colnames(fBM_mtx) %in% fBM_mdat$cellID[fBM_mdat$cell.labels == 'SCP']]
      # fBM_mdat = fBM_mdat[fBM_mdat$cell.labels != 'SCP',]
      # fBM = CreateSeuratObject(counts = fBM_mtx,meta.data = fBM_mdat)
      # rm(fBM_mtx)
      # fBM@misc$geneMap = fBM_genes
      # fBM <- NormalizeData(fBM)
      # fBM <- FindVariableFeatures(fBM, selection.method = "vst")
      # fBM <- ScaleData(fBM, features = rownames(fBM))
      # # CellCycle Scoring
      # s.genes <- cc.genes$s.genes
      # g2m.genes <- cc.genes$g2m.genes
      # fBM <- CellCycleScoring(fBM, s.features = s.genes, g2m.features = g2m.genes)
      # 
      # fBM@meta.data$Genotype = 'diploid'
      # #fBM@meta.data$Sex = fBM@meta.data$sex
      # fBM@meta.data$annot = as.character(fBM@meta.data$cell.labels)
      # fBM@meta.data$finalAnn = fBM@meta.data$annot
      # 
      # saveRDS(fBM,fBM_fp)  
      message('fBM sratObj does not exist! Please check and run ~/lustre_mt22/generalScripts/regenerate_fBM_sratObj.R if required')
      #system('Rscripts ~/lustre_mt22/generalScripts/regenerate_fBM_sratObj.R')
    }else{
      fBM = readRDS(fBM_fp)
      fBM@meta.data$batch = 'fBM'
      fBM@meta.data$Sex = fBM@meta.data$sex
    }
    
    ref_srat = merge_seurat_objects(ref_srat,fBM,keepAllGenes = F,genomeVersions = c('v38','v38'))
  }
  
  
  
  
  
  # Remove cycling cells if needed
  if(grepl('_G1',ref_dataset)){
    if(n_distinct(ref_srat$Phase) != 1){
      ref_srat = subset(ref_srat,subset = Phase == 'G1')
    }
  }else{
    # Make sure same cell types have the same labels
    # if(grepl('fLiver',ref_dataset)){
    #   ref_srat@meta.data$annot2 = ref_srat@meta.data$annot
    #   ref_srat@meta.data$annot2[ref_srat@meta.data$annot == 'Early.lymphoid_T.lymphocyte'] = 'T_cell_precursor'
    #   ref_srat@meta.data$annot2[ref_srat@meta.data$annot == 'pDC.precursor'] = 'pDC'
    #   ref_srat@meta.data$annot2[ref_srat@meta.data$annot == 'pre_B_cell'] = 'pre.B.cell'
    #   ref_srat@meta.data$annot2[ref_srat@meta.data$annot == 'pro_B_cell'] = 'pro.B.cell'
    #   ref_srat@meta.data$annot2[ref_srat@meta.data$annot == 'Mid_erythroid'] = 'Mid.Erythroid'
    #   ref_srat@meta.data$annot2[ref_srat@meta.data$annot == 'Naive_B_cell'] = 'B.cell'
    #   
    #   ref_srat@meta.data$annot2[ref_srat@meta.data$annot %in% c('Non_classical_Monocyte','Classical_Monocyte')] = 'Monocyte'
    #   
    # }
    
    
    
    
    # subsampling cycling cells such that there's the same %of cycling cells per cell type per genotype
    nCyclingCells = ref_srat@meta.data %>% group_by(Genotype,annot2) %>% summarise(n_cells = n(),
                                                                                   n_G1 = sum(Phase == 'G1'),
                                                                                   n_S = sum(Phase == 'S'),
                                                                                   n_G2M = sum(Phase == 'G2M'),
                                                                                   perc_S = 100*n_S/n_cells,
                                                                                   perc_G2M = 100*n_G2M/n_cells)
    # Removing categories with fewer than 30 cells
    nCyclingCells = nCyclingCells[nCyclingCells$n_cells >= 30,]
    
    # For each cell type, down sample G1 G2M and S cells down to match the proportion of the lowest group
    
    nCyclingCells = nCyclingCells %>% group_by(annot2) %>% mutate(n_cell_toKeep = min(n_cells),
                                                                  perc_S_tokeep = perc_S[n_cells == min(n_cells)],
                                                                  perc_G2M_tokeep = perc_G2M[n_cells == min(n_cells)],
                                                                  n_S_tokeep = round(n_cell_toKeep*perc_S_tokeep/100),
                                                                  n_G2M_tokeep = round(n_cell_toKeep*perc_G2M_tokeep/100),
                                                                  n_G1_tokeep = n_cell_toKeep - n_S_tokeep - n_G2M_tokeep)
    # some odd case need to be fized, however these are not very important cell types 
    nCyclingCells$n_G1_tokeep[nCyclingCells$n_G1_tokeep > nCyclingCells$n_G1] = nCyclingCells$n_G1[nCyclingCells$n_G1_tokeep > nCyclingCells$n_G1]
    nCyclingCells$n_G2M_tokeep[nCyclingCells$n_G2M_tokeep > nCyclingCells$n_G2M] = nCyclingCells$n_G2M[nCyclingCells$n_G2M_tokeep > nCyclingCells$n_G2M]
    nCyclingCells$n_S_tokeep[nCyclingCells$n_S_tokeep > nCyclingCells$n_S] = nCyclingCells$n_S[nCyclingCells$n_S_tokeep > nCyclingCells$n_S]
    #ref_srat@meta.data$n_G2M_tokeep = nCyclingCells$n_G2M_tokeep[match(paste0(ref_srat@meta.data$Genotype,'_',ref_srat@meta.data$annot2),paste0(nCyclingCells$Genotype,'_',nCyclingCells$annot2))]
    #ref_srat@meta.data$n_G1_tokeep = nCyclingCells$n_G1_tokeep[match(paste0(ref_srat@meta.data$Genotype,'_',ref_srat@meta.data$annot2),paste0(nCyclingCells$Genotype,'_',nCyclingCells$annot2))]
    #ref_srat@meta.data$n_S_tokeep = nCyclingCells$n_S_tokeep[match(paste0(ref_srat@meta.data$Genotype,'_',ref_srat@meta.data$annot2),paste0(nCyclingCells$Genotype,'_',nCyclingCells$annot2))]
    ref_srat@meta.data$n_cell_toKeep[ref_srat@meta.data$Phase == 'G1'] = nCyclingCells$n_G1_tokeep[match(paste0(ref_srat@meta.data$Genotype[ref_srat@meta.data$Phase == 'G1'],'_',ref_srat@meta.data$annot2[ref_srat@meta.data$Phase == 'G1']),paste0(nCyclingCells$Genotype,'_',nCyclingCells$annot2))]
    ref_srat@meta.data$n_cell_toKeep[ref_srat@meta.data$Phase == 'G2M'] = nCyclingCells$n_G2M_tokeep[match(paste0(ref_srat@meta.data$Genotype[ref_srat@meta.data$Phase == 'G2M'],'_',ref_srat@meta.data$annot2[ref_srat@meta.data$Phase == 'G2M']),paste0(nCyclingCells$Genotype,'_',nCyclingCells$annot2))]
    ref_srat@meta.data$n_cell_toKeep[ref_srat@meta.data$Phase == 'S'] = nCyclingCells$n_S_tokeep[match(paste0(ref_srat@meta.data$Genotype[ref_srat@meta.data$Phase == 'S'],'_',ref_srat@meta.data$annot2[ref_srat@meta.data$Phase == 'S']),paste0(nCyclingCells$Genotype,'_',nCyclingCells$annot2))]
    ref_srat@meta.data$n_cell_toKeep[is.na(ref_srat@meta.data$n_cell_toKeep)] = 0
    # Do the sub-sampling
    set.seed(1234)
    cellID_toKeep = ref_srat@meta.data %>% group_by(Genotype,annot2,Phase) %>% mutate(id = seq(1:n()),
                                                                                      toKeep = ifelse(id %in% sample(c(1:max(id)),unique(n_cell_toKeep)),T,F))
    cellID_toKeep = cellID_toKeep$cellID[cellID_toKeep$toKeep == T]
    ref_srat = subset(ref_srat,subset = cellID %in% cellID_toKeep)
  }
  
  
  ### Add cell number 
  ref_srat@meta.data$cellID = rownames(ref_srat@meta.data)
  ref_srat@meta.data$annot = paste0(ref_srat@meta.data$donorID,'_',ref_srat@meta.data$finalAnn_v2)
  nCell = table(ref_srat@meta.data$annot)
  ref_srat@meta.data$nCell = nCell[match(ref_srat@meta.data$annot,names(nCell))]
  ref_srat@meta.data$annot = paste0(ref_srat@meta.data$annot,' (n=',ref_srat@meta.data$nCell,')')
  
  
  #if(!grepl('BM',ref_dataset)){
  # Im suspecting that there might be issue with different genome versions etc --> investigate this further?
  ### Merge ref_srat with fAdr_SCPs
  ref_srat = merge_seurat_objects(ref_srat,fAdr,keepAllGenes = F,genomeVersions = c('v38','v38'))    
  #}
  
  
  
  
  
  
  
  
  
  
  
  message(sprintf('Running CSA using ref_dataset: %s and bulk dataset: %s',ref_dataset,paste(bulk_sources,collapse = '_')))
  
  
  
  
  
  #### Save scREF data in the requried format for CSA #####
  #out_csa = data.frame()
  
  fit_exp = list()
  for(i in c(1)){
    print(i)
    
    scREF_prefix = paste0(ref_dataset,'_',paste(bulk_sources,collapse = '.'),'_',i)
    
    if(i==1){
      ref_srat.sub = ref_srat
    }else{
      #### Sub-sampling the reference dataset for 5-fold cross validation
      ### Sub-sampling 80% of the reference dataset (by cell type)
      cellID_toKeep = ref_srat@meta.data %>% group_by(annot) %>% 
        mutate(id = seq(1,n()),selected = ifelse(id %in% sample(1:max(id),round(0.8*n())),T,F))
      cellID_toKeep = cellID_toKeep$cellID[cellID_toKeep$selected]
      
      ref_srat.sub = subset(ref_srat,subset = cellID %in% cellID_toKeep)
      
    }
    # Run CSA on ref_srat
    
    run_CSA(ref_srat,ref_dataset,scREF_prefix)  
    
    
    
    
    
    
    ### Redefine mdat
    mdat = bulk_samples
    
    # # Read in fit value
    # # For some weird reasons, there are 'NA...' columns in the results. However, I've checked that these are repeated results from some other columns --> can just be removed!
    fit <- normaliseExposures(tgt = paste0(ref_dataset,'/',scREF_prefix,"/",scREF_prefix,"_out"))
    samples_toKeep = intersect(colnames(fit$exposures),mdat$sampleID)
    mdat = mdat[mdat$sampleID %in% samples_toKeep,]
    fit = lapply(fit, function(x){x=x[,colnames(x) %in% samples_toKeep]})
    all((nrow(mdat) == ncol(fit$exposures)) == (ncol(fit$gof) == ncol(fit$raw)))
    
    
    # bulk_cat = source_cancer/healthy
    # subgroup = setbp1_in/out/none
    mdat$cancerType_broad = mdat$cancerType
    mdat$cancerType_broad[grepl('Healthy',mdat$cancerType_broad)] = 'Healthy'
    mdat$cancerType_broad[grepl('Down syndrome',mdat$cancerType_broad)] = 'MLDS'
    mdat$cancerType_broad[grepl('AML|Acute myeloid leukaemia|monocytic',mdat$cancerType_broad)] = 'AML'
    mdat$cancerType_broad[grepl('Mixed phenotype',mdat$cancerType_broad)] = 'MixedLineage'
    mdat$cancerType_broad[grepl('HGG|NBL|OS|ST|Therapy',mdat$cancerType_broad)] = 'Other_cancer'
    mdat$cancerType_broad[grepl('cord_blood',mdat$age)] = 'Healthy'
    
    
    mdat$bulk_cat = paste0(mdat$source,':',mdat$cancerType_broad)
    mdat$bulk_cat[grepl('fetal_adr_scRNA:fAdrenal_SCPs',mdat$bulk_cat)] = 'fAdr_SCPs'
    
    mdat$Subgroup = paste0(ifelse(mdat$setbp1_mut == 'none','none','SETBP1mut'),':',
                           ifelse(mdat$setbp1_mut_inHotspot == TRUE,'in',ifelse(mdat$setbp1_mut == 'none','none','out')))
    
    mdat$bulk_cat = factor(mdat$bulk_cat,levels = c('BEAT_AML:Healthy',"JMML:Healthy",
                                                    "JMML:JMML",'StJude:AML','BEAT_AML:AML','BEAT_AML:MixedLineage','BEAT_AML:MLDS', 'BEAT_AML:Other_cancer',
                                                    'StJude:BALL','StJude:ETV','StJude:HM','StJude:Other_cancer','StJude:TALL','fAdr_SCPs'))
    mdat = arrange(mdat, as.numeric(bulk_cat), sampleID)
    
    ### re-order fit$exposure and gof / raw to match the order of sample names in mdat
    fit = lapply(fit, function(x){x=x[,match(mdat$sampleID,colnames(x))]})
    # Check that colnames are in correct order
    sum(sapply(fit, function(x){all(colnames(x) == mdat$Sample)})) == 3
    
    hm = plot_CSA(fit,mdat)
    #pdf(paste0(ref_dataset,'/',scREF_prefix,'/',scREF_prefix,'_hm.pdf'),width = 305,height = 50)
    pdf(paste0(ref_dataset,'/',scREF_prefix,'/',scREF_prefix,'_hm.pdf'),width = 35,height = 20)  
    draw(hm)
    dev.off()
  }
}