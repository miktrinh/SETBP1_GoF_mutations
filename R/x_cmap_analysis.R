## Process SETBP1 Connectivity Map results
c('shared','somatic_all','lymphoid')
signature = 'shared'

library(tidyverse)
library(ggbeeswarm)
source('~/lustre_mt22/generalScripts/utils/misc.R')

if(signature == 'shared'){
  outDir = 'Results/06_cmap/'  
  gct_fp = file.path(outDir,'shared_my_analysis.sig_queryl1k_tool.67b339cd189d0e00137f9083//arfs/TAG/query_result.gct')
  gct_fp = file.path(outDir,'setbp1_myeloid_shared_my_analysis.sig_queryl1k_tool.67c1991c6ccd460013cf0d9e/arfs/TAG/query_result.gct')
  gct_fp = file.path(outDir,'setbp1_myeloid_shared_250325_my_analysis.sig_queryl1k_tool.67e2da0fed36750013648117/arfs/TAG/query_result.gct')
}else if(signature == 'somatic_all'){
  outDir = '~/lustre_mt22/SETBP1/Results/cmap/' 
  gct_fp = file.path(outDir,'somatic.all_my_analysis.sig_queryl1k_tool.67b33a48189d0e00137f9085/arfs/TAG/query_result.gct')
}else if(signature == 'lymphoid'){
  outDir = '~/lustre_mt22/SETBP1/Results/cmap/' 
  gct_fp = file.path(outDir,'lymphoid.unique_my_analysis.sig_queryl1k_tool.67b33ab6189d0e00137f9087/arfs/TAG/query_result.gct')
}


if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}

setwd(outDir)

top_n = 10
cellline = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/CLUE/cmap_cellline.csv')
cellline$group = ifelse(cellline$Culture.Properties == 'Suspension','Leukaemia',cellline$Type)
cellline = rbind(cellline,
                 data.frame(Cell.Line = c('MCF10A','HEK293','JURKAT','NEU','PHH','HUVEC','SKB','HPTEC'),
                            Culture.Properties = '?',
                            Description = '?',
                            Type = 'non-cancer',
                            group = 'non-cancer'))

cellline = rbind(cellline,
                 data.frame(Cell.Line = c('HUH7','NPC','MDAMB231','HELA','HCT116','U2OS','SKBR3','SW480','H1299','HS578T','NCIH508','SHSY5Y','NCIH2073','NCIH596','BT20'),
                            Culture.Properties = '?',
                            Description = '?',
                            Type = 'Cancer',
                            group = 'Cancer'))

cellline = rbind(cellline,
                 data.frame(Cell.Line = c('U937','HAP1','ASC','SKL','HL60','NOMO1'),
                            Culture.Properties = '?',
                            Description = '?',
                            Type = 'Leukaemia',
                            group = 'Leukaemia'))


cellline$group[cellline$group == 'Cancer'] = 'Solid cancer'










## Import cmap results
dat.gct <- read.delim(gct_fp,skip = 2)

## Extract compound only
cp = dat.gct[dat.gct$pert_type == 'trt_cp',]
cp$norm_cs = as.numeric(cp$norm_cs)

## extract top 20 similar and dissimilar entries
cp = cp[cp$moa != '-666',]
cp_toKeep = cp[cp$moa != '-666',]
cp_toKeep$norm_cs_ranking = order(cp_toKeep$norm_cs)
cp_toKeep_dissimilar = cp_toKeep[cp_toKeep$norm_cs_ranking %in% c(1:top_n),]
cp_toKeep_dissimilar$group = 'dissimilar'
cp_toKeep_similar = cp_toKeep[cp_toKeep$norm_cs_ranking %in% c(max(cp_toKeep$norm_cs_ranking):(max(cp_toKeep$norm_cs_ranking)-top_n)),]
cp_toKeep_similar$group = 'similar'
cp_toKeep_nfkb = cp_toKeep[grepl('nfkb|NFkB|NF-?κB|NFKB',cp_toKeep$moa),]
cp_toKeep_nfkb$group = ifelse(grepl('NFKB activator',cp_toKeep_nfkb$moa),'nfkb_activator',
                              'nfkb_inhibitor')

cp_toKeep = do.call(rbind,list(cp_toKeep_similar,cp_toKeep_dissimilar,cp_toKeep_nfkb))
cp_toKeep$group = ifelse(grepl('nfkb|NFkB|NF-?κB|NFKB',cp_toKeep$moa) &
                           grepl('NFKB activator',cp_toKeep$moa) &
                           cp_toKeep$group %in% c('dissimilar','similar') ,
                         pasteo('nfkb_activator_',cp_toKeep$group),
                         ifelse(grepl('nfkb|NFkB|NF-?κB|NFKB|NFKB',cp_toKeep$moa) & 
                                  grepl('NFKB inhibitor',cp_toKeep$moa) &
                                  cp_toKeep$group %in% c('dissimilar','similar'),
                         paste0('nfkb_inhibitor_',cp_toKeep$group),cp_toKeep$group))
table(cp_toKeep$group)
cp.sub = cp[cp$pert_iname %in% cp_toKeep$pert_iname,]
mtx = pivot_wider(cp.sub[,c('pert_iname','cell_iname','norm_cs')],id_cols = 'pert_iname',values_from = 'norm_cs', names_from = 'cell_iname')
mtx = column_to_rownames(mtx,'pert_iname')

for(i in 1:nrow(mtx)){
  for(j in 1:ncol(mtx)){
    mtx[i,j] = paste(unlist(mtx[i,j]),collapse = ':')
    if(!grepl(':',mtx[i,j])){
      mtx[i,j] = as.numeric(mtx[i,j])
    }else{
      # both -ve or +ve
      values = as.numeric(strsplit(mtx[i,j][[1]],split = ':')[[1]])
      if(all(values > 0)){
        mtx[i,j] = min(values)
      }else if(all(values < 0)){
        mtx[i,j] = max(values)
      }else{
        mtx[i,j] = NA
      }
    }
  }
}

#mtx2 = do.call(cbind,apply(mtx,2,unlist))
mtx2 = apply(mtx,2,unlist)
mtx2 = apply(mtx2, 2, as.numeric)
rownames(mtx2) = rownames(mtx)

library(ComplexHeatmap)
mtx2 = mtx2[match((cp_toKeep$pert_iname[order(cp_toKeep$norm_cs_ranking,decreasing = F)]),rownames(mtx2)),]
# Reorder to group MOA together
moa_dissimilar = unique(cp_toKeep_dissimilar$moa)
moa_similar = unique(cp_toKeep$moa)
cp_toKeep$moa = factor(cp_toKeep$moa,unique(c(moa_dissimilar,moa_similar)))

mtx2 = mtx2[order(cp_toKeep$moa[match(rownames(mtx2),cp_toKeep$pert_iname)]),]
rownames(mtx2) = paste0(rownames(mtx2),' (',cp_toKeep$moa[match(rownames(mtx2),cp_toKeep$pert_iname)],')')
dim(mtx2[!duplicated(mtx2),])
mtx2 = mtx2[!duplicated(mtx2),]
cp_toKeep$pert_name = paste0(cp_toKeep$pert_iname,' (',cp_toKeep$moa,')')
table(rownames(mtx2) %in% cp_toKeep$pert_name)
rownames(mtx2)[!rownames(mtx2) %in% cp_toKeep$pert_name]

mtx2 = mtx2[,colnames(mtx2) %in% cellline$Cell.Line]

## Make heatmap

colSplit = cellline$group[match(colnames(mtx2),cellline$Cell.Line)]
colSplit[is.na(colSplit)] = 'unknown'
colSplit = factor(colSplit,c('Leukaemia','Solid cancer','non-cancer','unknown'))
logitCols = c('#00008b','#adcae6','white','#ffd1d1','#d80000')

if(signature=='shared'){
  pdf('sharedSig_250325_heatmap.pdf',width = 12,height = 8.5)
  hm = Heatmap(mtx2,name='Similarity to shared up-reg module',na_col = grey(0.8),show_column_dend = F,show_row_dend = F,
               col = circlize::colorRamp2(seq(-2,2,length.out=length(logitCols)),logitCols),
               row_names_gp = gpar(fontsize=9),cluster_rows = F,cluster_columns = F,
               column_names_gp = gpar(fontsize=9),border = T,
               #split = cp_toKeep$group[match(rownames(mtx2),cp_toKeep$pert_name)],
               column_split = colSplit,  heatmap_legend_param = list(direction = "horizontal"))
  
  draw(hm,heatmap_legend_side = "bottom")
  dev.off()  
}else if (signature == 'MLDS'){
  pdf('MLDS_heatmap.pdf',width = 12,height = 7.5)
  hm = Heatmap(mtx2,name='Similarity to MLDS up-reg module',na_col = grey(1),show_column_dend = F,show_row_dend = T,
               col = circlize::colorRamp2(seq(-2,2,length.out=length(logitCols)),logitCols),
               row_names_gp = gpar(fontsize=9),cluster_rows = F,cluster_columns = F,
               column_names_gp = gpar(fontsize=9),border = T,
               split = cp_toKeep$group[match(rownames(mtx2),cp_toKeep$pert_name)],
               column_split = colSplit,  heatmap_legend_param = list(direction = "horizontal"))
  
  draw(hm,heatmap_legend_side = "bottom")
  dev.off() 
}



## As a box plot
cp.sub$dose = as.numeric(gsub(' uM','',cp.sub$pert_idose))
dd = cp.sub %>% group_by(pert_iname,moa) %>% summarise(med = median(norm_cs))
dd = dd[order(dd$med,decreasing = T),]
dd$moa_group = ifelse(dd$med >1,'high',
                      ifelse(dd$med < 1 & dd$med>-0.9,'mid','low'))
cp.sub$moa_group = dd$moa_group[match(cp.sub$pert_iname,dd$pert_iname)]
cp.sub$group = cp_toKeep$group[match(cp.sub$pert_iname,cp_toKeep$pert_iname)]
#cp.sub$pert_iname = paste0(cp.sub$pert_iname,' (',cp.sub$moa,')')

cp.sub$group[is.na(cp.sub$group)] = '-'
cp.sub$group = paste0(cp.sub$group,'::',cp.sub$moa_group)
table(cp.sub$group)
table(is.na(cp.sub$group))
#cp.sub$group['nfkb_inhibitor_dissimilar::mid'] = 'nfkb_inhibitor::mid'
cp.sub$group = factor(cp.sub$group,rev(c('dissimilar::low','dissimilar::mid',
                                     'nfkb_inhibitor_dissimilar::mid',
                                     'nfkb_inhibitor::low','nfkb_inhibitor::mid','nfkb_inhibitor::high',
                                     'nfkb_activator::mid',
                                     'similar::mid','similar::high')))
cp.sub_plot <- cp.sub %>%
  filter(!group %in% c('dissimilar::mid','similar::mid')) %>%
  arrange(norm_cs) %>%
  mutate(pert_iname = factor(pert_iname, levels = unique(pert_iname)))
table(cp.sub_plot$moa[cp.sub_plot$group == 'nfkb_activator::mid'])

plotDir='~/lustre_mt22/SETBP1/manuscriptDraft_0325/Plots'

plotFun_cmap_boxplot = function(noFrame=FALSE,noPlot=FALSE){
  
  p = ggplot(cp.sub_plot,aes(norm_cs,pert_iname))+
    geom_boxplot(aes(fill=group),outlier.shape = NA)+
    geom_quasirandom(height=1,grouponX=T,size=0.6)+
    #scale_fill_manual(values = c('high'=col25[2],'mid'=grey(0.8),'low'=col25[1]))+
    scale_fill_manual(values = c(
      # dissimilar
      "dissimilar::low"  = col25[1],
      "dissimilar::mid"  = "#6BAED6",
      
      # nfkb inhibitor – dissimilar
      "nfkb_inhibitor_dissimilar::mid" = col25[3],
      
      # nfkb inhibitor
      "nfkb_inhibitor::low"  = "#238B45",
      "nfkb_inhibitor::mid"  = "#74C476",
      "nfkb_inhibitor::high" = "#C7E9C0",
      
      # nfkb activator
      "nfkb_activator::mid" = "#9E9AC8",
      
      # similar
      "similar::mid"  = "#D95F5F",
      "similar::high" = col25[2]
    ))+
    theme_classic(base_size = 13)+
    geom_vline(xintercept = c(0),lwd=0.25,lty=2)+
    facet_grid(group~.,scales = 'free',space = 'free')+
    ggtitle('SETBP1 shared signature')+xlab('Normalised similarity score')+ylab('Compound')+
    theme(panel.border = element_blank(),
          axis.line = element_line(colour = 'black',linewidth = 0.2),
          strip.background=element_blank(),
          #legend.position = 'none',
          strip.text.x = element_text(size=10,colour = 'black'),
          strip.text.y = element_text(size=10,colour = 'black'),
          axis.ticks = element_line(colour = 'black',linewidth = 0.2),
          axis.text.x = element_text(size = 10,angle = 0, vjust = 0,hjust = 0.5,colour = 'black'),
          axis.text.y = element_text(size=10),
          axis.text = element_text(colour = 'black'),
          legend.text = element_text(colour = 'black',size = 10),
          legend.title = element_text(colour = 'black',size = 12))+
    scale_y_discrete(labels = function(x) substr(x, 1, 100)) +
    # Add MOA labels on the right
    geom_text(aes(
      x = max(cp.sub_plot$norm_cs) + 0.05 * diff(range(cp.sub_plot$norm_cs)),
      label = moa
    ),
    data = cp.sub_plot,
    hjust = 0,   # left-align so it's horizontal
    vjust = 0.5, # vertically centered on the y-axis level
    size = 3
    )
  
  print(p)
}

saveFig(file.path(plotDir,paste0('Fig2_SETBP1.shared_CMAP.boxplot_nfkb')),
        plotFun_cmap_boxplot,rawData=cp.sub[cp.sub$moa_group!='mid',],
        width = 8.5,height = 9,res = 500)

# Group drugs by MOA -----------------------------------------------------------
cp_agg = cp %>% group_by(pert_id,pert_iname,moa) %>% 
  summarise(norm_cs_median = median(norm_cs, na.rm = TRUE),
            n_cell_lines =n()) %>% 
  mutate(score = norm_cs_median,
         is_nfkb = grepl("nfkb|NFkB|NF-?κB|NFKB", moa, ignore.case = TRUE)) %>% 
  mutate(moa = ifelse(is.na(moa) | moa == "", NA_character_, moa)) %>%
  tidyr::separate_rows(moa, sep = "\\|") %>%
  filter(!is.na(moa) & moa != "")

n_distinct(cp_agg$moa)

cp_per_moa = cp_agg %>% group_by(moa) %>% summarise(med = median(score),n_drug=n())
cp_per_moa = cp_per_moa[order(cp_per_moa$med,decreasing = T),]
cp_per_moa$idx = seq(1:nrow(cp_per_moa))
View(cp_per_moa[cp_per_moa$n_drug >= 3,])
cp_agg.sub = cp_agg[cp_agg$moa %in% cp_per_moa$moa[cp_per_moa$n_drug >= 3],]
p <- ggplot(cp_agg.sub,aes(score,reorder(moa,score, FUN = median)))+
  geom_boxplot(outlier.shape = NA)+
  geom_quasirandom(size=0.6,aes(col=is_nfkb))+
  scale_color_manual(values = c('TRUE'='red','FALSE'=grey(0.7)))+
  scale_fill_manual(values = c('TRUE'='red','FALSE'=grey(0.7)))+
  geom_vline(xintercept = 0)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1))

pdf(file.path(plotDir,'cmap_shared_drugMOA.pdf'),width = 10,height = 40)
print(p)
dev.off()


dt = data.table::setDT(cp)
cp_agg <- dt[
  ,.(
    pert_iname = first(pert_iname),
    moa        = first(moa),
    norm_cs_mean = mean(norm_cs, na.rm = TRUE),
    norm_cs_median = median(norm_cs, na.rm = TRUE),
    n_cell_lines =.N
  ),
  by = pert_id
  ]

# Choose which aggregated score you want to rank by:
cp_agg[, `:=`(
  score = -norm_cs_mean,
  is_nfkb = grepl("nfkb|NFkB|NF-?κB|NFKB", moa, ignore.case = TRUE))] 

library(data.table)
