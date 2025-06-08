#' An R wrapper around the cellTypist implementation of logistic regression
#'
#' This has been tuned in such a way as tobe quasi-equivalent to running logistic regression as per the old logisticRegression.R script, which implements logistic regression as performed in our papers starting Young et al. Science, 2018.  What that means in practice is that I consider the output of logistic regression to be a table of logits predicting the probability of each cell (row) matching to the reference classes in the trained model (columns).
#'
#' This table of logits is the starting point for further downstream analysis as needed.  In particular, this script also includes functions for summarisation and visualisation.  
#'
#' Note that the cellTypist call defaults are purposely defined to be minimal, but can pass through parameters to the underlying cellTypist functions.  You will almost always want to do this for at least one or two parameters, such as specifying n_jobs when training so it doesn't take forever.
#'
#' Obviously, in order for this to work cellTypist needs to be install.


#' Train a logistic regression model with celltypist
#'
#' Passes (using files) the counts and cell type labels to celltypist.  There are many parameters to the celltypist.train function, which can be set using the ...
#'
#' @param cnts Sparse matrix to use to train the model.
#' @param outPath Where to save the model.  Should end in pkl.
#' @param labels The cell type labels to use for model training.
#' @param ... Extra parameters to pass to celltypist.train.
#' @return Nothing, just saves the model where we should.
trainCelltypistModel = function(cnts,outPath,labels=colnames(cnts),...){
  #Create three temp files and write out matrix
  base = tempfile()
  mtxFile = paste0(base,'.mtx')
  gnsFile = paste0(base,'_gns.tsv')
  labFile = paste0(base,'_cls.tsv')
  writeMM(cnts,mtxFile)
  write.table(rownames(cnts),gnsFile,sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE)
  write.table(labels,labFile,sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE)
  #Extra params
  theDots = list(...)
  #Default for transpose_input
  if(!'transpose_input' %in% names(theDots))
    theDots$transpose_input = TRUE
  if(length(theDots)==0){
    extraParams=''
  }else{
    if(any(lengths(theDots))!=1)
      stop("Any extra characters must be of length 1.")
    extraParams = paste0(names(theDots),'=',sprintf(ifelse(sapply(theDots,is.character),"'%s'",'%s'),as.character(unlist(theDots))))
    #Fix up bools into python format
    extraParams = gsub('FALSE','False',gsub('TRUE','True',extraParams))
    #Wrap it up
    extraParams = paste0(',',paste(extraParams,collapse=','))
  }
  #The python code to run
  tmpdir=tempdir()
  cmd = sprintf("import celltypist\nfit=celltypist.train('%s',labels='%s',genes='%s'%s)\nfit.write('%s')",mtxFile,labFile,gnsFile,extraParams,outPath)
  system(paste0('python -c "',cmd,'"'))
  unlink(c(mtxFile,gnsFile,labFile))
}


#' Run celltypist
#'
#' Export a count matrix via a file to python, run celltypist, then return the results.
#'
#' @param cnts The count matrix (in sparse matrix format) that we want to call things on.
#' @param ... Passed to the celltypist.annotate call.
#' @return A list containing the raw logit probabilities (logitMat) and the label predictions (labMat)
runCelltypist = function(cnts,...){
  #Create three temp files and write out matrix
  base = tempfile()
  mtxFile = paste0(base,'.mtx')
  gnsFile = paste0(base,'_gns.tsv')
  clsFile = paste0(base,'_cls.tsv')
  writeMM(cnts,mtxFile)
  write.table(rownames(cnts),gnsFile,sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE)
  write.table(colnames(cnts),clsFile,sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE)
  #Extra params
  theDots = list(...)
  #Default for transpose_input
  if(!'transpose_input' %in% names(theDots))
    theDots$transpose_input = TRUE
  if(length(theDots)==0){
    extraParams=''
  }else{
    if(any(lengths(theDots))!=1)
      stop("Any extra characters must be of length 1.")
    extraParams = paste0(names(theDots),'=',sprintf(ifelse(sapply(theDots,is.character),"'%s'",'%s'),as.character(unlist(theDots))))
    #Fix up bools into python format
    extraParams = gsub('FALSE','False',gsub('TRUE','True',extraParams))
    #Wrap it up
    extraParams = paste0(',',paste(extraParams,collapse=','))
  }
  #The python code to run
  tmpdir=tempdir()
  cmd = sprintf("import celltypist\npreds=celltypist.annotate('%s',gene_file='%s',cell_file='%s'%s)\npreds.to_table('%s')",mtxFile,gnsFile,clsFile,extraParams,tmpdir)
  system(paste0('python -c "',cmd,'"'))
  #Load output 
  logitMat = read.table(file.path(tmpdir,"decision_matrix.csv"),header=TRUE,sep=',')
  rownames(logitMat) = logitMat$X
  logitMat$X=NULL
  labMat = read.table(file.path(tmpdir,"predicted_labels.csv"),header=TRUE,sep=',')
  unlink(c(mtxFile,gnsFile,clsFile))
  return(list(logitMat=logitMat,labMat=labMat))
}


#' Collapse logit probabilities to clusters
#'
#' Given a logistic regression fit table at the individual cell level, collapses the cell probabilities into clusters.  This is done by first converting the input into the probability space (if needed) to prevent the averaging skewing things based on high information outliers in logit space.  The results are then optionally re-converted to logits.
#'
#' In addition to collapsing and averaging the heatmap, this function can also create a list of lists for use making the fancy heatmap.
#'
#' @param fit Logistic regression fit object where rows are cells and columns are the LR classes.  Can be in either logit or probability format.
#' @param classes A vector that is either the same length as the number of rows in \code{fit}, or is named by the cell names.
#' @param outputLogits Convert output back to logits?  If NULL, preserve whatever was inputted.
#' @param outputForFancyHeatmap Should we output a list containg paired objects to use plotting a fancy heatmap?
#' @return A table where the rows have been collapsed by averaging to cluster level.
collapseToClusters = function(fit,clusts,outputLogits=NULL,outputForFancyHeatmap=FALSE){
  #Convert input if logit space
  if(any(fit<0 | fit>1)){
    inputLogit=TRUE
    fit = (1+exp(-fit))**-1
  }else{
    inputLogit=FALSE
  }
  if(is.null(outputLogits))
    outputLogits = inputLogit
  #Get the class vector
  if(!is.null(names(clusts)) && all(rownames(fit) %in% names(clusts))){
    clVec = clusts[rownames(fit)]
  }else if(length(clusts)==nrow(fit)){
    clVec = clusts
  }else{
    stop("clusts not of appropriate length and names do not match rownames of fit")
  }
  #Do the averaging
  s = split(seq(nrow(fit)),clVec)
  out = lapply(s,function(e) colSums(fit[e,,drop=FALSE])/length(e))
  out = do.call(rbind,out)
  rownames(out) = names(s)
  if(outputLogits)
    out = log(out)-log(1-out)
  if(!outputForFancyHeatmap)
    return(out)
  out = list(mat = out)
  #This is a bit of a magic incantation, but the logic is that we first split by cluster, then within each cluster we split by the reference classes and sort each entry.  This are then coeerced to appear in the same order as in the averaged output
  if(outputLogits){
    tmp = log(fit) - log(1-fit)
  }else{
    tmp = fit
  }
  out$x = lapply(s,function(e) setNames(lapply(seq(ncol(tmp)),function(ee) sort(tmp[e,ee])),colnames(tmp))[colnames(out$mat)])[rownames(out$mat)]
  return(out)
}

#' Collapse classes
#'
#' Groups together reference class probabilities into one category by taking the maximum of all probabilities in a group.
#'
#' @param fit Logistic regression fit objecw where rows are cells and columns are the LR classes.  Can be in either logit or probability format.
#' @param groups A vector with either length equal to the number of columns in \code{fit} or a named vector where the names are existing column names of \code{fit}.  If a column of \code{fit} i- not present, it is assumed to not need collapsing. 
#' @param collapseFun The function used to collapse across columns.  You likely want \code{max} if softmax has not been applied to \code{fit} (the default) and \code{sum} if it has.
#' @param outputLogits Convert output back to logits?  If NULL, preserve whatever was inputted.
#' @return Collapsed version of \code{fit}
collapseClasses = function(fit,groups,collapseFun=max,outputLogits=NULL){
  #Convert input to prob space if not already
  if(all(fit>=0 & fit<=1)){
    inputLogit=FALSE
    #fit = log(fit)-log(1-fit)
  }else{
    inputLogit=TRUE
    fit = (1+exp(-fit))**-1
  }
  if(is.null(outputLogits))
    outputLogits=inputLogit
  #Get the class vector
  if(!is.null(names(groups)) && all(names(groups) %in% colnames(fit))){
    grpVec = colnames(fit)
    grpVec[match(names(groups),colnames(fit))]=groups
  }else if(length(groups)==ncol(fit)){
    grpVec = groups
  }else{
    stop("groups not of appropriate length and names do not match colnames of fit")
  }
  #Do the aggregation
  s = split(seq(ncol(fit)),grpVec)
  out = lapply(s,function(e) apply(fit[,e,drop=FALSE],1,collapseFun))
  out = do.call(cbind,out)
  colnames(out) = names(s)
  #Convert to logit space if that's what's required
  if(outputLogits)
    out = log(out)-log(1-out)
  #Convert to prob space if that's what you've bene asked for
  #if(!outputLogits)
  #  out = (1+exp(-out))**-1
  return(out)
}
 

#' Plot similarity score by group
#'
#' Given a collection of cells, their logistic regression scores against some referenc, and a way to group them together, plot the similarity in groups.
#'
#' Each group is represented by a bar indicating the similarity score for each cell in that grouping, with granularity determined by \code{Nchunks}.  This is built using ComplexHeatmaps.
#'
#' The assumption is that \code{fit} will be in probability space.  If they are logit space you need to either convert them to probability space or change the colorSpace parameters.
#'
#' Note that if you wish to pass row_split or otherwise work annotate or extend the row names of the collapsed heatmap, rows will appear in sorted order from \code{cellGroupings}.  That is, row names are equal to \code{sort(unique(cellGroupings))}.  You can also manually set the row ordering using the \code{row_order} parameter.
#'
#' The default is now that \code{colSpaceBreaks} is automatically set to \code{c(0,0.5,1)} if the input values are probabilities and \code{c(-5,0,5)} if they are logits.  You can change this by manually setting \code{colSpaceBreaks} to something other than \code{NULL}.
#'
#' @param fit Logistic regression fit object where rows are cells and columns are the LR classes.  Can be in either logit or probability format.
#' @param cellGroupings Vector with length matching the number of rows in \code{fit} defining how cells are to be grouped.
#' @param colSpaceBreaks The range of values for colours to represent, with breaks where colours are defined as per \code{\link{colorRamp2}}.  If NULL, automatically set based on input type.
#' @param colSpaceCols The colours for colour scale as per \code{\link{colorRamp2}}.
#' @param colSpaceBins The colour space values will be split into bins for computational efficiency.  This is the number of bins to create, spaced evenly across \code{range(colSpaceBreaks)}.
#' @param colSpaceTruncate If values fall outside the specified colour space, they will be replaced with NAs by default.  If this is set to TRUE, they are instead collapsed to the closest boundary of the colour space.
#' @param useLogits
#' @param bgColor Color used for background, which will only be visible for NA squares or in gaps between squares.
#' @param ... Extra options to pass to Heatmap.  
#' @return The output of the call to Heatmap.
plotByGroup = function(fit,cellGroupings,colSpaceBreaks=NULL,colSpaceCols =c('#2D3F90','white','#EA2335'),colSpaceBins=20,colSpaceTruncate=TRUE,bgColor='grey',...){
  #Detect the input type
  inputLogit = !(all(fit>=0 & fit<=1))
  #Auto-set if needed
  if(is.null(colSpaceBreaks)){
    if(inputLogit){
      colSpaceBreaks = c(-5,0,5)
    }else{
      colSpaceBreaks = c(0,.5,1)
    }
  }
  #Check sanity of input
  if(length(colSpaceBreaks)!=length(colSpaceCols))
    stop(sprintf('Length of colSpaceBreaks is %d, which does not match length of colSpaceCols, which is %d',length(colSpaceBreaks),length(colSpaceCols)))
  #Group data
  collapsedDat = collapseToClusters(fit,cellGroupings,outputForFancyHeatmap=TRUE)
  #Define the colour scheme
  col_fun = circlize::colorRamp2(colSpaceBreaks,colSpaceCols)
  colSpaceMin = min(colSpaceBreaks)
  colSpaceMax = max(colSpaceBreaks)
  #Fix colour space if required
  if(colSpaceTruncate){
    collapsedDat$mat[collapsedDat$mat < colSpaceMin] = colSpaceMin
    collapsedDat$mat[collapsedDat$mat > colSpaceMax] = colSpaceMax
    for(i in seq_along(collapsedDat$x)){
      for(j in seq_along(collapsedDat$x[[i]])){
        collapsedDat$x[[i]][[j]][collapsedDat$x[[i]][[j]] < colSpaceMin] = colSpaceMin
        collapsedDat$x[[i]][[j]][collapsedDat$x[[i]][[j]] > colSpaceMax] = colSpaceMax
      }
    }
  }
  #Define the cell plotting function
  cell_fun = function(j, i, x, y, width, height, fill) {
    #Define number of boxes to split into
    N=colSpaceBins
    #Do tons of teny tiny rectangles
    dat = collapsedDat$x[[i]][[j]]
    #dat = (1+exp(-sortDat[[i]][[j]]))**-1
    #White out the area
    grid.rect(x = x, y = y, width = width, height = height,gp = gpar(col = NA, fill = bgColor))
    #Draw N boxes per thingy
    #Have a bit less width available for this than the whole region
    subWidth = width*0.95
    subHeight = height*1.0
    sDat = split(dat,cut(dat,seq(colSpaceMin,colSpaceMax,length.out=N+1),include.lowest=TRUE))
    boxWidths = lengths(sDat)/length(dat)
    for(ii in seq_along(sDat)){
      grid.rect(x = (x-subWidth/2) + (sum(boxWidths[which(seq_along(boxWidths)<=(ii-1))]) + boxWidths[ii]/2)*subWidth,
                y = y,
                width = boxWidths[ii]*subWidth,
                height = subHeight,
                gp = gpar(fill=col_fun(mean(sDat[[ii]])),col=NA)
                )
    }
    #Draw the bounding box for the whole thing
    grid.rect(x = x, y = y, width = width, height = height,gp = gpar(col = bgColor, fill = NA))
  }
  hmParams = list(matrix = collapsedDat$mat,
                  name='Similarity',
                  col=col_fun,
                  cell_fun = cell_fun,
                  show_row_dend=FALSE,
                  show_column_dend=FALSE)
  theDots = list(...)
  for(nom in names(theDots))
    hmParams[[nom]] = theDots[[nom]]
  hm = do.call(Heatmap,hmParams)
  return(hm)
}


