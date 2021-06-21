### load DESeq package
suppressMessages(require("DESeq2"))
suppressMessages(library(genefilter,quietly=TRUE))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("gplots"))
suppressMessages(library(tidyverse))
###
suppressMessages(library(EDASeq))
suppressMessages(library(Cairo))
suppressMessages(library(fdrtool))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(pvclust))
suppressMessages(library(GSEABase))
suppressMessages(library(clusterProfiler))
suppressMessages(library(pathview))
suppressMessages(library(fgsea))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(annotate))
suppressMessages(library(enrichplot))
suppressMessages(library(topGO))
suppressMessages(library(gridExtra))
# suppressMessages(library(kableExtra))
# suppressMessages(library(knitr))
suppressMessages(library(IHW))
suppressMessages(library(scales))
suppressMessages(library(GetoptLong))
suppressMessages(library(DEGreport))
suppressMessages(library(ReactomePA))
suppressMessages(library(dendsort))
suppressMessages(library(circlize))

quantifyReads <- function(gr, bamlist, nthreads = 8, paired_end = T,allowMultiOverlap=F) {
  GenomicRanges::strand(gr) <- "*"
  saf <- data.frame(GeneID = as.character(gr), 
                    Chr = GenomicRanges::seqnames(gr),
                    Start = GenomicRanges::start(gr),
                    End = GenomicRanges::end(gr), 
                    Strand = GenomicRanges::strand(gr))
  cts <- Rsubread::featureCounts(bamlist, annot.ext = saf, nthreads = nthreads, 
                                 isPairedEnd = paired_end,
                                 allowMultiOverlap = allowMultiOverlap, 
                                 largestOverlap = T, 
                                 requireBothEndsMapped = T)
  cts$counts
}

plotPCAWithSampleNames = function(x, intgroup="condition", ntop=500){
  # fonction qui ameliore l'ACP en ajoutant le nom des echantillons et le % de chaque axe. (origine git igordot)
  suppressMessages(library(RColorBrewer))
  suppressMessages(library(genefilter))
  suppressMessages(library(lattice))
  
  # pca
  rv = rowVars(assay(x))
  select = order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(x)[select,]))
  
  # proportion of variance
  variance = pca$sdev^2 / sum(pca$sdev^2)
  variance = round(variance, 3) * 100
  
  # sample names
  names = colnames(x)
  
  # factor of groups
  fac = factor(apply(as.data.frame(colData(x)[, intgroup, drop=FALSE]), 1, paste, collapse=" : "))
  
  # colors
  if( nlevels(fac) >= 10 )
    colors = rainbow(nlevels(fac))
  else if( nlevels(fac) >= 3 )
    colors = brewer.pal(nlevels(fac), "Set1")
  else
    colors = c( "dodgerblue3", "firebrick3" )
  
  # plot
  xyplot(
    PC2 ~ PC1, groups=fac, data=as.data.frame(pca$x), pch=16, cex=1.5,
    aspect = "fill",
    col = colors,
    xlab = list(paste("PC1 (", variance[1], "%)", sep=""), cex=0.8),
    ylab = list(paste("PC2 (", variance[2], "%)", sep=""), cex=0.8),
    panel = function(x, y, ...) {
      panel.xyplot(x, y, ...);
      ltext(x=x, y=y, labels=names, pos=1, offset=0.8, cex=0.7)
    },
    main = draw.key(
      key = list(
        rect = list(col = colors),
        text = list(levels(fac)),
        rep = FALSE
      )
    )
  )
}

generateDESeqObject = function(counts=NULL,metadata=NULL){
  k27_rse_all <- SummarizedExperiment(assays = list(counts = counts),
                                      rowRanges = GRanges(rownames(counts)))
  k27_rse_all%<>% sort %>% chromVAR::filterPeaks(non_overlapping =T)
  
  seqlevelsStyle(rowRanges(k27_rse_all)) = "UCSC"
  
  gcview <- Biostrings::Views(genome, rowRanges(k27_rse_all))
  gcFrequency <- Biostrings::letterFrequency(gcview, 
                                             letters = "GC", 
                                             as.prob = TRUE) %>% 
    set_colnames("GC")
  
  mcols(k27_rse_all) <- cbind(mcols(k27_rse_all), gcFrequency)
  
  colData=data.frame(ID=metadata$sample,time=as.factor(metadata$time),
                     treatment=as.factor(metadata$treatment),
                     condition = as.factor(metadata$condition)) %>% 
    tibble::column_to_rownames("ID")
  colData$time = relevel(colData$time,ref = "3")
  colData$treatment = relevel(colData$treatment,ref = "nonDiff")
  colData$condition = relevel(colData$condition,ref = "nonDiff:3")
  
  colData(k27_rse_all) <- DataFrame(colData)
  
  
  eda_data <- newSeqExpressionSet(counts = as.matrix(counts(k27_rse_all)), 
                                  featureData = as.data.frame(mcols(k27_rse_all, 
                                                                    use.names = T)), 
                                  phenoData = colData['condition'])
  
  dataOffset <- EDASeq::withinLaneNormalization(eda_data, "GC", 
                                                which = "full", offset = T)
  dataOffset <- EDASeq::betweenLaneNormalization(dataOffset, 
                                                 which = "full", offset = T)
  
  EDASeq::biasPlot(eda_data, "GC", log = TRUE, ylim =c(0,10))
  
  EDASeq::biasPlot(dataOffset, "GC", log = TRUE, ylim = c(0, 10))
  
  EDASeqNormFactors <- exp(-1 * EDASeq::offst(dataOffset))
  EDASeqNormFactors <- EDASeqNormFactors/exp(rowMeans(log(EDASeqNormFactors)))
  
  counts(k27_rse_all) <- as.matrix(counts(k27_rse_all)) # deseq2 wants this to a vanilla matrix
  
  dds <- DESeqDataSet(k27_rse_all, design = ~ condition)
  
  #dds$CONDITION <- relevel(dds$CONDITION, ref = "")
  
  normalizationFactors(dds) <- EDASeqNormFactors
  
  return(DESeq(dds,quiet = T))
}

plotggplot_Heatmap=function(dds=dds,DE=DEgenes.names,export=F,genesOfInterest=genesOfInterest,
                            annot=NULL,condCol='condition',comp = "DE",
                            name="heatmap",colTopAnnot = colTopAnnot, 
                            samplesOrder = NA,split=0,textSize=10,titleSize = 12,
                            direction = NULL,direction_col = NULL,col_fill=NA){
  # dds = DESeq2 object, DE = DE genes table, genesOfInterest = list of genes of interest to annotate in the heatmap,
  # annot = DESeq2 object coldata (sample annotation group), colTopAnnot = column annotation,
  # samplesOrder = preferential samples order,  split = ENSG|SYMBOL where to find, Symbol within the gene_ids,
  # 
  
  # print("DEG Clustering...")
  bootstrap=FALSE
  
  if(is.null(annot)){
    annot=colData(dds)[condCol]
  }
  annot=data.frame(annot,stringsAsFactors = T)
  if(is.null(colTopAnnot)){
    colTopAnnot<-vector("list", ncol(annot))
    names(colTopAnnot)<-colnames(annot)
    colFun<-c(ggplotColours,rainbow)
    i<-1
    for(col in colnames(annot)){
      colTopAnnot[[col]]<-colFun[[i]](nlevels(annot[,col]))
      names(colTopAnnot[[col]])<-levels(annot[,col])
      i<-i+1
      if(i==4) i<-1
    }
  }
  
  ha<-HeatmapAnnotation(df=annot,col = colTopAnnot,name = condCol)
  
  DEgenes.names=DE
  colFun<-c(ggplotColours,rainbow)
  rld=rlogTransformation(dds)
  exprDatT=assay(rld)
  sampleHt<-colnames(exprDatT)
  haByComp<-ha
  exprDE=exprDatT[DEgenes.names,sampleHt,drop=FALSE]
  exprDE.scaled=rowScale(exprDE,center=T,scaled=T)
  if(split>0){
    rownames(exprDE.scaled)=sapply(rownames(exprDE.scaled), function(x){unlist(strsplit(x,split='|',fixed=T))[split]})
  }
  
  bootTemp=bootstrap
  if(nrow(exprDE.scaled>10)){bootTemp=FALSE}
  
  hclustGeneDE<-unsupervisedClustering(exprDE.scaled,transpose = F,nboot=nboot,bootstrap = bootTemp)
  hclustSampleDE<-unsupervisedClustering(exprDE.scaled,transpose = T,nboot=nboot,bootstrap = bootTemp,method.dist = "euclidean")
  ddSampleDE = as.dendrogram(hclustSampleDE)
  
  rowScaledExpr<-rowScale(exprDatT,center=TRUE,scale=TRUE)
  quantile.expr<-quantile(unlist(rowScaledExpr),seq(0,1,.01),na.rm = T)
  colHA=colorRamp2(c(quantile.expr[2],0,quantile.expr[100]),c("blue","white","red"))
  if(!is.null(direction)){
    direction = subset(direction,rownames(direction) %in% rownames(exprDE.scaled))
  }
  
  if(export==F & !is.null(direction)){
    htlist=Heatmap(exprDE.scaled,col = colHA,name=comp,
                   column_title = paste0("(n = ",nrow(exprDE),")"),
                   column_title_gp = gpar(fontsize = titleSize),
                   show_column_names = F,
                   heatmap_legend_param = list(title='Row Z score'),top_annotation = haByComp,
                   # cluster_columns = as.hclust(dendsort(as.dendrogram(hclustSampleDE))),
                   cluster_columns = reorder(ddSampleDE,wts = order(match(samplesOrder,colnames(exprDE.scaled)))),
                   cluster_rows = hclustGeneDE) +
      Heatmap(direction,name="direction",col = direction_col) + 
      rowAnnotation(link=anno_mark(at=which(rownames(exprDE.scaled)%in%genesOfInterest),
                                   labels=rownames(exprDE.scaled)[rownames(exprDE.scaled)%in%genesOfInterest],
                                   labels_gp = gpar(fontsize = textSize), padding = unit(3, "mm")))
    
    return(htlist)
  }
  if(export==F & is.null(direction)){
    htlist=Heatmap(exprDE.scaled,col = colHA,name=comp,
                   column_title = paste0("(n = ",nrow(exprDE),")"),
                   column_title_gp = gpar(fontsize = titleSize),
                   show_column_names = F,
                   heatmap_legend_param = list(title='Row Z score'),top_annotation = haByComp,
                   # cluster_columns = as.hclust(dendsort(as.dendrogram(hclustSampleDE))),
                   cluster_columns = reorder(ddSampleDE,wts = order(match(samplesOrder,colnames(exprDE.scaled)))),
                   cluster_rows = hclustGeneDE)+
      rowAnnotation(link=anno_mark(at=which(rownames(exprDE.scaled)%in%genesOfInterest),
                                   labels=rownames(exprDE.scaled)[rownames(exprDE.scaled)%in%genesOfInterest],
                                   labels_gp = gpar(fontsize = textSize,
                                                    col = col_fill), padding = unit(5, "mm")))
    
    return(htlist)
  }
}

degPlotLocal = function (dds, xs, res = NULL, n = 9, genes = NULL, 
                         group = NULL, batch = NULL, 
                         metadata = NULL, 
                         ann = c("geneID", "symbol"), 
                         slot = 1L, 
                         log2 = TRUE, 
                         xsLab = xs, 
                         ysLab = "abundance", 
                         color = "black", groupLab = group, batchLab = batch,
                         splitGeneName=0,reverseColor=T) 
  # modified degPlot function from DEGreport package
{
  if (class(dds)[1] %in% c("data.frame", "matrix")) 
    dds = SummarizedExperiment(assays = SimpleList(counts = as.matrix(dds)), 
                               colData = metadata)
  
  if (!("assays" %in% slotNames(dds))) 
    stop("dds object doesn't have assays slot")
  
  if (!is.null(res)) 
    res <- res[order(res$padj), ] %>% .[!is.na(res$padj),]
  
  if (is.null(genes))
    genes =  row.names(res)[1L:n]
  
  anno <- as.data.frame(rowData(dds))
  
  metadata = data.frame(colData(dds))
  if (class(dds)[1] == "DESeqDataSet")
    counts <- counts(dds, normalized = TRUE)
  else counts <- assays(dds)[[slot]]
  
  stopifnot(class(counts)[1] == "data.frame" | class(counts)[1] == "matrix")
  
  if (log2 & max(counts) < 500L)
    warning("Data seems to be already in log2. Please use log2 = FALSE.")
  if (log2){
    counts <- log2(counts + 0.2)
    ysLab <- paste("log2", ysLab)
  }
  
  newgenes <- genes
  if (ncol(anno) > 0) {
    name <- intersect(names(anno), ann)
    if (length(name) != 2L)
      message("No genes were mapped to rowData. check ann parameter values.")
    if (length(name) == 2L)
      newgenes <- anno[match(genes, anno[, ann[1L]]), ann[2L]]
    if (sum(is.na(newgenes)) > 0)
      warning(sum(is.na(newgenes)), " cannot be mapped to ", name[2L], 
              ". Those will be skipped.")
  }
  
  dd <- reshape::melt(as.data.frame(counts[genes, , drop = FALSE]) %>%
                        mutate(gene = factor(newgenes, levels=newgenes)))
  colnames(dd) = c("gene", "sample", "count")
  
  dd <-dd[!is.na(dd[["gene"]]),]
  
  dd$xs = as.factor(metadata[as.character(dd$sample), xs])
  
  if (!is.null(group)) {
    dd[, groupLab] = as.factor(metadata[as.character(dd$sample), group])
  }else {
    groupLab = "fake"
    dd[, groupLab] = "fake"
  }
  
  if(splitGeneName>0){
    
    dd[,'gene'] = sapply(as.vector(dd[,'gene']), function(x){
      unlist(strsplit(x,split="|",fixed = T))[2]
    })
  }
  
  if (!is.null(batch)) {
    dd[, batchLab] = as.factor(metadata[as.character(dd$sample), batch])
    
    p = ggplot(dd, aes_string(x = "xs", y = "count", color = groupLab,
                              shape = batchLab))
  }else{
    p = ggplot(dd, aes_string(x = "xs", y = "count", color = groupLab))
  }
  
  p = p +
    # geom_violin(alpha=0.3) +
    stat_smooth(fill = "grey80", method = 'loess') +
    geom_point(size = 1, alpha = 0.7,
               position = position_jitterdodge(dodge.width = 0.9)) +
    facet_wrap(~gene, scales = "free_y") +
    xlab(xsLab) +
    ylab(ysLab)
  if (length(unique(dd[, groupLab])) == 1L) {
    stopifnot(length(color) == 1)
    p = p +
      scale_color_manual(guide = FALSE, values = rev(color)) +
      scale_fill_manual(guide = FALSE, values = rev(color))
  }else{
    if (color == "black")
      color = "Set1"
    if (reverseColor){
      p = p +
        scale_color_brewer(palette = color,direction = -1) +
        scale_fill_brewer(palette = color,direction = -1)
    }else{
      p = p +
        scale_color_brewer(palette = color) +
        scale_fill_brewer(palette = color)
    }
    
  }
  p = p + theme_bw() +
    theme(strip.background = element_rect(fill = "white"),
          strip.text = element_text(colour = "black"))
  
  suppressWarnings(p)
}

plotEnrichFactor=function(enrGeneSet,n=20,title="Enriched Gene Sets",size = 12,legendSize = 12,titleSize = 15){
  library(ggplot2)
  library(forcats)
  library(enrichplot)
  # richFactor number of genes in common divided by number of genes in geneset
  enrGeneSet=clusterProfiler.dplyr::mutate(enrGeneSet,richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
  ggplot(enrGeneSet, showCategory = n, 
         aes(richFactor, fct_reorder(Description, richFactor))) + 
    geom_segment(aes(xend=0, yend = Description)) +
    geom_point(aes(color=p.adjust, size = Count)) +
    # scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
    scale_color_gradient(low = "red",high = "blue",guide = guide_colorbar(reverse=TRUE))+
    scale_size_continuous(range=c(2, 10)) +
    theme_minimal() +
    xlab("rich factor\n(count/geneset size)") +
    ylab(NULL) + 
    ggtitle(title)+ 
    theme(axis.text.y = element_text(size = size),
          axis.text.x = element_text(size = size),
          axis.title.x = element_text(size = size),
          legend.text = element_text(size = legendSize),
          plot.title = element_text(size = titleSize))
}