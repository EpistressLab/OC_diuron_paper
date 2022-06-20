
#####
# Functions
#####
genome <- BSgenome.Hsapiens.UCSC.hg19
importRoseEnhancerFile=function(fileName,onlySupers = T){
  se=read_tsv(file = fileName,comment = '#',col_names = T)
  if (onlySupers) {
    return(se%>%filter(isSuper==1))
  }else{
    return(se)
  }
  
}

## ----quantify-all-primary-se, cache=T, message=F, warning=F, results='hide'----
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

plotPCAWithSampleNames = function(x, intgroup="condition", ntop=500)
{
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

metadata = read_tsv('Z:/Robel/ChIP_96samples/aligned/OC_data/hg19_analysisSet/superAnalysis/metadata.csv')

hkgList = readxl::read_xlsx('Z:/Robel/scRNA_HOS/ReAnalysis/aad0501_Table_S16.xlsx') %>% 
  `colnames<-`('geneName') %>% 
  pull(geneName)
hkgList = hkgList[!hkgList %in% 'PRPS1L3'] %>% str_replace(pattern = 'HPRT',replacement = 'HPRT1')

genecode=rtracklayer::import.gff("Z:/Bioinfo/Annotations/human/gencode_GRCh37/gencode.v35lift37.basic.annotation.gtf.bgz")


tads_H1ESC_2015Dixon = rtracklayer::import.bed("Z:/Bioinfo/hg19.TADs/H1-ESC_Dixon2015-raw_TADs.txt")
tads_K562_Lieberman = rtracklayer::import.bed("Z:/Bioinfo/hg19.TADs/K562_Lieberman-raw_TADs.txt")
tads_KBM7_Lieberman = rtracklayer::import.bed("Z:/Bioinfo/hg19.TADs/KBM7_Lieberman-raw_TADs.txt")
tadsInfo = c(tads_KBM7_Lieberman,tads_H1ESC_2015Dixon,tads_K562_Lieberman)

hkProm_genecode = promoters(genecode[genecode$gene_name %in% hkgList],upstream = 2000,downstream = 500,use.names = T)
names(hkProm_genecode) = paste0(hkProm_genecode$gene_name,"|",hkProm_genecode$transcript_id)

#####

seList_filt = lapply(metadata_filt$roseFile,function(x){
  importRoseEnhancerFile(fileName = x,onlySupers = T)}) %>% 
  bind_rows()


regionsAll = bind_rows(seList_filt,hkgPromoters %>% 
                         data.frame() %>% 
                         rename(CHROM = seqnames,START = start,STOP = end,REGION_ID = combined))
regionsAll = GRanges(seqnames = regionsAll$CHROM,
                     ranges = IRanges(start = regionsAll$START,end = regionsAll$STOP),
                     regionID =regionsAll$REGION_ID)

ctsRaw_seFilt = quantifyReads(gr = regionsAll,bamlist = metadata_filt$ChIP_bam,
                              nthreads = 12,paired_end = T,allowMultiOverlap = T)
#####
# ==========     _____ _    _ ____  _____  ______          _____  
# =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
# =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
#   ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
#   ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
#   ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
#   Rsubread 2.0.1
# 
# //========================== featureCounts setting ===========================\\
# ||                                                                            ||
#   ||             Input files : 6 BAM files                                      ||
#   ||                           o OC.J3.-.pat2_merged.q20.noDup.bam              ||
#   ||                           o OC.J7.+.pat2_merged.q20.noDup.bam              ||
#   ||                           o OC.J10.+.pat2_merged.q20.noDup.bam             ||
#   ||                           o OC.J3.-.pat1_merged.q20.noDup.bam              ||
#   ||                           o OC.J7.+.pat1_merged.q20.noDup.bam              ||
#   ||                           o OC.J10.+.pat1_merged.q20.noDup.bam             ||
#   ||                                                                            ||
#   ||              Annotation : R data.frame                                     ||
#   ||      Dir for temp files : .                                                ||
#   ||                 Threads : 12                                               ||
#   ||                   Level : meta-feature level                               ||
#   ||              Paired-end : yes                                              ||
#   ||      Multimapping reads : counted                                          ||
#   || Multi-overlapping reads : counted                                          ||
#   ||   Min overlapping bases : 1                                                ||
#   ||                                                                            ||
#   ||          Chimeric reads : counted                                          ||
#   ||        Both ends mapped : required                                         ||
#   ||                                                                            ||
#   \\============================================================================//
# ||    Features : 4102                                                         ||
# ||    Meta-features : 4066                                                    ||
# ||    Chromosomes/contigs : 23

# BAM_file                                 |  total alignements  |  assignedAligments  |  assigned (%)
# OC.J3.-.pat2_merged.q20.noDup.bam        |  11599962           |  610951             | 5.3
# OC.J7.+.pat2_merged.q20.noDup.bam           14143284              1487262              10.5
# OC.J10.+.pat2_merged.q20.noDup.bam          13874169              2031869              14.6
# OC.J3.-.pat1_merged.q20.noDup.bam           11726795              897632               7.7
# OC.J7.+.pat1_merged.q20.noDup.bam           13850531              1503353              10.9
# OC.J10.+.pat1_merged.q20.noDup.bam          14718723              2240456              15.2
#####

counts <- ctsRaw_seFilt %>% set_colnames(metadata_filt$sample)
dds_filt = generateDESeqObject(counts = (counts %>% 
                                           data.frame %>% 
                                           rownames_to_column("SE") %>% 
                                           filter(!SE %in% (hkgPromoters %>% 
                                                              data.frame %>% 
                                                              str_glue_data("{seqnames}:{start}-{end}"))) %>% 
                                           column_to_rownames("SE") %>% 
                                           as.matrix),
                               metadata = metadata_filt)
dds_tmp = dds_filt
de_enh = results(dds_filt,name="condition_Diff.10_vs_nonDiff.3") %>% 
  data.frame() %>% 
  rownames_to_column('enhancer') %>% 
  filter(padj <= 0.1 & log2FoldChange <= -1)
de_gene = (results(dds_OC_diff_gene,contrast = c("condition","OC_J10_diff","OC_J3_min"),test = "Wald") %>% 
             data.frame %>% 
             rownames_to_column("gene") %>% 
             filter(padj <= 0.05 & log2FoldChange <= -1) %>% 
             pull(gene))
renames_down=sapply(rownames(dds_tmp), function(x){
  if(x %in% de_enh$enhancer){
    targets = GRanges(x)%>% 
      subsetByOverlaps(x = tads_H1ESC_2015Dixon) %>% 
      subsetByOverlaps(x = ncbi37.2) %>% 
      data.frame() %>% 
      filter(name %in% de_gene) %>% 
      arrange(name) %>% pull(name) %>% unique()
    if(length(targets)>0){
      return(str_glue("{x} -> {str_flatten(targets,collapse='|')}"))
    }
  }
})
rownames(dds_tmp)[match(names(Filter(Negate(is.null), renames_down)),rownames(dds_tmp))]=as.character(Filter(Negate(is.null), renames_down))

de_enh = results(dds_filt,name="condition_Diff.10_vs_nonDiff.3") %>% 
  data.frame() %>% 
  rownames_to_column('enhancer') %>% 
  filter(padj <= 0.1 & log2FoldChange >= 1)
de_gene = (results(dds_OC_diff_gene,contrast = c("condition","OC_J10_diff","OC_J3_min"),test = "Wald") %>% 
             data.frame %>% 
             rownames_to_column("gene") %>% 
             filter(padj <= 0.05 & log2FoldChange >= 1) %>% 
             pull(gene))
renames_up=sapply(rownames(dds_tmp), function(x){
  if(x %in% de_enh$enhancer){
    targets = GRanges(x)%>% 
      subsetByOverlaps(x = tads_H1ESC_2015Dixon) %>% 
      subsetByOverlaps(x = ncbi37.2) %>% 
      data.frame() %>% 
      filter(name %in% de_gene) %>% 
      arrange(name) %>% pull(name) %>% unique()
    if(length(targets)>0){
      return(str_glue("{x} -> {str_flatten(targets,collapse='|')}"))
    }
  }
})
rownames(dds_tmp)[match(names(Filter(Negate(is.null), renames_up)),rownames(dds_tmp))]=as.character(Filter(Negate(is.null), renames_up))



png("Z:/Robel/OC_diff_diuron_paper/figures_genecode/nfatc1_browserplot_chr_bw.png",height = 9,width = 10,units = "cm",res=600)
pageCreate(width = 10, height = 9, default.units = "cm",showGuides = T)
## Set the coordinates
# plotGG(plot = corr_test,x=10.2,y = 0,just = c("top","left"),width = 9.8,height = 8,default.units = "cm")
# plotGG(plot = corr_test,x=13,y = 8.5,just = c("top","left"),width = 9.5,height = 5.5,default.units = "cm")
# chr18:77,100,195-77,510,622
params <- pgParams(
  chrom = "chr18",
  chromstart = 77100000, chromend = 77510600,
  assembly = "hg19",x=1,
  width = 9,just=c('bottom',"left"),default.units="cm",range=c(0,250)
)
genes = plotGenes(params = params,geneOrder = c("NFATC1"),geneHighlights = data.frame(gene=c("NFATC1"),color=c("red")),
                  y = 8,height = 0.8,strandLabels = F)
annoGenomeLabel(fontsize = 10,params = params,y = 9,height=0.8,genes)
signal1 = plotSignal(data = "Z:/Robel/OC_diff_diuron_paper/deeptools/OC.J3.-.pat1_merged_chr.bw",
                     params = params,y=7,height = 0.6,fill = rgb(maxColorValue = 255,red = 67,green = 147,blue = 195),
                     linecolor = rgb(maxColorValue = 255,red = 67,green = 147,blue = 195)
)
# annoYaxis(signal1,fontsize=8,params = params,fontcolor="black",axisLine = T)
signal2 = plotSignal(data = "Z:/Robel/OC_diff_diuron_paper/deeptools/OC.J3.-.pat2_merged_chr.bw",
                     params = params,y=6,height = 0.6,fill = rgb(maxColorValue = 255,red = 67,green = 147,blue = 195),
                     linecolor = rgb(maxColorValue = 255,red = 67,green = 147,blue = 195)
)
# annoYaxis(signal2,fontsize=8,params = params,fontcolor="black",axisLine = T)
signal3 = plotSignal(data = "Z:/Robel/OC_diff_diuron_paper/deeptools/OC.J7.+.pat1_merged_chr.bw",
                     params = params,y=5,height = 0.6,fill = "#FF7F00",
                     linecolor = "#FF7F00"
)
# annoYaxis(signal3,fontsize=8,params = params,fontcolor="black",axisLine = T)
signal4 = plotSignal(data = "Z:/Robel/OC_diff_diuron_paper/deeptools/OC.J7.+.pat2_merged_chr.bw",
                     params = params,y=4,height = 0.6,fill = "#FF7F00",
                     linecolor = "#FF7F00"
)
# annoYaxis(signal4,fontsize=8,params = params,fontcolor="black",axisLine = T)
signal5 = plotSignal(data = "Z:/Robel/OC_diff_diuron_paper/deeptools/OC.J10.+.pat1_merged_chr.bw",
                     params = params,y=3,height = 0.6,fill = rgb(maxColorValue = 255,red = 255,green = 0,blue = 0),
                     linecolor = rgb(maxColorValue = 255,red = 255,green = 0,blue = 0)
)
# annoYaxis(signal5,fontsize=8,params = params,fontcolor="black",axisLine = T)
signal6 = plotSignal(data = "Z:/Robel/OC_diff_diuron_paper/deeptools/OC.J10.+.pat2_merged_chr.bw",
                     params = params,y=2,height = 0.6,fill = rgb(maxColorValue = 255,red = 255,green = 0,blue = 0),
                     linecolor = rgb(maxColorValue = 255,red = 255,green = 0,blue = 0)
)
annoYaxis(signal6,fontsize=8,params = params,fontcolor="black",axisLine = T)

annoHighlight(signal1,chrom = "chr18",chromstart = 77370732,chromend = 77403785 ,fill = "yellow",alpha = 0.3,
              params = params,y = 7,height = 5.6)

plotText(label = "SE",x = 7.25,y = 1,params = params,just = c("top","center"),fontsize = 10)
plotText(label = "d.3 - Diff pat.1",y=6.2,just = c("top","left"),x = 1.2,params = params,fontsize = 10)
plotText(label = "d.3 - Diff pat.2",y=5.2,just = c("top","left"),x = 1.2,params = params,fontsize = 10)
plotText(label = "d.7 + Diff pat.1",y=4.2,just = c("top","left"),x = 1.2,params = params,fontsize = 10)
plotText(label = "d.7 + Diff pat.2",y=3.2,just = c("top","left"),x = 1.2,params = params,fontsize = 10)
plotText(label = "d.10 + Diff pat.1",y=2.2,just = c("top","left"),x = 1.2,params = params,fontsize = 10)
plotText(label = "d.10 + Diff pat.2",y=1.2,just = c("top","left"),x = 1.2,params = params,fontsize = 10)
plotText(label = "H3K27ac enrichment around NFATC1",y=0.2,just = c("top","centre"),
         x = 4.5,params = params,fontsize = 10,fontface="bold")
pageGuideHide()
dev.off()



png("Z:/Robel/OC_diff_diuron_paper/figures_genecode/cd93_browserplot_chr_bw.png",height = 9,width = 10,units = "cm",res=600)
pageCreate(width = 10, height = 9, default.units = "cm",showGuides = T)
## Set the coordinates
# plotGG(plot = corr_test,x=10.2,y = 0,just = c("top","left"),width = 9.8,height = 8,default.units = "cm")
# plotGG(plot = corr_test,x=13,y = 8.5,just = c("top","left"),width = 9.5,height = 5.5,default.units = "cm")
# chr18:77,100,195-77,510,622
params <- pgParams(
  chrom = "chr20",
  chromstart = 23050000, chromend = 23154000,
  assembly = "hg19",x=1,
  width = 9,just=c('bottom',"left"),default.units="cm",range=c(0,200)
)
genes = plotGenes(params = params,geneOrder = c("CD93"),geneHighlights = data.frame(gene=c("CD93"),color=c("red")),
                  y = 8,height = 0.8,strandLabels = F)
annoGenomeLabel(fontsize = 10,params = params,y = 9,height=0.8,genes)
signal1 = plotSignal(data = "Z:/Robel/OC_diff_diuron_paper/deeptools/OC.J3.-.pat1_merged_chr.bw",
                     params = params,y=7,height = 0.6,fill = rgb(maxColorValue = 255,red = 67,green = 147,blue = 195),
                     linecolor = rgb(maxColorValue = 255,red = 67,green = 147,blue = 195)
)
# annoYaxis(signal1,fontsize=8,params = params,fontcolor="black",axisLine = T)
signal2 = plotSignal(data = "Z:/Robel/OC_diff_diuron_paper/deeptools/OC.J3.-.pat2_merged_chr.bw",
                     params = params,y=6,height = 0.6,fill = rgb(maxColorValue = 255,red = 67,green = 147,blue = 195),
                     linecolor = rgb(maxColorValue = 255,red = 67,green = 147,blue = 195)
)
# annoYaxis(signal2,fontsize=8,params = params,fontcolor="black",axisLine = T)
signal3 = plotSignal(data = "Z:/Robel/OC_diff_diuron_paper/deeptools/OC.J7.+.pat1_merged_chr.bw",
                     params = params,y=5,height = 0.6,fill = "#FF7F00",
                     linecolor = "#FF7F00"
)
# annoYaxis(signal3,fontsize=8,params = params,fontcolor="black",axisLine = T)
signal4 = plotSignal(data = "Z:/Robel/OC_diff_diuron_paper/deeptools/OC.J7.+.pat2_merged_chr.bw",
                     params = params,y=4,height = 0.6,fill = "#FF7F00",
                     linecolor = "#FF7F00"
)
# annoYaxis(signal4,fontsize=8,params = params,fontcolor="black",axisLine = T)
signal5 = plotSignal(data = "Z:/Robel/OC_diff_diuron_paper/deeptools/OC.J10.+.pat1_merged_chr.bw",
                     params = params,y=3,height = 0.6,fill = rgb(maxColorValue = 255,red = 255,green = 0,blue = 0),
                     linecolor = rgb(maxColorValue = 255,red = 255,green = 0,blue = 0)
)
# annoYaxis(signal5,fontsize=8,params = params,fontcolor="black",axisLine = T)
signal6 = plotSignal(data = "Z:/Robel/OC_diff_diuron_paper/deeptools/OC.J10.+.pat2_merged_chr.bw",
                     params = params,y=2,height = 0.6,fill = rgb(maxColorValue = 255,red = 255,green = 0,blue = 0),
                     linecolor = rgb(maxColorValue = 255,red = 255,green = 0,blue = 0)
)
annoYaxis(signal6,fontsize=8,params = params,fontcolor="black",axisLine = T)

annoHighlight(signal1,chrom = "chr20",chromstart = 23122846,chromend = 23137416 ,fill = "yellow",alpha = 0.3,
              params = params,y = 7,height = 5.6)

plotText(label = "SE",x = 8,y = 1,params = params,just = c("top","center"),fontsize = 10)
plotText(label = "d.3 - Diff pat.1",y=6.2,just = c("top","left"),x = 1.2,params = params,fontsize = 10)
plotText(label = "d.3 - Diff pat.2",y=5.2,just = c("top","left"),x = 1.2,params = params,fontsize = 10)
plotText(label = "d.7 + Diff pat.1",y=4.2,just = c("top","left"),x = 1.2,params = params,fontsize = 10)
plotText(label = "d.7 + Diff pat.2",y=3.2,just = c("top","left"),x = 1.2,params = params,fontsize = 10)
plotText(label = "d.10 + Diff pat.1",y=2.2,just = c("top","left"),x = 1.2,params = params,fontsize = 10)
plotText(label = "d.10 + Diff pat.2",y=1.2,just = c("top","left"),x = 1.2,params = params,fontsize = 10)
plotText(label = "H3K27ac enrichment around CD93",y=0.2,just = c("top","centre"),
         x = 4.5,params = params,fontsize = 10,fontface="bold")
pageGuideHide()
dev.off()


genecode_gene_cln = genecode %>% data.frame() %>% filter(type=="gene") %>% separate(gene_id,sep = "\\.",c("gene_id_clean","version"),remove = F)
pegasus=read_tsv("Z:/Bioinfo/PEGASUS_enh_promoter_synteny/hg19/hg19_CNEs_PEGASUS.data.gz")




peg_preds = apply(de_se_d10, 1, function(x){
  se = x["enhancer"] %>% GRanges() %>% `seqlevelsStyle<-`("UCSC")
  
  peg_cnes = pegasus[queryHits(findOverlaps(query = (pegasus %>% 
                                                       mutate(seqnames=chromosome) %>% 
                                                       select(c(seqnames,start,end)) %>% 
                                                       plyranges::as_granges()),subject = se)),]
  if(nrow(peg_cnes)==0){
    return(c(PEG_CNE="NO_CNE",Target_PEG = "NO_CNE"))
  }else{
    overlapping_cnes = peg_cnes %>% 
      mutate(CNE_peg = paste0(chromosome,":",start,'-',end)) %>% 
      pull(CNE_peg) %>% 
      paste(collapse = "; ")
    genes = genecode_gene_cln %>% 
      filter(gene_id_clean%in%str_split_fixed(pattern = "_",string = peg_cnes$target,n = Inf)) %>% 
      pull(gene_name) %>% 
      paste(collapse = "; ")
    
    return(c(PEG_CNE = overlapping_cnes,Target_PEG = genes))
  }
}
) %>% t()

# Clément et al, NAR, 2020
# PEG_CNE = Conserved non coding elements from Pegasus overlapping the super-enhancers identified. 
# PEG_targets =  And the gene(s) showing the highest linkage score (synteny) with the CNE, in a fixed window (1 Mb)
de_se_d10$PEG_CNE = peg_preds[,"PEG_CNE"]
de_se_d10$PEG_targets = peg_preds[,"Target_PEG"]
de_se_d10$target_symbols = str_replace_all(str_split_fixed(de_se_d10$target,pattern = "-> ",n=2)[,2],replacement = "; ",pattern = "\\|\\|")



de_se_d10$peg_n = sapply(de_se_d10$PEG_targets, function(x){
  
  if(x!="NO_CNE"){
    return(length(str_split_fixed(x, pattern = "; ", n = Inf)))
    }else{return("NO_CNE")}})
de_se_d10$local_n = sapply(de_se_d10$target_symbols, function(x){
  if(x!=""){return(length(str_split_fixed(x, pattern = "; ", n = Inf)))
    }else{return(NA)}})
de_se_d10$overlap = sapply(c(1:nrow(de_se_d10)), function(x){
  if(de_se_d10$target_symbols[x]!=""){
    if(de_se_d10$PEG_targets[x]=="NO_CNE"){
      return("NO_CNE")
    }else{
      return(paste(intersect(str_split_fixed(de_se_d10$target_symbols[x],pattern = "; ",n = Inf),str_split_fixed(de_se_d10$PEG_targets[x],pattern = "; ",n = Inf)),
                   collapse = ","))
    }
  }else{
    return(NA)
  }
  })
de_se_d10$overlap_n = sapply(de_se_d10$overlap, function(x){
  if(!is.na(x) & x!="NO_CNE"){
    overlap_x = str_split_fixed(x,pattern = ",",n=Inf)
    if(overlap_x==""){
      return(NA)
    }else{
        length(overlap_x)
      }
    }else{return(NA)}})
de_se_d10$overlap_n[which(de_se_d10$PEG_CNE=="NO_CNE")] = "NO_CNE"

de_se_d10$target_local_n = as.numeric(sapply(de_se_d10$target_symbols, function(x){length(unlist(strsplit(x,split="; ")))}))
# de_se_d10$overlap_perc[which(de_se_d10$target_local_n>0 & !is.na(de_se_d10$overlap_n))] = de_se_d10$overlap_n[which(de_se_d10$target_local_n>0 & !is.na(de_se_d10$overlap_n))]/de_se_d10$target_local_n[which(de_se_d10$target_local_n>0 & !is.na(de_se_d10$overlap_n))]

de_se_d10_overlap_filt = de_se_d10 %>% filter(!is.na(overlap_n) & overlap_n!="NO_CNE") %>% mutate(overlap_perc = as.numeric(overlap_n)/target_local_n*100)

#######
# DAY 10 enhancers
#######
peg_preds = apply(de_enh_d10, 1, function(x){
  se = x["enhancer"] %>% GRanges() %>% `seqlevelsStyle<-`("UCSC")
  
  peg_cnes = pegasus[queryHits(findOverlaps(query = (pegasus %>% 
                                                       mutate(seqnames=chromosome) %>% 
                                                       select(c(seqnames,start,end)) %>% 
                                                       plyranges::as_granges()),subject = se)),]
  if(nrow(peg_cnes)==0){
    return(c(PEG_CNE="NO_CNE",Target_PEG = "NO_CNE"))
  }else{
    overlapping_cnes = peg_cnes %>% 
      mutate(CNE_peg = paste0(chromosome,":",start,'-',end)) %>% 
      pull(CNE_peg) %>% 
      paste(collapse = "; ")
    genes = genecode_gene_cln %>% 
      filter(gene_id_clean%in%str_split_fixed(pattern = "_",string = peg_cnes$target,n = Inf)) %>% 
      pull(gene_name) %>% 
      paste(collapse = "; ")
    
    return(c(PEG_CNE = overlapping_cnes,Target_PEG = genes))
  }
}
) %>% t()

# Clément et al, NAR, 2020
# PEG_CNE = Conserved non coding elements from Pegasus overlapping the super-enhancers identified. 
# PEG_targets =  And the gene(s) showing the highest linkage score (synteny) with the CNE, in a fixed window (1 Mb)
de_enh_d10$PEG_CNE = peg_preds[,"PEG_CNE"]
de_enh_d10$PEG_targets = peg_preds[,"Target_PEG"]
de_enh_d10$target_symbols = str_replace_all(str_split_fixed(de_enh_d10$target,pattern = "-> ",n=2)[,2],replacement = "; ",pattern = "\\|\\|")



de_enh_d10$peg_n = sapply(de_enh_d10$PEG_targets, function(x){
  
  if(x!="NO_CNE"){
    return(length(str_split_fixed(x, pattern = "; ", n = Inf)))
  }else{return("NO_CNE")}})
de_enh_d10$local_n = sapply(de_enh_d10$target_symbols, function(x){
  if(x!=""){return(length(str_split_fixed(x, pattern = "; ", n = Inf)))
  }else{return(NA)}})
de_enh_d10$overlap = sapply(c(1:nrow(de_enh_d10)), function(x){
  if(de_enh_d10$target_symbols[x]!=""){
    if(de_enh_d10$PEG_targets[x]=="NO_CNE"){
      return("NO_CNE")
    }else{
      return(paste(intersect(str_split_fixed(de_enh_d10$target_symbols[x],pattern = "; ",n = Inf),str_split_fixed(de_enh_d10$PEG_targets[x],pattern = "; ",n = Inf)),
                   collapse = ","))
    }
  }else{
    return(NA)
  }
})
de_enh_d10$overlap_n = sapply(de_enh_d10$overlap, function(x){
  if(!is.na(x) & x!="NO_CNE"){
    overlap_x = str_split_fixed(x,pattern = ",",n=Inf)
    if(overlap_x==""){
      return(NA)
    }else{
      length(overlap_x)
    }
  }else{return(NA)}})
de_enh_d10$overlap_n[which(de_enh_d10$PEG_CNE=="NO_CNE")] = "NO_CNE"

de_enh_d10$target_local_n = as.numeric(sapply(de_enh_d10$target_symbols, function(x){length(unlist(strsplit(x,split="; ")))}))
# de_se_d10$overlap_perc[which(de_se_d10$target_local_n>0 & !is.na(de_se_d10$overlap_n))] = de_se_d10$overlap_n[which(de_se_d10$target_local_n>0 & !is.na(de_se_d10$overlap_n))]/de_se_d10$target_local_n[which(de_se_d10$target_local_n>0 & !is.na(de_se_d10$overlap_n))]

de_enh_d10_overlap_filt = de_enh_d10 %>% filter(!is.na(overlap_n) & overlap_n!="NO_CNE") %>% mutate(overlap_perc = as.numeric(overlap_n)/target_local_n*100)







#######
# DAY 7
#######


peg_preds = apply(de_se_d7, 1, function(x){
  se = x["enhancer"] %>% GRanges() %>% `seqlevelsStyle<-`("UCSC")
  
  peg_cnes = pegasus[queryHits(findOverlaps(query = (pegasus %>% 
                                                       mutate(seqnames=chromosome) %>% 
                                                       select(c(seqnames,start,end)) %>% 
                                                       plyranges::as_granges()),subject = se)),]
  if(nrow(peg_cnes)==0){
    return(c(PEG_CNE="NO_CNE",Target_PEG = "NO_CNE"))
  }else{
    overlapping_cnes = peg_cnes %>% 
      mutate(CNE_peg = paste0(chromosome,":",start,'-',end)) %>% 
      pull(CNE_peg) %>% 
      paste(collapse = "; ")
    genes = genecode_gene_cln %>% 
      filter(gene_id_clean%in%str_split_fixed(pattern = "_",string = peg_cnes$target,n = Inf)) %>% 
      pull(gene_name) %>% 
      paste(collapse = "; ")
    
    return(c(PEG_CNE = overlapping_cnes,Target_PEG = genes))
  }
}
) %>% t()

de_se_d7$PEG_CNE = peg_preds[,"PEG_CNE"]
de_se_d7$PEG_targets = peg_preds[,"Target_PEG"]
de_se_d7$target_symbols = str_replace_all(str_split_fixed(de_se_d7$target,pattern = "-> ",n=2)[,2],replacement = "; ",pattern = "\\|\\|")



de_se_d7$peg_n = sapply(de_se_d7$PEG_targets, function(x){
  
  if(x!="NO_CNE"){
    return(length(str_split_fixed(x, pattern = "; ", n = Inf)))
  }else{return("NO_CNE")}})

de_se_d7$local_n = sapply(de_se_d7$target_symbols, function(x){
  if(x!=""){return(length(str_split_fixed(x, pattern = "; ", n = Inf)))
  }else{return(NA)}})

de_se_d7$overlap = sapply(c(1:nrow(de_se_d7)), function(x){
  if(de_se_d7$target_symbols[x]!=""){
    if(de_se_d7$PEG_targets[x]=="NO_CNE"){
      return("NO_CNE")
    }else{
      return(paste(intersect(str_split_fixed(de_se_d7$target_symbols[x],pattern = "; ",n = Inf),str_split_fixed(de_se_d7$PEG_targets[x],pattern = "; ",n = Inf)),
                   collapse = ","))
    }
  }else{
    return(NA)
  }
})
de_se_d7$overlap_n = sapply(de_se_d7$overlap, function(x){
  if(!is.na(x) & x!="NO_CNE"){
    overlap_x = str_split_fixed(x,pattern = ",",n=Inf)
    if(overlap_x==""){
      return(NA)
    }else{
      length(overlap_x)
    }
  }else{return(NA)}})
de_se_d7$overlap_n[which(de_se_d7$PEG_CNE=="NO_CNE")] = "NO_CNE"

de_se_d7$target_local_n = as.numeric(sapply(de_se_d7$target_symbols, function(x){length(unlist(strsplit(x,split="; ")))}))
# de_se_d7$overlap_perc[which(de_se_d7$target_local_n>0 & !is.na(de_se_d7$overlap_n))] = 
#   de_se_d7$overlap_n[which(de_se_d7$target_local_n>0 & !is.na(de_se_d7$overlap_n))]/
#   de_se_d7$target_local_n[which(de_se_d7$target_local_n>0 & !is.na(de_se_d7$overlap_n))]

de_se_d7_overlap_filt = de_se_d7 %>% filter(!is.na(overlap_n) & overlap_n!="NO_CNE") %>% mutate(overlap_perc = as.numeric(overlap_n)/target_local_n*100)



#####
#day 7 enhancers
#####
peg_preds = apply(de_enh_d7, 1, function(x){
  se = x["enhancer"] %>% GRanges() %>% `seqlevelsStyle<-`("UCSC")
  
  peg_cnes = pegasus[queryHits(findOverlaps(query = (pegasus %>% 
                                                       mutate(seqnames=chromosome) %>% 
                                                       select(c(seqnames,start,end)) %>% 
                                                       plyranges::as_granges()),subject = se)),]
  if(nrow(peg_cnes)==0){
    return(c(PEG_CNE="NO_CNE",Target_PEG = "NO_CNE"))
  }else{
    overlapping_cnes = peg_cnes %>% 
      mutate(CNE_peg = paste0(chromosome,":",start,'-',end)) %>% 
      pull(CNE_peg) %>% 
      paste(collapse = "; ")
    genes = genecode_gene_cln %>% 
      filter(gene_id_clean%in%str_split_fixed(pattern = "_",string = peg_cnes$target,n = Inf)) %>% 
      pull(gene_name) %>% 
      paste(collapse = "; ")
    
    return(c(PEG_CNE = overlapping_cnes,Target_PEG = genes))
  }
}
) %>% t()

de_enh_d7$PEG_CNE = peg_preds[,"PEG_CNE"]
de_enh_d7$PEG_targets = peg_preds[,"Target_PEG"]
de_enh_d7$target_symbols = str_replace_all(str_split_fixed(de_enh_d7$target,pattern = "-> ",n=2)[,2],replacement = "; ",pattern = "\\|\\|")



de_enh_d7$peg_n = sapply(de_enh_d7$PEG_targets, function(x){
  
  if(x!="NO_CNE"){
    return(length(str_split_fixed(x, pattern = "; ", n = Inf)))
  }else{return("NO_CNE")}})

de_enh_d7$local_n = sapply(de_enh_d7$target_symbols, function(x){
  if(x!=""){return(length(str_split_fixed(x, pattern = "; ", n = Inf)))
  }else{return(NA)}})

de_enh_d7$overlap = sapply(c(1:nrow(de_enh_d7)), function(x){
  if(de_enh_d7$target_symbols[x]!=""){
    if(de_enh_d7$PEG_targets[x]=="NO_CNE"){
      return("NO_CNE")
    }else{
      return(paste(intersect(str_split_fixed(de_enh_d7$target_symbols[x],pattern = "; ",n = Inf),str_split_fixed(de_enh_d7$PEG_targets[x],pattern = "; ",n = Inf)),
                   collapse = ","))
    }
  }else{
    return(NA)
  }
})
de_enh_d7$overlap_n = sapply(de_enh_d7$overlap, function(x){
  if(!is.na(x) & x!="NO_CNE"){
    overlap_x = str_split_fixed(x,pattern = ",",n=Inf)
    if(overlap_x==""){
      return(NA)
    }else{
      length(overlap_x)
    }
  }else{return(NA)}})
de_enh_d7$overlap_n[which(de_enh_d7$PEG_CNE=="NO_CNE")] = "NO_CNE"

de_enh_d7$target_local_n = as.numeric(sapply(de_enh_d7$target_symbols, function(x){length(unlist(strsplit(x,split="; ")))}))
# de_se_d7$overlap_perc[which(de_se_d7$target_local_n>0 & !is.na(de_se_d7$overlap_n))] = 
#   de_se_d7$overlap_n[which(de_se_d7$target_local_n>0 & !is.na(de_se_d7$overlap_n))]/
#   de_se_d7$target_local_n[which(de_se_d7$target_local_n>0 & !is.na(de_se_d7$overlap_n))]

de_enh_d7_overlap_filt = de_enh_d7 %>% filter(!is.na(overlap_n) & overlap_n!="NO_CNE") %>% mutate(overlap_perc = as.numeric(overlap_n)/target_local_n*100)

