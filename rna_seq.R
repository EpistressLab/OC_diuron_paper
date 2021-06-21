#####################################################################################################
###### Setting up necessary elements ################################################################
#####################################################################################################

source('toGithub/functions.R')
# List of genes of interest in osteoclastogenesis
genesOfInterest=read.table('genesOfInterest.txt',stringsAsFactors = F)
genesOfInterest=unique(genesOfInterest$V1)
# Checking gene names are official gene symbols.
genesOfInterest=limma::alias2SymbolTable(genesOfInterest)
mostInterestingGenes=c('C20orf123','CTSK','SIGLEC15','MMP9','NFATC1','JDP2','RAB38','ACP5')

mostInterestingGenes=limma::alias2SymbolTable(mostInterestingGenes)

# importing raw count matrices and building colData
diuron_gene=as.matrix(read.csv('diuron_gene_count_matrix.csv',row.names = 'gene_id'))
newColData=data.frame(id=colnames(diuron_gene),
                      day=factor(c('10','10','10','10','10','10','3','3','3','3','3','3'),levels=c('3','10')),
                      treatment=factor(c('DMSO','DMSO','DMSO','Diuron','Diuron','Diuron','DMSO','DMSO','DMSO','Diuron','Diuron','Diuron'),
                                       levels = c('DMSO','Diuron')),
                      row.names = 'id')

# One condition column for some complicated comparisons
newColData$condition =factor(paste0(newColData$treatment,"_",newColData$day))
newColData$condition = relevel(newColData$condition,ref = "Diuron_10")

# DDS object creation with full design and LRT test
newDesign_diuron_gene=DESeqDataSetFromMatrix(countData,colData = newColData,design = ~ treatment + day+ treatment:day)
newDesign_diuron_gene=DESeq(newDesign_diuron_gene,reduced = ~treatment + day,test = 'LRT')

# DDS for potentially complicated pair-wise comparisons
dds_combined = DESeqDataSetFromMatrix(countData,colData = newColData,design = ~ condition)
dds_combined = DESeq(dds_combined)

rld = rlogTransformation(newDesign_diuron_gene,blind = F)
mrld = assay(rld)
plotPCAWithSampleNames(rld,intgroup=c('day','treatment'))


# Gene names of count data are ENSG000|SYMBOL so list of names of genes of interest 
# has to correspond to what's found within the DESeq2 object
genesOfInterest_ext = rownames(newDesign_diuron_gene)[
  which(grepl(pattern = paste0(paste0("\\|",genesOfInterest,"\\b$"),collapse='|'),
              rownames(newDesign_diuron_gene)))]
mostInterestingGenes_ext = rownames(newDesign_diuron_gene)[
  which(grepl(pattern = paste0(paste0("\\|",mostInterestingGenes,"\\b$"),collapse='|'),
              rownames(newDesign_diuron_gene)))]

# treatment_Diuron_vs_DMSO effect using LRT
resultsLRT=results(newDesign_diuron_gene)%>%
  data.frame()%>%
  rownames_to_column(var='SYMBNAME')%>%
  filter(padj < 0.05 )

# Columns annotation and coloring for heatmap
annot = data.frame(colData(newDesign_diuron_gene)[,c('day','treatment')],stringsAsFactors = T)
annot=data.frame(annot,stringsAsFactors = T)
colTopAnnot<-vector("list", ncol(annot))
names(colTopAnnot)<-colnames(annot)

colTopAnnot$day[["3"]] = "#e5f5e0"
colTopAnnot$day[["10"]] = "#31a354"
colTopAnnot$treatment[["DMSO"]] = "#4393c3"
colTopAnnot$treatment[["Diuron"]] = "#d6604d"

# Preferential Choice of samples order submission, eventhough clustering will take the upper hand
samplesOrder = c("OC_J3_DMSO_rep1","OC_J3_DMSO_rep2","OC_J3_DMSO_rep3",
                 "OC_J10_DMSO_rep1","OC_J10_DMSO_rep2","OC_J10_DMSO_rep3",
                 "OC_J3_Diuron_rep1","OC_J3_Diuron_rep2","OC_J3_Diuron_rep3",
                 "OC_J10_Diuron_rep1","OC_J10_Diuron_rep2","OC_J10_Diuron_rep3")

# Heatmap of the effect of diuron treatment on OC differentiation (LRT test)
ht1=plotggplot_Heatmap(newDesign_diuron_gene,annot = annot,
                       DE = resultsLRT$SYMBNAME,
                       genesOfInterest = genesOfInterest,
                       condCol = c('day','treatment'), 
                       name = '(treatment effect over time)',
                       comp = 'fullModel',
                       colTopAnnot = colTopAnnot,
                       samplesOrder = samplesOrder,
                       split = 2)

# Gene clusters with similar dynamic among DE genes
clustering_sig_genes_LRT=resultsLRT%>%arrange(padj)
cluster_rlog=mrld[clustering_sig_genes_LRT$SYMBNAME,]
clusters <- degPatterns(cluster_rlog, metadata = newColData, time="day", col="treatment",plot=F)
degPlotCluster(clusters[["normalized"]],time = "day",color = "treatment",lines = F,textSize = 15)

# cluster groups of genes of interest
table(clusters$df[which(grepl(pattern = paste0(paste0("\\.",genesOfInterest,"\\b"),collapse='|'),clusters$df$genes)),])

# Individual counts of genes of interest
degPlotLocal(dds = newDesign_diuron_gene,
             genes =genesOfInterest_ext,
             xs = 'day',n=40,group = 'treatment',xsLab = "time",
             ysLab = "normalized counts \n(Median of ratios method)",
             splitGeneName = 2,
             reverseColor = T)


##### Functional Annotation #####

# getting ENTREZID for DE genes
geneList=bitr((sapply(resultsLRT$SYMBNAME, function(x){unlist(strsplit(x,split = '|',fixed = T))[2]})),
              fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)[,'ENTREZID']
# GO enrichment
resultsLRT_GO=enrichGO(gene = geneList,OrgDb = org.Hs.eg.db,ont = 'BP',readable = T,maxGSSize = 100)
# Removing redundent GO terms
resultsLRT_GO=simplify(resultsLRT_GO)
# subset=c("GO:0022617","GO:0090382","GO:0048016","GO:0033173","GO:0006029","GO:0030574","GO:0007045","GO:0061718","GO:0051453","GO:000673","GO:1904645")

plotEnrichFactor(resultsLRT_GO,title = 'Enriched GO BP gene sets',n=30)
# # plotEnrichFactor(resultsLRT_GO[subset],title = 'Enriched GO BP gene sets',n=30)

# Reactome.db enrichment
resultsLRT_Pathway=enrichPathway(gene=geneList,organism='human')
# # unlist(mget(get('GO:0022617',org.Hs.egGO2ALLEGS),org.Hs.egSYMBOL))[which(unlist(mget(get('GO:0022617',org.Hs.egGO2ALLEGS),org.Hs.egSYMBOL))%in%limma::alias2Symbol(resultsLRT$SYMBNAME))]
plotEnrichFactor(resultsLRT_Pathway,title = 'Enriched Reactome db gene sets',n=30)

################################################################
################################################################
##### Pair-wise comparisons #####
################################################################
################################################################


################################################################
##### day 10 DMSO vs day 10 Diuron ######

results = lfcShrink(dds = dds_combined,res = results(dds_combined,name = "condition_DMSO_10_vs_Diuron_10",test = "Wald"),coef = "condition_DMSO_10_vs_Diuron_10",type = "apeglm") %>% data.frame() %>% rownames_to_column("SYMBNAME") %>% filter(abs(log2FoldChange) > 1 & padj < 0.05)
samplesOrder = c("OC_J3_DMSO_rep1","OC_J3_DMSO_rep2","OC_J3_DMSO_rep3","OC_J3_Diuron_rep1","OC_J3_Diuron_rep2","OC_J3_Diuron_rep3","OC_J10_DMSO_rep1","OC_J10_DMSO_rep2","OC_J10_DMSO_rep3","OC_J10_Diuron_rep1","OC_J10_Diuron_rep2","OC_J10_Diuron_rep3")
ht1=plotggplot_Heatmap(newDesign_diuron_gene,
                       DE = results$SYMBNAME,
                       genesOfInterest = genesOfInterest,
                       condCol = c('day','treatment'), 
                       name = '(treatment effect at day 10)',
                       comp = 'fullModel',colTopAnnot = colTopAnnot,
                       samplesOrder = samplesOrder,split = 2)
draw(ht1)

geneList=bitr((sapply(results$SYMBNAME, function(x){unlist(strsplit(x,split = '|',fixed = T))[2]})),
              fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)[,'ENTREZID']
results_GO=enrichGO(gene = geneList,OrgDb = org.Hs.eg.db,ont = 'BP',readable = T)
results_GO=simplify(results_GO)
plotEnrichFactor(results_GO,title = 'Enriched GO BP gene sets',n=30)

GO_Interest=c('GO:0045453','GO:0046849','GO:1901136','GO:0061041','GO:0030316','GO:0003414',
              'GO:0045124','GO:0003416',"GO:0006027","GO:0032963","GO:0032623","GO:0050702",
              "GO:0071222","GO:0043030","GO:0032088","GO:0007249","GO:0060348","GO:0048771","GO:0030198")

# extracting log2 fold-changes for DE genes 
upFchanges=results%>%
  filter(abs(log2FoldChange) >= 1 & padj < 0.05)%>%
  pull(log2FoldChange)

# naming the log2 fold changes with gene names
names(upFchanges)=results%>%
  filter(abs(log2FoldChange) >= 1 & padj < 0.05)%>%
  pull(SYMBNAME) %>% 
  sapply(function(x){unlist(strsplit(x,split = '|',fixed = T))[2]})

# interesting GO: terms among significantly enriched are listed
GO_Interest_2=c('GO:0030198','GO:0006027','GO:0048771','GO:1901136','GO:0032963','GO:0032623','GO:0060348','GO:0046849')
ego=results_GO
ego@result=results_GO@result[which(results_GO@result$ID%in%GO_Interest_2),]
# cairo_ps(filename = "figures_genecode/Du_10DMvs10Du_cnetPlot_test.ps",fallback_resolution = 600,width = 9)
# cnetplot with GO: terms genes and their log2 fold changes
cnetplot(ego,categorySize = 'p.adjust',showCategory=9,foldChange = upFchanges)+
  theme(text=element_text(size=12))
# dev.off()

# GO_Interest=c('GO:0045453','GO:0046849','GO:1901136','GO:0061041','GO:0030316','GO:0003414','GO:0045124','GO:0003416',"GO:0006027","GO:0032963")
# upFchanges=results%>%filter(abs(log2FoldChange) >= 1 & padj < 0.05)%>%pull(log2FoldChange)
# names(upFchanges)=results%>%filter(abs(log2FoldChange) >= 1 & padj < 0.05)%>%pull(SYMBNAME) %>% sapply(function(x){unlist(strsplit(x,split = '|',fixed = T))[2]})
# ego=results_GO
# ego@result=results_GO@result[which(results_GO@result$ID%in%GO_Interest),]
# 
# cnetplot(ego,categorySize = 'p.adjust',showCategory=12,foldChange = upFchanges)

#####################################################################################################
##### day 10 DMSO vs day 3 DMSO ######
#####################################################################################################
results = lfcShrink(dds = newDesign_diuron_gene,res = results(newDesign_diuron_gene,name = "day_10_vs_3",test="Wald"),coef = "day_10_vs_3",type = "apeglm") %>% data.frame()%>%rownames_to_column(var='SYMBNAME')%>%filter(abs(log2FoldChange) > 1 & padj < 0.05)


samplesOrder = c("OC_J3_DMSO_rep1","OC_J3_DMSO_rep2","OC_J3_DMSO_rep3","OC_J3_Diuron_rep1","OC_J3_Diuron_rep2","OC_J3_Diuron_rep3","OC_J10_DMSO_rep1","OC_J10_DMSO_rep2","OC_J10_DMSO_rep3","OC_J10_Diuron_rep1","OC_J10_Diuron_rep2","OC_J10_Diuron_rep3")
ht1=plotggplot_Heatmap(newDesign_diuron_gene,DE = results$SYMBNAME,genesOfInterest = genesOfInterest,condCol = c('day','treatment'), name = '(treatment effect at day 10)',comp = 'fullModel',colTopAnnot = colTopAnnot,samplesOrder = samplesOrder,split = 2)
# cairo_ps(filename = "figures_genecode/Du_10DMvs3DM_heatmap.ps",fallback_resolution = 600)
draw(ht1)
# dev.off()

geneList=bitr((sapply(results$SYMBNAME, function(x){unlist(strsplit(x,split = '|',fixed = T))[2]})),fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)[,'ENTREZID']
results_GO=enrichGO(gene = geneList,OrgDb = org.Hs.eg.db,ont = 'BP',readable = T)
results_GO=simplify(results_GO)
plotEnrichFactor(results_GO,title = 'Enriched GO BP gene sets',n=30)

GO_Interest=c('GO:0045453','GO:0046849','GO:1901136','GO:0061041','GO:0030316','GO:0003414','GO:0045124','GO:0003416',"GO:0006027","GO:0032963","GO:0032623","GO:0050702","GO:0071222","GO:0043030","GO:0032088","GO:0007249","GO:0060348")
upFchanges=results%>%
  filter(abs(log2FoldChange) >= 1 & padj < 0.05)%>%
  pull(log2FoldChange)
names(upFchanges)=results%>%
  filter(abs(log2FoldChange) >= 1 & padj < 0.05)%>%
  pull(SYMBNAME) %>% 
  sapply(function(x){unlist(strsplit(x,split = '|',fixed = T))[2]})
ego=results_GO
ego@result=results_GO@result[which(results_GO@result$ID%in%GO_Interest),]
cairo_ps(filename = "figures_genecode/Du_10DMvs10Du_cnetPlot.ps",fallback_resolution = 600,width = 9)
cnetplot(ego,categorySize = 'p.adjust',showCategory=9,foldChange = upFchanges)
dev.off()

results_up = results%>% 
  filter(log2FoldChange > 1 & padj < 0.05)
geneList=bitr((sapply(results_up$SYMBNAME, function(x){unlist(strsplit(x,split = '|',fixed = T))[2]})),fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)[,'ENTREZID']
results_GO=enrichGO(gene = geneList,OrgDb = org.Hs.eg.db,ont = 'BP',readable = T)
results_GO=simplify(results_GO)
cairo_ps(filename = "figures_genecode/Du_10DMvs10Du_up_heatmap.ps",fallback_resolution = 600)
plotEnrichFactor(results_GO,title = 'Enriched GO BP gene sets',n=30)
dev.off()

results_down = results%>% 
  filter(log2FoldChange < -1 & padj < 0.05)
geneList=bitr((sapply(results_down$SYMBNAME, function(x){unlist(strsplit(x,split = '|',fixed = T))[2]})),fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)[,'ENTREZID']
results_GO=enrichGO(gene = geneList,OrgDb = org.Hs.eg.db,ont = 'BP',readable = T)
results_GO=simplify(results_GO)
cairo_ps(filename = "figures_genecode/Du_10DMvs10Du_up_heatmap.ps",fallback_resolution = 600)
plotEnrichFactor(results_GO,title = 'Enriched GO BP gene sets',n=30)
dev.off()

GO_Interest=c('GO:0045453','GO:0046849','GO:1901136','GO:0061041','GO:0030316','GO:0003414','GO:0045124','GO:0003416',"GO:0006027","GO:0032963")
upFchanges=results%>%filter(abs(log2FoldChange) >= 1 & padj < 0.05)%>%pull(log2FoldChange)
names(upFchanges)=results%>%filter(abs(log2FoldChange) >= 1 & padj < 0.05)%>%pull(SYMBNAME) %>% sapply(function(x){unlist(strsplit(x,split = '|',fixed = T))[2]})
ego=results_GO
ego@result=results_GO@result[which(results_GO@result$ID%in%GO_Interest),]

cnetplot(ego,categorySize = 'p.adjust',showCategory=12,foldChange = upFchanges)

# GSEA
res = lfcShrink(dds = newDesign_diuron_gene,res = results(newDesign_diuron_gene,name = "day_10_vs_3",test="Wald"),
                coef = "day_10_vs_3",type = "apeglm") %>% 
  data.frame()%>%rownames_to_column(var='SYMBNAME')

dataTot=res
dataTot$ENSID = sapply(dataTot$SYMBNAME, function(x){unlist(strsplit(x,split = '|',fixed = T))[1]})
dataTot$entrezID = idCorrTable$entrezgene_id[match(dataTot$ENSID,idCorrTable$ENSID)]
dataTot=dataTot[-c(which(is.na(dataTot$padj))),]
geneID=clusterProfiler::bitr(dataTot$SYMBNAME, fromType="SYMBOL",toType=c("ENTREZID"),OrgDb=org.Hs.eg.db)
doublons=which(duplicated(geneID$SYMBOL))
if(length(doublons)>0){geneID=geneID[-doublons,]}
dataTot=dataTot[is.element(dataTot%>%pull(SYMBNAME), geneID$SYMBOL),]

geneListTot=dataTot$log2FoldChange
names(geneListTot)=dataTot$entrezID
# Sorting gene list with log2FoldChange
geneListTot=sort(geneListTot,decreasing = T)

# GSEA on GO:terms

ego2=gseGO(geneList=geneListTot,OrgDb=org.Hs.eg.db,ont="BP",nPerm=1000,minGSSize=20,maxGSSize=500,pvalueCutoff=0.05,verbose=FALSE,keyType = "ENTREZID")
# ego2=simplify(ego2)
ego2Symbol=setReadable(ego2,org.Hs.eg.db,keyType = "ENTREZID")
gseInterest=c('GO:0060350','GO:0001501','GO:0061448','GO:0060348','GO:1901136','GO:0051216','GO:0032963','GO:0019886','GO:0060349','GO:0060351','GO:0030198','GO:0038095','GO:0002474','GO:0036230','GO:0098868','GO:0003416','GO:0002431','GO:0050766','GO:0042476','GO:0010743','GO:0001501')
#   c('GO:0038095','GO:0002474','GO:0036230','GO:0098868','GO:0003416','GO:0002431','GO:0050766')
# gseUp=c('GO:0042476','GO:0010743','GO:0001501')
ego=ego2
ego@result=ego2@result[which(ego2@result$ID%in%gseInterest),]
# kable(ego2Symbol%>%data.frame(),caption = 'GSEA GO gene sets in (D8_DMSO Vs D1_DMSO)')%>%kable_styling(full_width = F)%>%
#   scroll_box(width = "100%", height = "200px")
write.csv(x=ego2Symbol%>%data.frame(),file = 'gsea_GO_D8DMSO_vs_D1DMSO.csv')
ridgeplot(ego)
gseaplot2(ego2,geneSetID = gseInterest)

# GSEA on Reactome.db

ePath=ReactomePA::gsePathway(geneList = geneListTot,organism = 'human',verbose = F,nPerm = 1000,minGSSize = 20,maxGSSize = 500,pvalueCutoff = 0.1)
epath_up=c('R-HSA-1442490','R-HSA-1474290','R-HSA-1793185','R-HSA-1474228','R-HSA-983231','R-HSA-1266695','R-HSA-1592389')
epath_down=c('R-HSA-1236978','R-HSA-2871837','R-HSA-168253','R-HSA-909733','R-HSA-6783783','R-HSA-447115','HSA-1169091','R-HSA-2454202','R-HSA-5689901','R-HSA-8939902','R-HSA-5686938')
# ePath=setReadable(ePath,OrgDb = org.Hs.eg.db)
eptha_sub=ePath
eptha_sub@result=ePath@result[which(ePath@result$ID%in%c(epath_down,epath_up)),]
cairo_ps(filename = "figures_genecode/Du_10DMvs3Dm_gseaReact.ps",fallback_resolution = 600,width = 10)
ridgeplot(eptha_sub)
dev.off()

#####################################################################################################
#####################################################################################################
########## Diuron data set ##################
#####################################################################################################
#####################################################################################################

OC_diff_gene=as.matrix(read.csv('diffDataSet_gene_RAW_count_matrix.csv',row.names = 'gene_id'))
colData=data.frame(sample=colnames(OC_diff_gene),
                   treatment=rep(c('D10_Diff','D10_noDiff','D3_noDiff','D7_Diff','D7_noDiff'),
                                 each=3),
                   row.names = 'sample')
colData$treatment = factor(colData$treatment,levels = c("D3_noDiff","D7_noDiff","D10_noDiff","D7_Diff","D10_Diff"))
newColData=data.frame(row.names=colnames(OC_diff_gene),
                      day=rep(c('10','10','3','7','7'),
                              each=3),
                      treatment=rep(c('Diff','noDiff','noDiff','Diff','noDiff'),each=3))

countData=OC_diff_gene[,rownames(colData)]

arrange(newColData,day)

newColData$treatment=factor(newColData$treatment,levels = c("noDiff","Diff"))
newColData$treatment = relevel(newColData$treatment,ref = "noDiff")
newColData$day=factor(newColData$day,levels = c("3","7","10"))
newColData$day = relevel(newColData$day,ref = "3")
mm_full = model.matrix(~ day + treatment : day,newColData)
all.zero <- apply(mm_full, 2, function(x) all(x==0))
mm_full <- mm_full[,-which(all.zero)]
mm_reduced = model.matrix(~ treatment + day,newColData) # only the non combinatorial, so the nonDiffs
dds_diff_full = DESeqDataSetFromMatrix(countData = countData,colData = newColData,design = ~ treatment + day)

dds_diff_full=DESeq(dds_diff_full,full = mm_full,reduced = mm_reduced,test = 'LRT',betaPrior = F)

# a second DDS object with a simplified design
dds_diff_simp = DESeqDataSetFromMatrix(countData = countData,colData,design = ~treatment)
dds_diff_simp = DESeq(dds_diff_simp)
# countsNorm=counts(dds_diff_full,normalized=T)
rld=rlogTransformation(dds_diff_full,blind = F)
mrld = assay(rld)
# cairo_ps(filename = "figures_genecode/diff_PCA.ps",fallback_resolution = 600)
plotPCAWithSampleNames(dds_diff_full,intgroup = c("day","treatment"))

# Setting of colors for heatmap
annot = data.frame(colData(dds_diff_full)[,c('day','treatment')],stringsAsFactors = T)
annot=data.frame(annot,stringsAsFactors = T)
colTopAnnot<-vector("list", ncol(annot))
names(colTopAnnot)<-colnames(annot)

colTopAnnot$day[["3"]] = "#e5f5e0"
colTopAnnot$day[["7"]] = "#74c476"
colTopAnnot$day[["10"]] = "#238b45"
colTopAnnot$treatment[["noDiff"]] = "#4393c3"
colTopAnnot$treatment[["Diff"]] = "#d6604d"


#####################################################################################################
########## LRT on full design ##################

#
# samplesOrder = c("OC_J3_min_rep1","OC_J3_min_rep2","OC_J3_min_rep3","OC_J7_min_rep1","OC_J7_min_rep2","OC_J7_min_rep3","OC_J10_min_rep1","OC_J10_min_rep2","OC_J10_min_rep3","OC_J7_diff_rep1","OC_J7_diff_rep2","OC_J7_diff_rep3","OC_J10_diff_rep1","OC_J10_diff_rep2","OC_J10_diff_rep3")
# 
# results = results(dds_diff_full) %>% data.frame() %>% rownames_to_column("SYMBNAME") %>% filter(padj < 0.001)
# ht1=plotggplot_Heatmap(dds_diff_full,DE = results$SYMBNAME,genesOfInterest = genesOfInterest,condCol = c('day','treatment'), name = '(treatment effect at day 10)',comp = 'fullModel',colTopAnnot = colTopAnnot,samplesOrder = samplesOrder,split = 2)
# # cairo_ps(filename = "figures_genecode/diff_LRT0.001_heatmap.ps",fallback_resolution = 600)
# draw(ht1)
# dev.off()


#####################################################################################################
########## Pair-wise: day7 treatment effect ##################

# 2 step pair-wise comparison, Genes DE in both comparisons:
# - day 7 treated vs day 7 untreated cells, and 
# - day 7 treated vs day3 untreated
results = lfcShrink(dds = dds_diff_full,
                    res = results(dds_diff_full,name = "day7.treatmentDiff",test="Wald"),type = 'ashr',quiet = T) %>% 
  data.frame() %>% rownames_to_column("SYMBNAME") %>% 
  filter(abs(log2FoldChange)>1 & padj<0.05 & SYMBNAME %in% (lfcShrink(dds = dds_diff_full,
                                                                      res = results(dds_diff_full,
                                                                                    contrast = list(c("day7","day7.treatmentDiff")),
                                                                                    test="Wald"),type = 'ashr',quiet = T) %>% 
                                                              data.frame() %>% 
                                                              rownames_to_column("SYMBNAME") %>% 
                                                              filter(abs(log2FoldChange)>1 & padj<0.05) %>% 
                                                              pull(SYMBNAME)))

ht1=plotggplot_Heatmap(dds_diff_full,DE = results$SYMBNAME,genesOfInterest = genesOfInterest,condCol = c('day','treatment'), name = '(treatment effect at day 10)',comp = 'fullModel',colTopAnnot = colTopAnnot,samplesOrder = samplesOrder,split = 2)
# cairo_ps(filename = "figures_genecode/diff_2step_day7_heatmap.ps",fallback_resolution = 600)
draw(ht1)
dev.off()

geneList=bitr((sapply(results$SYMBNAME, function(x){unlist(strsplit(x,split = '|',fixed = T))[2]})),fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)[,'ENTREZID']
results_GO=enrichGO(gene = geneList,OrgDb = org.Hs.eg.db,ont = 'BP',readable = T)
# results_GO=simplify(results_GO)
plotEnrichFactor(results_GO,title = 'Enriched GO BP gene sets',n=30)
GO_Interest=c("GO:0006029","GO:0050746","GO:0006026","GO:0030641","GO:0042475","GO:0046849","GO:0014065","GO:0006027",
              "GO:0046850","GO:0051453","GO:0090382","GO:0051452","GO:0090383","GO:0014031","GO:0048771",
              "GO:0042476","GO:0070371","GO:0034103","GO:0032963","GO:0051384")
# cairo_ps(filename = "figures_genecode/diff_2step_day7_GOBP.ps",fallback_resolution = 600)
plotEnrichFactor(results_GO[GO_Interest[-c(2,3,5,10)]],title = 'Enriched GO BP gene sets',n=30)
# dev.off()

#####################################################################################################
########## Pair-wise: day 10 treatment effect ##################

# 2 step pair-wise comparison, Genes DE in both comparisons:
# - day 10 treated vs day 10 untreated cells, and 
# - day 10 treated vs day 3 untreated

results = lfcShrink(dds = dds_diff_full,res = 
                      results(dds_diff_full,name = "day10.treatmentDiff",test="Wald"),type = 'ashr',quiet = T) %>% 
  data.frame() %>% rownames_to_column("SYMBNAME") %>% 
  filter(abs(log2FoldChange)>1 & padj<0.05 & SYMBNAME %in% (lfcShrink(dds = dds_diff_full,
                                                                      res = results(dds_diff_full,
                                                                                    contrast = list(c("day10","day10.treatmentDiff")),
                                                                                    test="Wald"),type = 'ashr',quiet = T) %>% 
                                                              data.frame() %>% rownames_to_column("SYMBNAME") %>% 
                                                              filter(abs(log2FoldChange)>1 & padj<0.05) %>% 
                                                              pull(SYMBNAME)))
ht1=plotggplot_Heatmap(dds_diff_full,DE = results$SYMBNAME,genesOfInterest = genesOfInterest,condCol = c('day','treatment'), name = '(treatment effect at day 10)',comp = 'fullModel',colTopAnnot = colTopAnnot,samplesOrder = samplesOrder,split = 2)
cairo_ps(filename = "figures_genecode/diff_2step_day10_heatmap.ps",fallback_resolution = 600)
draw(ht1)
dev.off()
include_graphics("figures_genecode/diff_2step_day10_heatmap.ps")

geneList=bitr((sapply(results$SYMBNAME, function(x){unlist(strsplit(x,split = '|',fixed = T))[2]})),fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)[,'ENTREZID']
results_GO=enrichGO(gene = geneList,OrgDb = org.Hs.eg.db,ont = 'BP',readable = T)
# results_GO=simplify(results_GO)
plotEnrichFactor(results_GO,title = 'Enriched GO BP gene sets',n=30)
GO_Interest=c("GO:0006029","GO:0050746","GO:0006026","GO:0030641","GO:0042475","GO:0046849","GO:0014065","GO:0006027",
              "GO:0046850","GO:0051453","GO:0090382","GO:0051452","GO:0090383","GO:0014031","GO:0048771",
              "GO:0042476","GO:0070371","GO:0034103","GO:0032963","GO:0051384","GO:0019471","GO:0048813","GO:0046427","GO:0022617")
cairo_ps(filename = "figures_genecode/diff_2step_day10_GOBP.ps",fallback_resolution = 600)
plotEnrichFactor(results_GO[GO_Interest[which(GO_Interest %in% results_GO$ID)]],title = 'Enriched GO BP gene sets',n=30)
dev.off()

cairo_ps(filename = "figures_genecode/diff_genesOfInterest_exp.ps",fallback_resolution = 600,width = 10)
degPlotLocal(dds = dds_diff_full,genes =genesOfInterest_ext,xs = 'day',n=40,group = 'treatment',xsLab = "time",ysLab = "normalized counts \n(Median of ratios method)",splitGeneName = 2,reverseColor = T)
dev.off()
include_graphics("figures_genecode/diff_genesOfInterest_exp.ps")