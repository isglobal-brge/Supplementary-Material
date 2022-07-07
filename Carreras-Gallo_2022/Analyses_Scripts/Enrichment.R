######################################
## ENRICHMENT SUBSTANCE CONSUMPTION ##
######################################

#Load libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(ggupset)
library(missMethyl)

`%notin%` <- Negate(`%in%`)

#Function for doing the enrichment in missMethyl package
call_gometh <- function(cpgs, topcpgs, collection){
  res <- gometh(
    sig.cpg = rownames(cpgs),
    all.cpg = rownames(topcpgs),
    collection = collection,
    array.type = "EPIC",
    plot.bias = FALSE,
    prior.prob = TRUE,
    anno = NULL,
    equiv.cpg = TRUE,
    fract.counts = TRUE,
    genomic.features = "ALL",
    sig.genes = FALSE)
    
    #Sort by P-Value
    res <- res[order(res$P.DE),]
    return(res[which(res$P.DE<0.01),])
}

#Load DisGeNET database
gda <- read.delim("/PROJECTES/GENOMICS/TruDiagnostic/Analysis/Cannabis/curated_gene_disease_associations.tsv.gz")
disease2gene <- gda[, c("diseaseId", "geneId")]
disease2name <- gda[, c("diseaseId", "diseaseName")]

#Create a function
enrich <- function(substance, model, pval){
  #Set wd
  mainDir <- paste0("/PROJECTES/GENOMICS/TruDiagnostic/Analysis/Uniform_analysis/",substance,"/",model)
  setwd(mainDir)
  
  #Load topcpgs
  load(list.files(pattern="topcpgs.+Rdata"))

  #Create enrichment directory if it does not exist
  subDir <- "Enrichment"
  dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
  setwd(file.path(mainDir, subDir))
  
  #Select top CpG sites
  cpgs <- topcpgs[which(topcpgs$P.Value<pval),]
  nrow(cpgs)
  
  ### missMethyl enrichment (GO and KEGG)
  
  res_go <- call_gometh(cpgs, topcpgs, "GO")
  res_kegg <- call_gometh(cpgs, topcpgs, "KEGG")
  writexl::write_xlsx(res_go,path=paste0("GO_gometh_",pval,".xlsx"))
  writexl::write_xlsx(res_kegg,path=paste0("KEGG_gometh_",pval,".xlsx"))
  
  #Transform gene symbols to Entrez and Ensembl
  not <- ""
  genes <- unique(strsplit(paste(cpgs$HGNC_GeneSymbol, collapse=";"),";")[[1]])
  genes <- genes[genes %notin% not]

  entrez <- bitr(genes, fromType = "SYMBOL",
                 toType = "ENTREZID",
                 OrgDb = "org.Hs.eg.db")
  ensembl <- bitr(genes, fromType = "SYMBOL",
                  toType = "ENSEMBL",
                  OrgDb = "org.Hs.eg.db")

  ### DigGeNET enrichment

  ans.dis <- enricher(entrez$ENTREZID, TERM2GENE=gda[, c("diseaseId", "geneId")],
                      TERM2NAME=gda[, c("diseaseId", "diseaseName")])
  as.data.frame(ans.dis)

  list <- list()
  if (nrow(ans.dis)>0){
    p1 <- dotplot(ans.dis, showCategory=20) + ggtitle(paste0("DisGeNET enrichment for ",substance," effect"))
    list <- append(list,list(p1))
  }


  ### KEGG enrichment

  ans.kegg <- enrichKEGG(gene = entrez$ENTREZID,
                         organism = 'hsa',
                         pvalueCutoff = 0.05)
  as.data.frame(ans.kegg)

  if (nrow(ans.kegg)>0){
    p2 <- dotplot(ans.kegg, showCategory=20) + ggtitle(paste0("KEGG enrichment for ",substance," effect"))
    list <- append(list,list(p2))
  }


  ### GO enrichment

  ans.go <- enrichGO(gene = ensembl$ENSEMBL, ont = "BP",
                     OrgDb ="org.Hs.eg.db",
                     keyType = "ENSEMBL",
                     readable = TRUE,
                     pvalueCutoff = 0.05)
  as.data.frame(ans.go)

  if (nrow(ans.go)>0){
    p3 <- dotplot(ans.go, showCategory=20) + ggtitle(paste0("GO enrichment for ",substance," effect"))
    list <- append(list,list(p3))
  }

  #Figure with all the plots
  if (length(list)==1){
    png(paste0("enrichment_",pval,".png"), width=1000, height = 700, units="px")
    print(list[[1]])
    dev.off()
  }

  if (length(list)==2){
    png(paste0("enrichment_",pval,".png"), width=700, height = 900, units="px")
    print(cowplot::plot_grid(plotlist=list, ncol=1,labels="AUTO"))
    dev.off()
  }

  if (length(list)==3){
    png(paste0("enrichment_",pval,".png"), width=700, height = 1300, units="px")
    print(cowplot::plot_grid(plotlist=list, ncol=1,labels="AUTO"))
    dev.off()
  }
  return(paste0("Total CpGs: ",nrow(cpgs),"\nTotal genes: ",length(genes),"\nEntrez genes: ",nrow(entrez),"\nEnsembl genes: ",nrow(ensembl),"\n"))
}

for (model in c("5levels","Never_vs_Ever","Never_vs_min3week")){
  enrich(substance="Marijuana",model=model,pval=1e-4)
  enrich(substance="Alcohol",model=model,pval=1e-4)
}

for (model in c("7levels","None_vs_Any","None_vs_min1day")){
  enrich(substance="Tobacco",model=model,pval=1e-4)
}
