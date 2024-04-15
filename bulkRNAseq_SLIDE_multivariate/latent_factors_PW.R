# Libraries
library('TxDb.Mmulatta.UCSC.rheMac3.refGene')
library('biomaRt')
library('stringr')
library('orthogene')
library('ReactomePA')
library('clusterProfiler')
library('AnnotationDbi')
library('enrichplot')
library('org.Hs.eg.db')
library('RColorBrewer')
library('WebGestaltR')
library('ggplot2')

# Set paths
proj_path <- "/ix/cigcore/proj/mauricio"
dataPath <- file.path(proj_path, "result", "SLIDE_tpm/SLIDE02212024")
experiments <- c("Day28_2Vs1/0.01_0.5_out_spec0.2", "Day28_2Vs3/delta0.02/spec0.25/0.02_2_out")

mapDET <- list()
ORA <- list()

for(e1 in experiments){
  # result path
  cp_path <- file.path(dataPath, e1, "pathway_analysis_CP-ORA")
  dir.create(cp_path, recursive = TRUE)
  ## Load data
  ## list of features from each Z
  features <- readRDS(file.path(dataPath, e1, "plotSigGenes_data.RDS"))
  if(e1 == "Day28_2Vs3"){
    features <- readRDS(file.path(dataPath, e1, "delta0.02/plotSigGenes_data.RDS"))
  }
  ## allFactors <- Reduce(function(x, y) merge(x, y, all=TRUE), features)
  allFactors <- features[!duplicated(features$names), ]
  rownames(allFactors) <- allFactors$names
  
  ortho_df <- orthogene::convert_orthologs(allFactors, gene_input="names", 
                                           gene_output = "columns", input_species="mmulatta",
                                           output_species = "hsapiens", non121_strategy = "keep_both_species")
  gene_list <- ortho_df$ortholog_gene
  ## projectName <- paste0(g1,"_",a1)
  mapDET[[e1]] <- ortho_df
  
  # Get gene names and convert to entrez id
  entrez_genes <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  
  # Reactome pathway over-representation analysis
  reactome <- enrichPathway(gene = entrez_genes$ENTREZID, organism = "human", qvalueCutoff=0.05, readable=TRUE)
  reactome_df <- as.data.frame(reactome)
  ORA[[paste0(e1, "_Reactome")]] <- reactome
  write.csv(as.data.frame(reactome), file.path(cp_path, "ORA_Reactome.csv"))
  
  # Plot enrichment result - dotplot
  if(nrow(reactome_df) != 0){
    dotplot(reactome, showCategory = 20)
    ggsave(file.path(cp_path, "ORA_Reactome_dotplot_selected.pdf"), device = "pdf",
           width = 15, height = 15, units = "cm")
  } else {
    print("No enriched pathways found using Reactome database")
  }
  
  # KEGG pathway over-representation analysis
  kegg <- enrichKEGG(gene = entrez_genes$ENTREZID, organism = "hsa",  qvalueCutoff=0.05,
                     use_internal_data = FALSE)
  kegg2 <- setReadable(kegg, OrgDb = "org.Hs.eg.db", keyType="ENTREZID")
  kegg_df <- as.data.frame(kegg2)
  ORA[[paste0(e1, "_KEGG")]] <- kegg
  kegg2@result$Description <- gsub(pattern = " - Homo sapiens (human)", replacement = "", 
                                   kegg2@result$Description, fixed = T)    
  write.csv(as.data.frame(kegg2), file.path(cp_path, "ORA_KEGG.csv"))
  
  # Plot enrichment result - dotplot
  if(nrow(kegg_df) != 0){
    dotplot(kegg2, showCategory = 20)
    ggsave(file.path(cp_path,  "ORA_KEGG_dotplot.pdf"), device = "pdf")
  } else {
    print("No enriched pathways found using KEGG database")
  }
}

