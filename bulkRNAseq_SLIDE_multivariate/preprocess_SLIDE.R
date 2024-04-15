# Library
library(matrixStats)

# Set paths
proj_path <- "/ix/cigcore/proj/mauricio"
data_path <- file.path(proj_path, "data")
result_path <- file.path(proj_path, "result")
xVarPath <- file.path(result_path, "12-12-2023-DE_Box1vs2vs3_All")

# dataDirs <- c(file.path(data_path, "202212_OBS03"), file.path(data_path, "2022_Kiran", "Mauricio Data"), 
#               file.path(data_path, "2021_fastq"))
allCounts <- read.csv(file.path(xVarPath, "tpm_names.csv"), row.names = 1)
allCounts <- allCounts[allCounts$gene_name != "", ]
t2g <- as.data.frame(allCounts$gene_name)
rownames(t2g) <- rownames(allCounts)
colnames(t2g) <- "gene_name"

timePoints <- c("PostEarly", "Late")
sInfos <- list()
counts <- list()
for(t1 in timePoints){
    metadataPath <- file.path(data_path, paste0("samples_", t1, ".csv"))
    sInfo <- read.csv(metadataPath, row.names = 1)
    if(t1 == "PostEarly"){
      sInfos[["Day1"]] <- sInfo[sInfo$Days_Sample_Collected == 1,]
      counts[["Day1"]] <- allCounts[, sInfos[["Day1"]]$Samples]
      sInfos[["Day3"]] <- sInfo[sInfo$Days_Sample_Collected != 1,]
      counts[["Day3"]] <- allCounts[, sInfos[["Day3"]]$Samples]
    } else {
      sInfos[["Day28"]] <- sInfo
      counts[["Day28"]] <- allCounts[, sInfo$Samples]
    }
}

s2Vs1 <- c(1,2)
s2Vs3 <- c(2,3)
groups <- c("s2Vs1", "s2Vs3")
## sInfosGroup <- list()
## countsGroup <- list()
## varMatrices <- list()

for(s1 in names(sInfos)){
  for(g1 in groups){
    sInfo2 <- sInfos[[s1]]
    sGrouped <- sInfo2[sInfo2$Group %in% get(g1), c("Samples", "Group")] 
    sGrouped$Group[sGrouped$Group == 2] <- 0
    sGrouped$Group[sGrouped$Group == 3] <- 1
    ## sInfosGroup[[paste0(g1, s1)]]  <- sGrouped
    write.csv(sGrouped, file.path(data_path, paste0(g1, s1, ".csv")), row.names = FALSE)
    
    ## Processing TPM files to get high variable genes
    cGroup <- as.matrix(counts[[s1]][,sGrouped$Samples])
    
    ## Calculate variances in genes
    varsData <- rowVars(cGroup)
    ## Add variance column to TPM matrix
    countsMat <- cbind(cGroup, "Vars" = varsData)
    ## Ordering based on variances
    orderCounts <- countsMat[order(countsMat[,"Vars"], decreasing = TRUE), ]
    orderCounts <- orderCounts[, -ncol(orderCounts)]
    
    ## Imputing or removing rows with too many zeroes
    ## Replacing 0.00 with NAs
    orderCounts[orderCounts == 0.000] <- NA 
    naNum <- colSums(is.na(orderCounts))
    
    ## Identifying columns with too many NAs
    tooManyNAs <- naNum > nrow(orderCounts)/2
    orderCounts <- orderCounts[ ,!tooManyNAs]
    orderCounts[is.na(orderCounts)] <- 0.00 
    
    ## Adding gene names and filtering based on annotation and unique gene ids
    orderCounts <- as.data.frame(cbind(orderCounts, "gene_name"= t2g[rownames(orderCounts), "gene_name"]))
    orderCounts <- orderCounts[!duplicated(orderCounts$gene_name), ]
    
    ## Adding gene_names as row names
    rownames(orderCounts) <- orderCounts$gene_name
    orderCounts <- orderCounts[,-ncol(orderCounts)]
    
    ## Subsetting top 1000 genes
    topVar <- t(orderCounts[1:2000, ])
    ## varMatrices[[paste0(s1, g1)]] <- topVar
    write.csv(topVar, file.path(data_path, paste0("topTpm", s1, g1, "_2000.csv")))
  }
}

