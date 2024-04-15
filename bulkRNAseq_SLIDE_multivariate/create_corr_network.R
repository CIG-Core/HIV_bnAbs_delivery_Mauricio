## this code makes correlation networks ##
# load packages
library(tidyverse)
library(qgraph) ## for making the network

create_corr_network = function(path, select_LF, plot_label = NULL, geneList) {
  
  yaml_input = yaml::yaml.load_file(path)
  
  dir_name = yaml_input$out_path
  
  x_path = yaml_input$x_path
  y_path = yaml_input$y_path
  
  x_mat = as.matrix(read.csv(x_path, row.names = 1, check.names = F))
  y_mat = as.matrix(read.csv(y_path, row.names = 1, check.names = F))
  
  df = cbind.data.frame(y_mat, x_mat)
  
  
  ############### uncomment if you want to use the ER path instead
  check_for_file = function(path, file_pattern) {

    f = list.files(path, recursive = T, full.names = T, pattern = file_pattern)
    if( length(f) > 0 ) {
      return(f[1])
    } else {
      cat("\n Folder must have exactly one file with desired name: ", pattern, "\n")
      return(NULL)
    }
  }
  
  # get_x_and_y = function(path) {
  #   x = check_for_file(x_path, file_pattern = ".csv")
  #   y = check_for_file(y_path, file_pattern = ".csv")
  #   
  #   load_matrix_csv = function(path) {
  #     return(read.csv(path, row.names = 1))
  #   }
  #   
  #   if (all(!is.null(x), !is.null(y))) {
  #     return(cbind.data.frame(load_matrix_csv(y), load_matrix_csv(x)))
  #   } else {
  #     df = check_for_file(path, file_pattern = df_filename_RDS_pattern)
  #     if (!is.null(df)) {
  #       return(readRDS(df))
  #     } else {
  #       return(NULL)
  #     }
  #   }
  # }
  
  # df = get_x_and_y(path)
  # if (is.null(df)) {
  #   cat("\n Failed to load dataframes. Check path \n")
  #   return()
  # } else {
  #   names(df)[1] = "y"
  # }
  
  sig_genes_data = check_for_file(dir_name, file_pattern = "plotSigGenes_data.RDS")
  sig_genes_data = readRDS(sig_genes_data)
  
  all_LFs = c()
  
  if ( !is.null(sig_genes_data) ) {
    all_LFs = stringr::str_to_upper(sig_genes_data$names)
    
  } else {
    cat("\n Couldn't find genes in significant latent factors. Check path. \n")
    return()
  }
  
  color_code = stringr::str_to_lower(sig_genes_data$color)
  
  x = df[, -1]
  
  if ( dim(x)[2] <= 1) {
    cat("\n Latent factor genes not found in X data. Check path \n")
    return()
  }
  
  dir_name = paste0(dir_name, "LF_correlation_plots/")
  
  if ( !dir.exists(dir_name) ) {
    dir.create(dir_name)
  }
  
  original_wd = getwd()
  
  sg_temp = sig_genes_data %>% filter(sig_genes_data$lf_num %in% select_LF)
  
  sg_temp <- sg_temp[abs(sg_temp$AUCs) > 0.4, ]
  sg_temp$shape <- "NA"
  sg_temp[sg_temp$names %in% geneList, "shape"] <- "triangle"
  sg_temp[sg_temp$shape == "NA", "shape"] <- "circle"
  
  x_temp = x[, colnames(x) %in% sg_temp$names]
  
  x_cor = cor(as.matrix(x_temp))
  
  col_auc = round(apply(x_temp, 2, function(xs) glmnet:::auc(y_mat, as.matrix(xs))), 2)
  
  node_color = ifelse(col_auc > 0.55, "cyan", ifelse(col_auc < 0.45, "lightgreen", "lightgray"))
  ## node_color = ifelse(sg_temp$lf_num == select_LF[1], "#40006D", "#59A14F")
  
  node_shape = ifelse(is.na(sg_temp$color), "circle", ifelse(sg_temp$color == "Red", "triangle", "square"))
  ## node_shape = ifelse(sg_temp$lf_num == select_LF[1], "circle",  "triangle")
  names(node_shape) <- sg_temp$names
  node_shape2 <- node_shape[colnames(x_cor)]
  
  node_groups <- node_shape2
  node_groups[node_groups == "square"] <- "Box Z"
  node_groups[node_groups == "triangle"] <- "Box Y"
  
  # change wd cause weird
  setwd(dir_name)
  plot_title = ifelse(is.null(plot_label), "", plot_label)
  
  # plot_title = paste0("LF #", LF, "\n", plot_title, "\nTriangle = Up\nSquare = Down")
  plot_title = paste0("LF # ", str_c(select_LF, collapse = " and "))
  LF = paste0("LF", str_c(select_LF, collapse = "_"))
  
  # this won't get plotted. the only important part is to have the 'groups' argument be whatever you want to group
  # by - e.g if you want the shapes clustered together or if you want nodes with a specific color clustered together
  pl_shape = qgraph(x_cor, filename = LF, groups = node_color,
                    legend.mode = "groups",
                    color = node_color,
                    shape = sg_temp$shape,
                    layout = "spring",
                    # minimum=0.1, # use minimum instead of threshold
                    labels = colnames(x_cor),
                    label.scale.equal=TRUE,label.prop=0.8, fade=FALSE,
                    ## shape="ellipse",
                    DoNotPlot = TRUE,
                    posCol="darkred", negCol="magenta2",filetype='pdf',height=10,width=10)
  
  # change the 'groups' argument here to whatever you want the legend to show
  pl = qgraph(x_cor, filename = paste0("grouped_LF_shape_cluster_", LF),
              color = node_color,
              shape = sg_temp$shape,
              edge.width = 1,
              node.width = 0.75,
              node.height = 0.75, fade = FALSE, # turn off fade if you want the edge color to be constant
              minimum=0.5,
              labels = colnames(x_cor), label.cex = 1.5,
              # title = plot_title,
              layout = pl_shape$layout, # use the plot layout above
              label.scale.equal=TRUE,label.prop=0.95,shape="ellipse", 
              posCol="darkred", negCol="magenta2",filetype='pdf',height=5,width=8)
  
  setwd(original_wd)
}