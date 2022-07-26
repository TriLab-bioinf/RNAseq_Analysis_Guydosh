# Plotting functions
draw_enrichment_barplot <- function(my_enricher_obj, 
                                    my_pathway_counts = 10, 
                                    file_name_prefix = "enricher_barplot", 
                                    my_width = 11, my_height = 8){
  my_enricher_obj@result$p.adjust <- as.numeric(format(my_enricher_obj@result$p.adjust,  digits=3))
  sp <- ggbarplot(my_enricher_obj@result[my_pathway_counts, ], 
                  x = "Description", 
                  y = "Count",
                  fill = "p.adjust",          # change fill color by cyl
                  color = "white",            # Set bar border colors to white
                  sort.val = "desc",          # Sort the value in dscending order
                  sort.by.groups = FALSE,     # Don't sort inside each group
                  x.text.angle = 90,          # Rotate vertically x axis texts
                  rotate = TRUE,
                  ggtheme = theme_minimal(),
                  ylab = c("Gene counts"),
  ) + gradient_fill("RdYlBu")
  
  print(sp)
  
  ggsave2(filename = paste0(file_name_prefix,".pdf"), plot = sp, width = my_width, height = my_height)
  
  return(sp)
}

draw_GSEA_barplot <- function(my_gsea_obj, 
                              my_pathway_counts = 10, 
                              file_name_prefix = "gsea_barplot", 
                              my_width = 11, my_height = 8){
  my_gsea_obj@result$p.adjust <- as.numeric(format(my_gsea_obj@result$p.adjust,  digits=3))
  sp <- ggbarplot(head(my_gsea_obj@result, n = my_pathway_counts), 
                  x = "Description", 
                  y = "NES",
                  fill = "p.adjust",          # change fill color by cyl
                  color = "white",            # Set bar border colors to white
                  sort.val = "desc",          # Sort the value in dscending order
                  sort.by.groups = FALSE,     # Don't sort inside each group
                  x.text.angle = 90,          # Rotate vertically x axis texts
                  rotate = TRUE,
                  ggtheme = theme_minimal(),
                  ylab = c("Normalized Enrichment Score (NES)"),
                  lab.size = 3
  ) + gradient_fill("RdYlBu")
  print(sp)
  ggsave2(filename = paste0(file_name_prefix,".pdf"), plot = sp, width = my_width, height = my_height)
  return(sp)
}

# Function for overrepresentation analysis (ONLY FOR HUMAN GENES SO FAR)
run_overrepresentation_analysis <- function(gene_set, dds_res, 
                                            analysis_name = "dds.result", 
                                            gs_name = "MSigDB", 
                                            adj_p_cutoff = 0.05, 
                                            type = "general",
                                            log2FC_cutoff = 0
                                            ){
  # Inputs: gene_set=Gene set from MSigDB, 
  # genes=vector with ensembl gene IDs.
  # analysis_name=prefix describing DE comparison for output file names
  # gs_name=identifier of the MSigDB used (gene_set)
  
  # Check if dds_res is empty and return if that is the case
  if(dim(subset(dds_res, padj < 0.05))[1] == 0){ return('')}
  
  # Convert NAs in padj column to "1"s
  dds_res$padj[is.na(dds_res$padj)] <- 1
  
  # Create output directories
  my_path = paste0("./ClusterProfiler/",analysis_name,"/")
  dir.create(path = my_path, recursive = TRUE)
  
  # Save Log2FC vals in log2fc_symb.df using gene symbols as IDs. 
  log2fc_symb.df <- as.data.frame(dds_res[,c("log2FoldChange","padj")])
  log2fc_symb.df$symbols <- replace_gene_acc_by_symbol_ids(rownames(log2fc_symb.df),
                                                           return_all = TRUE)
  
  geneList.filtered <- rownames(dds_res)[abs(dds_res$log2FoldChange) > log2FC_cutoff 
                                         & dds_res$padj < adj_p_cutoff]
  
  
  if(type == "reactome"){
      # Replace ensembl IDs by Entrez IDs for compatibility with gsePathway
      entrez_ids.df <- clusterProfiler::bitr(geneID = geneList.filtered, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")
      msig_h.oa <- enrichPathway(gene = entrez_ids.df$ENTREZID, qvalueCutoff = adj_p_cutoff, readable = TRUE)
    } else {
    msig_h.oa <- enricher(gene = geneList.filtered, TERM2GENE=gene_set,  pvalueCutoff = adj_p_cutoff)
  }
    
  if(dim(msig_h.oa@result)[1] > 0){
      # Convert ENSEMBL gene IDs to gene symbols
      msig_h.oa.gs <- setReadable(msig_h.oa, org.Hs.eg.db, keyType = "ENSEMBL")
      
      
      if(dim(msig_h.oa.gs)[1] > 0){
        
        # Write output to a file (OR stands for Overrepresenation analysis). Change file names for diff conditions.
        write.table(as.data.frame(msig_h.oa.gs), file = paste0(my_path,"OR_msigdb.",gs_name,".",analysis_name,".txt"), sep = "\t")
        
        # Dotplot
        p1 <- dotplot(msig_h.oa.gs, showCategory=36) + 
          ggtitle(paste0(analysis_name," dotplot for ORA")) + 
          theme(axis.text=element_text(size=3))
        
        ggsave2(filename = paste0(my_path,"OA_msigdb.",gs_name,".",analysis_name,"_dotplot.pdf"), plot = p1, width = 11, height = 8)
        
        # Barplot
        p3 <- draw_enrichment_barplot(my_enricher_obj = msig_h.oa.gs, 
                                      my_pathway_counts = 36, 
                                      file_name_prefix = paste0(my_path,
                                                                "OA_msigdb.",
                                                                gs_name,".",
                                                                analysis_name,
                                                                "_barplot"
                                                                ), 
                                      my_width = 11, 
                                      my_height = 11)
        print(p1)
        print(p3)
        
        ##############################################
        if(dim(msig_h.oa)[1] > 1){ # We need at least 2 nodes for a network or tree plot
          distance_matrix <- pairwise_termsim(msig_h.oa, showCategory = 400)
          
          # Get median Log2FC per identified pathway
          my_median_log2fc <- c()
          for (g in distance_matrix@result$geneID){
            g.vec <- strsplit(g, "/")[[1]]
            log2fc.median <- median(subset(log2fc_symb.df, rownames(log2fc_symb.df) %in% g.vec | log2fc_symb.df$symbols %in% g.vec)[,"log2FoldChange"])
            my_median_log2fc <- c(my_median_log2fc,log2fc.median)
          }
          # Add median Log2FC column
          if(length(my_median_log2fc) == 0){my_median_log2fc <- 0}
          distance_matrix@result$median.log2fc <- my_median_log2fc
          
          # Network plot
          p6 <- emapplot(distance_matrix, 
                         repel = T, 
                         showCategory = 200, 
                         legend_n = 5, 
                         min_edge = 0.4 , 
                         color = "median.log2fc", 
                         cex_label_category = 0.4,
                         node_label = "category", label_format = 20)
          ggsave2(filename = paste0(my_path,
                                    "OA_msigdb.",
                                    gs_name,".",
                                    analysis_name,"_network.pdf"), 
                  plot = p6, width = 11, height = 8)
          
          # Treeplots
          number_of_categories = min(80, as.vector(table(distance_matrix@result$p.adjust < 0.05)[['TRUE']]))
          p12 <- treeplot(distance_matrix, showCategory = 80, 
                          nCluster = round(2 * sqrt(number_of_categories), digits = 0), 
                          color = "median.log2fc", nWords = 0)
          ggsave2(filename = paste0(my_path,
                                    "OA_msigdb.",
                                    gs_name,".",
                                    analysis_name,"_tree.pdf"), 
                  plot = p12, width = 11, height = max(8, 8 * (number_of_categories/40)))
          
          print(p6)
          print(p12)
        }
        #####################################
      }
    }  
  
}

# Function for enrichment analysis (ONLY FOR HUMAN GENES SO FAR)
run_enrichment_analysis <- function(gene_set, 
                                    geneList, 
                                    analysis_name = "dds.result", 
                                    gs_name = "MSigDB",
                                    adj_p_cutoff = 0.05, 
                                    type = "general", no_plots = FALSE){
  # Inputs: gene_set=Gene set from MSigDB, 
  # geneList=a list with ensembl_ids as names and decremented order of Log2FCs as values
  # analysis_name=prefix describing DE comparison for output file names
  # gs_name=identifier of the MSigDB used (gene_set)
  
  # Create output directories
  my_path = paste0("./ClusterProfiler/",analysis_name,"/")
  dir.create(path = my_path, recursive = TRUE)
  
  # Run Gene Set Enrichment Analysis. Reactome type display reactome IDs and Descriptions properly
  if(type == "reactome"){
    # Replace ensembl IDs by Entrez IDs for compatibility with gsePathway
    entrez_ids.df <- clusterProfiler::bitr(geneID = names(geneList), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")
    geneList.entrez <- geneList[c(entrez_ids.df$ENSEMBL)]
    names(geneList.entrez) <- entrez_ids.df$ENTREZID
    msig_h.gsea <- gsePathway(geneList = geneList.entrez, pvalueCutoff = adj_p_cutoff, verbose = TRUE, eps = 0)
    if(dim(msig_h.gsea@result)[1] > 0){ 
      # Convert Entrez gene IDs to gene symbols
      msig_h.gsea.gs <- setReadable(msig_h.gsea, org.Hs.eg.db, keyType = "ENTREZID")
      if(no_plots){return(msig_h.gsea.gs)}
    } 
  } else {
    msig_h.gsea <- GSEA(geneList, TERM2GENE = gene_set, pvalueCutoff = adj_p_cutoff, verbose = TRUE, eps = 0)
    if(dim(msig_h.gsea@result)[1] > 0){ 
      # Convert Entrez gene IDs to gene symbols
      msig_h.gsea.gs <- setReadable(msig_h.gsea, org.Hs.eg.db, keyType = "ENSEMBL")
      if(no_plots){return(msig_h.gsea.gs)}
    }  
    
  } 
  
  if(dim(msig_h.gsea@result)[1] > 0){ 
    # Write output to a file (OR stands for Overrepresenation analysis). Change file names for diff conditions.
    write.table(as.data.frame(msig_h.gsea.gs), file = paste0(my_path,"GSEA_msigdb.",gs_name,".",analysis_name,".txt"), sep = "\t")
    
    # Generate plots summarizing enrichment results (You might need to tweak the height and width parameters of the plot to make it fit within the page limits and look nice, otherwise it might shrink to fit within the page)
    
    # Upset plots
    # filter out signatures with NES < 0 since they are likely due to RNAseL degradation to simplify the plot
    msig_h.gsea_bt_0 <- filter(msig_h.gsea, NES > 0) # filter out signatures with NES < 0
    msig_h.gsea_lt_0 <- filter(msig_h.gsea, NES < 0) # filter out signatures with NES < 0
    
    if(dim(msig_h.gsea_bt_0@result)[1] > 0){
      p_ups_up <- upsetplot(msig_h.gsea_bt_0, n=36) + scale_x_upset(n_intersections = 30) 
      ggsave2(filename = paste0(my_path,"GSEA_msigdb.",gs_name,".",analysis_name,"_UP_upset_plot.pdf"), 
             plot = p_ups_up, width = 22, 
             height = 4 * (max(2,min(36, length(msig_h.gsea_bt_0@result$ID))/15)))
    }
    if(dim(msig_h.gsea_lt_0@result)[1] > 0){
      p_ups_down <- upsetplot(msig_h.gsea_lt_0, n=36) + scale_x_upset(n_intersections = 30) 
      ggsave2(filename = paste0(my_path,"GSEA_msigdb.",gs_name,".",analysis_name,"_DOWN_upset_plot.pdf"), 
             plot = p_ups_down, width = 22, 
             height = 4 * (max(2,min(36, length(msig_h.gsea_lt_0@result$ID))/15)))
    }
    
    # Dotplot
    p2 <- dotplot(msig_h.gsea.gs, showCategory=36, 
                  font.size = 7, 
                  color = "p.adjust", 
                  label_format = 100) + 
      ggtitle(paste0(my_prefix," dotplot for GSEA"))
    
    ggsave2(filename = paste0(my_path,"GSEA_msigdb.",gs_name,".",analysis_name,"_dotplot.pdf"), 
            plot = p2, width = 15, height = 8)
    
    # Barplot
    p4 <- draw_GSEA_barplot(my_gsea_obj =  msig_h.gsea.gs, my_pathway_counts = 36, 
                            file_name_prefix = paste0(my_path,"GSEA_msigdb.",gs_name,".",analysis_name,"_barplot"), 
                                                my_width = 11, 
                                                my_height =11
                            )
    
    # Network plots and tree plots for all, up and downregulated pathways
    if(dim(msig_h.gsea)[1] > 1){ # We need at least 2 nodes for a network or tree plot
      distance_matrix <- pairwise_termsim(msig_h.gsea, showCategory = 400)
      
      # Network plot
      p6 <- emapplot(distance_matrix, 
                     repel = T, 
                     showCategory = 200, 
                     legend_n = 5, 
                     min_edge = 0.2 , 
                     color = "NES", 
                     cex_label_category = 0.4,
                     node_label = "category", label_format = 20)
      ggsave2(filename = paste0(my_path,"GSEA_msigdb.",gs_name,".",analysis_name,"_network.pdf"), 
              plot = p6, width = 11, height = 8)
      
      # Treeplots
      number_of_categories = min(80, as.vector(table(distance_matrix@result$p.adjust < 0.05)[['TRUE']]))
      p12 <- treeplot(distance_matrix, showCategory = 80, 
                      nCluster = 2 * sqrt(number_of_categories), 
                      color = "NES", nWords = 0)
      ggsave2(filename = paste0(my_path,"GSEA_msigdb.",gs_name,".",analysis_name,"_tree.pdf"), 
              plot = p12, width = 11, height = 8 * (number_of_categories/40))
    }
    
    if(dim(msig_h.gsea_bt_0)[1] > 1){ # We need at least 2 nodes for a network or tree plot
      distance_matrix.up <- pairwise_termsim(msig_h.gsea_bt_0, showCategory = 400)
      
      # Network plot
      p8 <- emapplot(distance_matrix.up, 
                     repel = T, 
                     showCategory = 200, 
                     legend_n = 5, 
                     min_edge = 0.2 , 
                     color = "NES", 
                     cex_label_category = 0.4,
                     node_label = "category", label_format = 20)
      ggsave2(filename = paste0(my_path,"GSEA_msigdb.",gs_name,".",analysis_name,"_network_UP.pdf"), 
              plot = p8, width = 11, height = 8)
      
      # Treeplots
      number_of_categories.up = min(80, length(distance_matrix.up@result$ID))
      p14 <- treeplot(distance_matrix.up, showCategory = 80, 
                      nCluster = 2 * sqrt(number_of_categories.up), 
                      color = "NES", nWords = 0)
      ggsave2(filename = paste0(my_path,"GSEA_msigdb.",gs_name,".",analysis_name,"_tree_UP.pdf"), 
              plot = p14, width = 11, height = 8 * (number_of_categories.up/40))
      
    }
    
    if(dim(msig_h.gsea_lt_0)[1] > 1){ # We need at least 2 nodes for a network or tree plot
      distance_matrix.down <- pairwise_termsim(msig_h.gsea_lt_0, showCategory = 400)
      
      # Network plot
      p10 <- emapplot(distance_matrix.down, 
                     repel = T, 
                     showCategory = 200, 
                     legend_n = 5, 
                     min_edge = 0.2 , 
                     color = "NES", 
                     cex_label_category = 0.4,
                     node_label = "category", label_format = 20)
      ggsave2(filename = paste0(my_path,"GSEA_msigdb.",gs_name,".",analysis_name,"_network_DOWN.pdf"), 
              plot = p10, width = 11, height = 8)
      
      # Treeplots
      number_of_categories.down = min(80, length(distance_matrix.down@result$ID))
      p16 <- treeplot(distance_matrix.down, showCategory = 80, 
                      nCluster = 2 * sqrt(number_of_categories.down), 
                      color = "NES", nWords = 0)
      ggsave2(filename = paste0(my_path,"GSEA_msigdb.",gs_name,".",analysis_name,"_tree_DOWN.pdf"), 
              plot = p16, width = 11, height = 8 * (number_of_categories.down/40))
      
    }
    
    
    
    print(p2)
    print(p4)
    print(p_ups_up)
    print(p_ups_down)
    print(p6)
    print(p8)
    print(p10)
    print(p12)
    print(p14)
    print(p16)
  }
  return(msig_h.gsea.gs)
}

go_overrep <- function(dds_res, my_file, golevel = 3, 
                       onthology = "BP", qvalue = 0.05, 
                       fdr = 0.05, log2FC_cutoff = 0){
  
  # Check if dds_res is empty and return if that is the case
  if(dim(dds_res)[2] == 0){ return('')}
  
  # Output directory
  my_path = paste0("./GO/",my_file,"/")
  
  # Save Log2FC vals in log2fc_symb.df using gene symbols as IDs. 
  log2fc_symb.df <- as.data.frame(dds_res[,c("log2FoldChange","padj")])
  log2fc_symb.df$symbols <- replace_gene_acc_by_symbol_ids(rownames(log2fc_symb.df),
                                                             return_all = TRUE)
  
  # Convert Ensembl IDs to EntrezIDs
  rownames(dds_res) <- replace_gene_acc_by_gb_ids(rownames(dds_res),
                                                  return_all = TRUE)
  # Discard genes missing Entrez Gene IDs
  keep <- grep('ENSG', rownames(dds_res), invert = TRUE)
  dds_res <- dds_res[keep,]
  
  geneList.filtered <- rownames(dds_res)[abs(dds_res$log2FoldChange) > log2FC_cutoff 
                                         & dds_res$padj < qvalue]
  
  # Check if dds_res is empty and return if that is the case
  if(length(geneList.filtered) == 0){ return('')}
  
  ego <- enrichGO(gene = geneList.filtered,
                  OrgDb= "org.Hs.eg.db",
                  ont = onthology, # Either MF, BP, CC or ALL
                  pAdjustMethod = "BH",
                  qvalueCutoff  = fdr,
                  keyType = 'ENTREZID',
                  readable      = TRUE,
                  minGSSize = 10,
                  maxGSSize = 500
        )
  
  write.table(ego, file = paste0(my_path,"GO_OVERREP.",onthology,".",my_file,".txt"), sep = "\t", col.names=NA)
  
  # Network plots and tree plots for all, up and downregulated pathways
  if(dim(ego)[1] > 1){ # We need at least 2 nodes for a network or tree plot
    distance_matrix <- pairwise_termsim(ego, showCategory = 400)
    
    # Get median Log2FC per identified pathway
    my_median_log2fc <- c()
    for (g in distance_matrix@result$geneID){
      g.vec <- strsplit(g, "/")[[1]]
      log2fc.median <- median(subset(log2fc_symb.df, log2fc_symb.df$symbols %in% g.vec)[,"log2FoldChange"])
      my_median_log2fc <- c(my_median_log2fc,log2fc.median)
    }
    # Add median Log2FC column to distance matrix
    if(length(my_median_log2fc) == 0){my_median_log2fc <- 0}
    distance_matrix@result$median.log2fc <- my_median_log2fc
    
    # Add median Log2FC column to ego object
    if(length(my_median_log2fc) == 0){my_median_log2fc <- 0}
    ego@result$median.log2fc <- my_median_log2fc
    
    # Network plot
    p6 <- emapplot(distance_matrix, 
                   repel = T, 
                   showCategory = 200, 
                   legend_n = 5, 
                   min_edge = 0.4 , 
                   color = "median.log2fc", 
                   cex_label_category = 0.4,
                   node_label = "category", label_format = 20)
    ggsave2(filename = paste0(my_path,"GO_OVERREP.",onthology,".",my_file,"_network.pdf"), 
            plot = p6, width = 11, height = 8)
    
    # Treeplots
    number_of_categories = min(80, as.vector(table(distance_matrix@result$p.adjust < 0.05)[['TRUE']]))
    p12 <- treeplot(distance_matrix, showCategory = 80, 
                    nCluster = 2 * sqrt(number_of_categories), 
                    color = "median.log2fc", nWords = 0)
    
    ggsave2(filename = paste0(my_path,"GO_OVERREP.",onthology,".",my_file,"_tree.pdf"), 
            plot = p12, width = 11, height = 8 * (number_of_categories/40))
    
    print(p6)
    print(p12)
  }
  
  return(ego)
}

go_classification <- function(dds_res, my_file, golevel = 3, onthology = "BP", log2FC_cutoff = 0, alpha = 1){
  
  # Output directory
  my_path = paste0("./GO/",my_file,"/")
  
  # Convert Ensembl IDs to EntrezIDs
  rownames(dds_res) <- replace_gene_acc_by_gb_ids(rownames(dds_res),
                                                  return_all = TRUE)
  # Discard genes missing Entrez Gene IDs
  keep <- grep('ENSG', rownames(dds_res), invert = TRUE)
  dds_res <- dds_res[keep,]
  
  geneList.filtered <- rownames(dds_res)[abs(dds_res$log2FoldChange) >= log2FC_cutoff & dds_res$padj <= alpha]
  
  print(paste( "Total number of Entrez gene IDs = ",length(geneList.filtered)))
  
  ggo <- groupGO(gene     = geneList.filtered,
                 OrgDb    = "org.Hs.eg.db",
                 ont      = onthology,  # Either MF, BP or CC
                 level    = golevel, # The higher the number, the more specific the GO terms.
                 readable = TRUE, # if readable is TRUE, the gene IDs will mapping to gene symbols.
                 keyType = 'ENTREZID'
  )
  write.table(ggo, file = paste0(my_path,"GO_CLASSIF.",onthology,".",my_file,".txt"), sep = "\t",col.names=NA)
  
  return(ggo)
}

go_gsea <- function(dds_res, my_file, golevel = 3, 
                    onthology = "BP", qvalue = 0.05, 
                    fdr = 0.05){
  
  # Check if dds_res is empty and return if that is the case
  if(dim(dds_res)[2] == 0){ return('')}
  
  # Skip analysis if there are no genes with padj < 0.05
  geneList.filtered <- rownames(dds_res)[abs(dds_res$log2FoldChange) > log2FC_cutoff 
                                         & dds_res$padj < qvalue]
  
  # Check if dds_res is empty and return if that is the case
  if(length(geneList.filtered) == 0){ return('')}
  
  # Output directory
  my_path = paste0("./GO/",my_file,"/")
  
  # Convert Ensembl IDs to EntrezIDs
  rownames(dds_res) <- replace_gene_acc_by_gb_ids(rownames(dds_res),
                                                  return_all = TRUE)
  # Discard genes missing Entrez Gene IDs
  keep <- grep('ENSG', rownames(dds_res), invert = TRUE)
  dds_res <- dds_res[keep,]
  
  # Sort genes by Log2FC (decreasing order)
  geneList = list()
  dds_res.sorted <- dds_res[order(dds_res$log2FoldChange, decreasing = TRUE ), ]
  geneList <- dds_res.sorted$log2FoldChange
  names(geneList) <- rownames(dds_res.sorted)
  
  go.gsea <- gseGO(gene = geneList,
                   OrgDb= "org.Hs.eg.db",
                   ont = onthology, # Either MF, BP, CC or ALL
                   pAdjustMethod = "BH",
                   pvalueCutoff  = fdr,
                   keyType = 'ENTREZID',
                   minGSSize = 10,
                   maxGSSize = 500,
                   eps = 0
                  )
  go.gsea.gs <- setReadable(go.gsea , org.Hs.eg.db, keyType = "ENTREZID")
  write.table(go.gsea.gs, file = paste0(my_path,"GO_GSEA.",onthology,".",my_file,".txt"), sep = "\t", col.names=NA)
  
  # Generate plots summarizing enrichment results (You might need to tweak the height and width parameters of the plot to make it fit within the page limits and look nice, otherwise it might shrink to fit within the page)
  
  if(length(go.gsea.gs@result$NES) > 0){
    
    # Upset plots
    # filter out signatures with NES < 0 
    go.gsea_bt_0 <- filter(go.gsea.gs, NES > 0) # filter out signatures with NES < 0
    go.gsea_lt_0 <- filter(go.gsea.gs, NES < 0) # filter out signatures with NES < 0
    
    if(dim(go.gsea_bt_0@result)[1] > 0){
      p_ups_up <- upsetplot(go.gsea_bt_0, n=10) + scale_x_upset(n_intersections = 30) 
      ggsave2(filename = paste0(my_path,"GSEA_GO.",onthology,".",my_file,"_UP_upset_plot.pdf"), 
              plot = p_ups_up, width = 22, 
              height = 4 * (max(2,min(36, length(go.gsea_bt_0@result$ID))/15)))
    }
    if(dim(go.gsea_lt_0@result)[1] > 0){
      p_ups_down <- upsetplot(go.gsea_lt_0, n=36) + scale_x_upset(n_intersections = 30) 
      ggsave2(filename = paste0(my_path,"GSEA_GO.",onthology,".",my_file,"_DOWN_upset_plot.pdf"), 
              plot = p_ups_down, width = 22, 
              height = 4 * (max(2,min(36, length(go.gsea_lt_0@result$ID))/15)))
    }
    # Dotplot
    p2 <- dotplot(go.gsea.gs, showCategory=36, 
                  font.size = 7, 
                  color = "p.adjust", 
                  label_format = 100) + 
      ggtitle(paste0(my_prefix," dotplot for GSEA"))
    
    ggsave2(filename = paste0(my_path,"GSEA_GO.",onthology,".",my_file,"_dotplot.pdf"), 
            plot = p2, width = 15, height = 8)
    
    # Barplot
    p4 <- draw_GSEA_barplot(my_gsea_obj =  go.gsea.gs, my_pathway_counts = 36, 
                            file_name_prefix = paste0(my_path,"GSEA_GO.",onthology,".",my_file,"_barplot"), 
                            my_width = 11, 
                            my_height =11
    )
    
    # Network plots and tree plots for all, up and downregulated pathways
    if(dim(go.gsea.gs)[1] > 1){ # We need at least 2 nodes for a network or tree plot
      distance_matrix <- pairwise_termsim(go.gsea.gs, showCategory = 400)
      
      # Network plot
      p6 <- emapplot(distance_matrix, 
                     repel = T, 
                     showCategory = 200, 
                     legend_n = 5, 
                     min_edge = 0.2 , 
                     color = "NES", 
                     cex_label_category = 0.4,
                     node_label = "category", label_format = 20)
      ggsave2(filename = paste0(my_path,"GSEA_GO.",onthology,".",my_file,"_network.pdf"), 
              plot = p6, width = 11, height = 8)
      
      # Treeplots
      number_of_categories = min(80, as.vector(table(distance_matrix@result$p.adjust < 0.05)[['TRUE']]))
      p12 <- treeplot(distance_matrix, showCategory = 80, 
                      nCluster = 2 * sqrt(number_of_categories), 
                      color = "NES", nWords = 0)
      ggsave2(filename = paste0(my_path,"GSEA_GO.",onthology,".",my_file,"_tree.pdf"), 
              plot = p12, width = 11, height = 8 * (number_of_categories/40))
    }
    
    if(dim(go.gsea_bt_0)[1] > 1){ # We need at least 2 nodes for a network or tree plot
      distance_matrix.up <- pairwise_termsim(go.gsea_bt_0, showCategory = 400)
      
      # Network plot
      p8 <- emapplot(distance_matrix.up, 
                     repel = T, 
                     showCategory = 200, 
                     legend_n = 5, 
                     min_edge = 0.2 , 
                     color = "NES", 
                     cex_label_category = 0.4,
                     node_label = "category", label_format = 20)
      ggsave2(filename = paste0(my_path,"GSEA_GO.",onthology,".",my_file,"_network_UP.pdf"), 
              plot = p8, width = 11, height = 8)
      
      # Treeplots
      number_of_categories.up = min(80, length(distance_matrix.up@result$ID))
      p14 <- treeplot(distance_matrix.up, showCategory = 80, 
                      nCluster = 2 * sqrt(number_of_categories.up), 
                      color = "NES", nWords = 0)
      ggsave2(filename = paste0(my_path,"GSEA_GO.",onthology,".",my_file,"_tree_UP.pdf"), 
              plot = p14, width = 11, height = 8 * (number_of_categories.up/40))
      
    }
    
    if(dim(go.gsea_lt_0)[1] > 1){ # We need at least 2 nodes for a network or tree plot
      distance_matrix.down <- pairwise_termsim(go.gsea_lt_0, showCategory = 400)
      
      # Network plot
      p10 <- emapplot(distance_matrix.down, 
                      repel = T, 
                      showCategory = 200, 
                      legend_n = 5, 
                      min_edge = 0.2 , 
                      color = "NES", 
                      cex_label_category = 0.4,
                      node_label = "category", label_format = 20)
      ggsave2(filename = paste0(my_path,"GSEA_GO.",onthology,".",my_file,"_network_DOWN.pdf"), 
              plot = p10, width = 11, height = 8)
      
      # Treeplots
      number_of_categories.down = min(80, length(distance_matrix.down@result$ID))
      p16 <- treeplot(distance_matrix.down, showCategory = 80, 
                      nCluster = 2 * sqrt(number_of_categories.down), 
                      color = "NES", nWords = 0)
      ggsave2(filename = paste0(my_path,"GSEA_GO.",onthology,".",my_file,"_tree_DOWN.pdf"), 
              plot = p16, width = 11, height = 8 * (number_of_categories.down/40))
      
    }
    
    print(p2)
    print(p4)
    print(p_ups_up)
    print(p_ups_down)
    print(p6)
    print(p8)
    print(p10)
    print(p12)
    print(p14)
    print(p16)
    
    #return(go.gsea.gs)
  }  
  return(go.gsea.gs)
}

PAPER_plot_gsea <- function(gsea_result, comparison_id, 
                          analysis_type = "GSEA",
                          my_x_label = "Gene sets",
                          my_path = "./PAPER"){
  
  dir.create(path = my_path, showWarnings = FALSE)
  my_file <- comparison_id
  
  if(str_count(string = my_file, pattern = "_PAPER") == 0){
    my_file <- paste0(my_file,"_PAPER")
  }
  
  if(str_count(string = my_file, pattern = analysis_type) == 0){
    my_file <- paste0(my_file,"_", analysis_type)
  }
  
  REPORT_TOP_CUTOFF = 10 # per NES score direction (UP or DOWN)
  
  if(length(gsea_result@result$NES) > 0){
    
    # Upset plots
    # filter out signatures with NES < 0 
    gsea_result_bt_0 <- filter(gsea_result, NES > 0) # filter out signatures with NES > 0
    gsea_result_lt_0 <- filter(gsea_result, NES < 0) # filter out signatures with NES < 0
    
    
    # Select representative gene set per clusters based on NES * gene_number
    p_top_go_up <- get_gsea_clusters(gsea_result_bt_0, my_report_top_cutoff = REPORT_TOP_CUTOFF) # UP gene sets
    p_top_go_down <- get_gsea_clusters(gsea_result_lt_0, my_report_top_cutoff = REPORT_TOP_CUTOFF) # DOWN gene sets
    
    # Merge UP and DOWN gene sets
    if(!is.null(dim(p_top_go_up)[1]) & !is.null(dim(p_top_go_down)[1])){
      my_top_go_per_cluster <- bind_rows(list(p_top_go_up, p_top_go_down))
    } else if (!is.null(dim(p_top_go_up)[1])){
      my_top_go_per_cluster <- p_top_go_up
    } else if (!is.null(dim(p_top_go_down)[1])){
        my_top_go_per_cluster <- p_top_go_down
    } else {return(paste("The following gave no results:", my_file))}
    
    # Make barplot for PAPER
    p.treebased_barplot <- ggbarplot(my_top_go_per_cluster, 
                                     x = "label", 
                                     y = "color",
                                     fill = my_top_go_per_cluster$my_custom_color, #"count",          # change fill color by cyl
                                     color = "white",            # Set bar border colors to white
                                     sort.by.groups = FALSE,     # Don't sort inside each group
                                     x.text.angle = 90,          # Rotate vertically x axis texts
                                     rotate = TRUE,
                                     ggtheme = theme_classic(),
                                     ylab = c("Normalized Enrichment Score (NES)"),
                                     lab.size = 1
    ) + #+ gradient_fill("RdYlBu")
      xlab(my_x_label) + 
      scale_y_continuous(breaks=seq(
                                    round(min(my_top_go_per_cluster$color), digits = 0) - 1,
                                    round(max(my_top_go_per_cluster$color), digits = 0) + 1,
                                    1
                                    )
                        ) # + scale_x_discrete(labels = function(x) str_wrap(x, width = 65, ))
    
    # Save plot
    ggsave2(filename = paste0("barplot_",my_file,".pdf"), plot = p.treebased_barplot, path = my_path, width = 12 )
    
  }  
  
}

# Network plots and tree plots for all, up and downregulated pathways
get_gsea_clusters <- function(my_gsea_result, my_report_top_cutoff = 10){
  
  if(dim(my_gsea_result)[1] > 1){ # We need at least 2 nodes for a network or tree plot
    dm <- pairwise_termsim(my_gsea_result, showCategory = 400)
    
    # Treeplots
    number_of_categories = min(80, as.vector(table(dm@result$p.adjust < 0.05)[['TRUE']]))
    tree.p <- treeplot(dm, showCategory = number_of_categories,
                       color = 'NES',
                       nCluster = min(my_report_top_cutoff, number_of_categories), # 2 * sqrt(number_of_categories), 
                       nWords = 0)
    
    # ggsave2(filename = paste0(my_path,"GSEA_GO.",onthology,".",my_file,"_tree.pdf"), 
    #        plot = tree.p, width = 11, height = 8 * (number_of_categories/40))
  

  
    # keep just top GO per cluster based on abs(NES). Cluster ID = group, NES = color
    my_data <- tree.p$data
    my_data <- mutate(my_data, score = color * count)
    
    # Split hits into up and downs so there is a fare representation of both in the final plot for the paper
    my_top_go_per_cluster <- subset(my_data, isTip == TRUE) %>% group_by(group) %>% slice_max(order_by = abs(score), n = 1, with_ties = FALSE)
    
    # Add bar's colors 
    my_top_go_per_cluster$my_custom_color <- ifelse(my_top_go_per_cluster$color > 0, "darkgreen", "red")
    
    # Sort barchar by NES
    my_top_go_per_cluster <- arrange(my_top_go_per_cluster, desc(color) )
    
    # Edit hallmark gene set labels
    if(str_count(string = my_top_go_per_cluster$label[1], pattern = "HALLMARK_") > 0){
      my_top_go_per_cluster$label <- str_replace_all(string = my_top_go_per_cluster$label, pattern = "_", replacement = " ") %>% str_replace(pattern = "HALLMARK ", "") %>% tolower()
    } else if (str_count(string = my_top_go_per_cluster$label[1], pattern = "_TARGET_GENES") > 0 ){
      my_top_go_per_cluster$label <- str_replace_all(string = my_top_go_per_cluster$label, pattern = "_TARGET_GENES", replacement = " target genes")
    } else {
      my_top_go_per_cluster$label <- str_replace_all(string = my_top_go_per_cluster$label, pattern = "_", replacement = " ")
    }
    
    return(my_top_go_per_cluster)
  } else {print("get_gsea_clusters: GSEA RESULT SHOULD HAVE MORE THAN ONE SIGNIFICANT GENE SET!")}
}

add_lo2fc <- function(dm){ # dm = distance matrix
  # Get median Log2FC per identified pathway
  my_median_log2fc <- c()
  for (g in dm@result$geneID){
    g.vec <- strsplit(g, "/")[[1]]
    log2fc.median <- median(subset(log2fc_symb.df, log2fc_symb.df$symbols %in% g.vec)[,"log2FoldChange"])
    my_median_log2fc <- c(my_median_log2fc,log2fc.median)
  }
  # Add median Log2FC column
  if(length(my_median_log2fc) == 0){my_median_log2fc <- 0}
  dm@result$median.log2fc <- my_median_log2fc
  return(dm)
}

plot_simplified_network <- function(enrichment_result, top_categories){
  
  REPORT_TOP_CUTOFF <- top_categories
  
  #distance_matrix <- pairwise_termsim(enrichment_result, showCategory = 400)
  #distance_matrix <- add_lo2fc(dm = distance_matrix)
  
  # add NES column with median.log2fc values to be compatible with get_gsea_clusters function
  if(is.null(enrichment_result@result$NES)){
    enrichment_result@result$NES <- enrichment_result@result$median.log2fc
  }
  
  ########################## FIX BEGIN
  # process up and down pathways separately,
  # so both changes are represented in the simplified plot
  
  # filter out signatures with NES < 0 
  enrichment_result_bt_0 <- filter(enrichment_result, NES > 0) # filter out signatures with NES > 0
  enrichment_result_lt_0 <- filter(enrichment_result, NES < 0) # filter out signatures with NES < 0
  
  
  # Select representative gene set per clusters based on NES * gene_number
  p_top_go_up <- get_gsea_clusters(my_gsea_result = enrichment_result_bt_0, my_report_top_cutoff = REPORT_TOP_CUTOFF) # UP gene sets
  p_top_go_down <- get_gsea_clusters(my_gsea_result = enrichment_result_lt_0, my_report_top_cutoff = REPORT_TOP_CUTOFF) # DOWN gene sets
  
  # Merge UP and DOWN gene sets
  if(!is.null(dim(p_top_go_up)[1]) & !is.null(dim(p_top_go_down)[1])){
    my_top_go_per_cluster <- bind_rows(list(p_top_go_up, p_top_go_down))
  } else if (!is.null(dim(p_top_go_up)[1])){
    my_top_go_per_cluster <- p_top_go_up
  } else if (!is.null(dim(p_top_go_down)[1])){
    my_top_go_per_cluster <- p_top_go_down
  } else {return(paste("The following gave no results:", my_file))}
  
  
  ########################## FIX END
  
  
  #my_top_go_per_cluster <- get_gsea_clusters(my_gsea_result = enrichment_result, 
  #                                           my_report_top_cutoff = top_categories)
  top_cluster_ids <- my_top_go_per_cluster$label
  
  ego.filtered = enrichment_result %>% filter(Description %in% top_cluster_ids)
  distance_matrix.filter <- pairwise_termsim(ego.filtered, showCategory = 400)
  #distance_matrix.filter <- add_lo2fc(distance_matrix.filter)
  
  # network
  my_max_range <- max(abs(distance_matrix.filter$median.log2fc))
  p.net.filter <- emapplot(distance_matrix.filter, 
                           repel = T, 
                           showCategory = 200, 
                           legend_n = 5, 
                           min_edge = 0.4 , 
                           color = "median.log2fc", 
                           cex_label_category = 0.4,
                           node_label = "category",
                           shadowtext = FALSE,
                           label_style = "ggforce",
                           label_format = 10) + 
    scale_edge_width(range = c(0.1,2)) +
    
    #scale_colour_gradientn(colours = myPalette(100), limits=c(-1 * my_max_range , my_max_range))
    
    return(p.net.filter)
}
