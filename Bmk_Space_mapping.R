Create_object <- function(FilePath,
                          barcode_pos_file,   
                          out_path,
                          png_path,
                          min.cells=10, 
                          min.features=100,
                          dims = 1:30, 
                          resolution = 0.5,
                          point_size = 3.2,
                          width = 12, 
                          height = 5,
                          alpha = 1,
                          alpha_continuous = c(0.9, 1),
                          output_polt = TRUE,
                          Cluster = FALSE,
                          label = FALSE,
                          UMI_stat = FALSE,
                          nFeature_stat = FALSE,
                          Gene_stat = FALSE,
                          Custom_gene = FALSE,
                          dark_background = FALSE,
                          gene_list = c('a','b'),
                          top_gene = 1,
                          min.pct = 0.25,   
                          logfc.threshold = 0.25, 
                          markpic_width = 7,
                          markpic_height = 5,
                          Single_cluster = FALSE,
                          color_gradientn = c('#1200ff', '#b892ff', '#faf7ff', '#ffa68a', '#ff1b07'),
                          col =  c("#F56867","#FEB915","#C798EE","#59BE86","#7495D3","#D1D1D1","#6D1A9C","#15821E","#3A84E6","#997273","#787878","#DB4C6C","#9E7A7A","#554236","#AF5F3C","#93796C","#F9BD3F","#DAB370","#877F6C","#268785","#ed1299", "#09f9f5", "#246b93", "#cc8e12", "#d561dd", "#c93f00", "#ddd53e","#4aef7b", "#e86502", "#9ed84e", "#39ba30", "#6ad157", "#8249aa", "#99db27", "#e07233", "#ff523f","#ce2523", "#f7aa5d", "#cebb10", "#03827f", "#931635", "#373bbf", "#a1ce4c", "#ef3bb6", "#d66551","#1a918f", "#ff66fc", "#2927c4", "#7149af" ,"#57e559" ,"#8e3af4" ,"#f9a270" ,"#22547f", "#db5e92","#edd05e", "#6f25e8", "#0dbc21", "#280f7a", "#6373ed", "#5b910f" ,"#7b34c1" ,"#0cf29a" ,"#d80fc1")
                          
)

{
  #加载包
  suppressMessages({
    if(!require(Seurat))install.packages("Seurat")
    if(!require(tidyverse))install.packages("tidyverse")
    if(!require(ggplot2))install.packages("ggplot2")
    if(!require(patchwork))install.packages("patchwork")
    if(!require(reshape2))install.packages("reshape2")
    if(!require(ggdark))install.packages("ggdark")
    if(!require(cluster))install.packages("cluster")
    if(!require(png))install.packages("png")
    if(!require(ggpubr))install.packages("ggpubr")
  })
  
  #创建输出文件目录
  if(!dir.exists(out_path)){dir.create(path=out_path)}
  #创建object对象
  if(!exists('object')){
    expr <- Read10X(FilePath,cell.column = 1)
    object <- CreateSeuratObject(counts = expr,
                                 assay = "Spatial",
                                 min.cells=min.cells,
                                 min.features=min.features)
    object <- SCTransform(object, 
                          assay = "Spatial",
                          verbose = FALSE,
                          return.only.var.genes = FALSE)
    object <- object %>%
      RunPCA(verbose = F,assay = "SCT") %>% 
      FindNeighbors(reduction = "pca",dims = dims,verbose = F) %>% 
      FindClusters(verbose = F,resolution = resolution) %>%
      RunUMAP(reduction = "pca",dims = dims) %>%
      RunTSNE(reduction = "pca",dims = dims)
  } 
  
  #图片缩放
  cal_zoom_rate <- function(width, height){
    std_width = 1000
    std_height = std_width / (46 * 31) * (46 * 36 * sqrt(3) / 2.0)
    if (std_width / std_height > width / height){
      scale = width / std_width
    }
    else{
      scale = height / std_height
    }
    return(scale)
  }
  png <- readPNG(png_path)
  zoom_scale = cal_zoom_rate(dim(png)[1], dim(png)[2])
  
  #barcode pos file
  barcode_pos <- read.table(gzfile(barcode_pos_file),header = F) %>% 
    rename(Barcode = V1 , pos_w = V2, pos_h = V3) %>% 
    mutate(across(c(pos_w, pos_h), ~ .x*zoom_scale))
  
  #viol function
  viol <- function(Count_file, title){
    p <- ggplot(Count_file, aes(x = 'Identity', y = Count_file[,2]))+
      geom_violin(fill = '#ed0000')+
      geom_boxplot(width = 0.1)+
      geom_point(size = 0.5, position = position_jitter(width = 0.3), color = '#002240')+
      theme_classic()+
      labs(x = NULL, title = title, y = NULL)+
      theme(plot.title = element_text(hjust = 0.5))
    return(p)
  }
  #space_heatmap function
  Space_heatmap <- function(Count_file, title = 'heatmap'){
    p <- ggplot(Count_file, aes(x = pos_w, y = (dim(png)[2] - pos_h)))+
      background_image(png)+
      geom_point(shape = 16, alpha = alpha, size = point_size, aes(color = Count_file[,2]))+
      coord_cartesian(xlim = c(0, dim(png)[2]), ylim = c(0, dim(png)[1]),expand = FALSE)+
      scale_color_gradientn(colours = color_gradientn)+
      #scale_color_gradientn(colours = c('#555fa9','#8ad0a1','#daeb97','#f9f5b2','#ef7e4b','#a4144c'))+
      theme_void()+
      labs(colour = title)

    return(p)
  }
  

  #cluster picture
  if (Cluster == TRUE) {
    
    #cluster file 
    Cluster_df <- object@meta.data %>%
      dplyr::select(.,seurat_clusters) %>%
      rownames_to_column(.,var = "Barcode") %>%
      arrange(.,seurat_clusters) 
    
    #barcode pos & cluster file 
    barcode_pos_cluster <- left_join(Cluster_df, barcode_pos, by = 'Barcode') %>% 
      rename(cluster =seurat_clusters)
    barcode_pos_cluster$cluster <- as.factor(barcode_pos_cluster$cluster)
    #umap file
    UMAP_df <- object[['umap']]@cell.embeddings %>% 
      as.data.frame() %>% 
      rownames_to_column(var = 'Barcode') %>% 
      left_join(Cluster_df, by = 'Barcode') %>% 
      rename(cluster = seurat_clusters) %>% 
      select(Barcode, cluster, everything())
    UMAP_df$cluster <- as.factor(UMAP_df$cluster)
    
    #p1
    cluster_pic <- ggplot(barcode_pos_cluster, aes(x = pos_w, y = (dim(png)[2] - pos_h)))+
      background_image(png)+
      geom_point(shape = 16, alpha = alpha, size = point_size, aes(color = cluster))+
      coord_cartesian(xlim = c(0, dim(png)[2]), ylim = c(0, dim(png)[1]),expand = FALSE)+
      scale_color_manual(values = col)+
      theme_void()
    
    #p2
    umap_pic <- ggplot(UMAP_df, aes(x = UMAP_1, y = UMAP_2))+
      geom_point(aes(color = cluster))+
      scale_color_manual(values = col)+
      labs(colour = "cluster")+
      theme_bw(base_size = 20)
    
    umap_cluster <- umap_pic + cluster_pic + plot_layout(widths = c(1.8,2))
    
    if(label == TRUE){
      label_func <- function(df){
        label_data <- df %>% 
          group_by(cluster) %>% 
          do(model = pam(.[colnames(df[-1:-2])], 1, stand = T)) %>% 
          ungroup() %>% group_by(cluster) %>% 
          do(map_df(.$model, broom::tidy)) %>% ### 整理模型数据
          ungroup() %>% select(cluster, colnames(df[3]), colnames(df[4])) %>% 
          data.frame() %>% 
          dplyr::rename(x.center = colnames(df[3]), y.center = colnames(df[4]), cluster = cluster)
        label_data$cluster <- c(0:(max(as.numeric(df$cluster))-1))
        label_data$cluster <- as.factor(label_data$cluster)
        return(label_data)
      }
      #p1_label data
      cluster_label_data <- label_func(barcode_pos_cluster)
      #p1_label pic
      cluster_label_pic <- cluster_pic + 
        geom_label(data = cluster_label_data, aes(label = cluster, x = x.center, y = (dim(png)[1]-y.center) ,fill = cluster), 
                   colour = 'white', fontface = "bold")+
        scale_fill_manual(values = col) + guides(fill = 'none')
      
      #p2_labe_data
      UMAP_label_data <- label_func(UMAP_df)
      #p2_label pic 
      umap_label_pic <- umap_pic + 
        geom_label(data = UMAP_label_data, aes(label = cluster, x = x.center, y = y.center,fill = cluster), 
                   colour = 'white', fontface = "bold")+
        scale_fill_manual(values = col) + guides(fill = 'none')
      
      #all
      umap_cluster_label <- umap_label_pic + cluster_label_pic + plot_layout(widths = c(1.8,2))
    }
    
    if(output_polt == TRUE){
      #p all out
      pdf(file = paste(out_path, 'umap_cluster.pdf', sep = '/'), width = width, height = height)
      print(umap_cluster)
      dev.off()
      ggsave(umap_cluster, filename = paste(out_path, 'umap_cluster.png', sep = '/'), width = width, height = height, dpi = 300) 
      
      # p1 out 
      pdf(file = paste(out_path, 'cluster.pdf', sep = '/'), width = width/3.8*2, height = height)
      print(cluster_pic)
      dev.off()
      ggsave(cluster_pic, filename = paste(out_path, 'cluster_pic.png', sep = '/'), width = width/3.8*2, height = height, dpi = 300) 
      
      # p2 out 
      pdf(file = paste(out_path, 'umap.pdf', sep = '/'), width = width/3.8*1.8, height = height)
      print(umap_pic)
      dev.off()
      ggsave(umap_pic, filename = paste(out_path, 'umap.png', sep = '/'), width = width/3.8*1.8, height = height, dpi = 300) 
    }
    
    if(label == TRUE){
      #p label all out
      pdf(file = paste(out_path, 'umap_cluster_label.pdf', sep = '/'), width = width, height = height)
      print(umap_cluster_label)
      dev.off()
      ggsave(umap_cluster_label, filename = paste(out_path, 'umap_cluster_label.png', sep = '/'), width = width, height = height, dpi = 300) 
      
      
      # p1 label out 
      pdf(file = paste(out_path, 'cluster_label.pdf', sep = '/'), width = width/3.8*2, height = height)
      print(cluster_label_pic)
      dev.off()
      ggsave(cluster_label_pic, filename = paste(out_path, 'cluster_label.png', sep = '/'), width = width/3.8*2, height = height, dpi = 300) 
      
      
      
      # p2 label out 
      pdf(file = paste(out_path, 'umap_label.pdf', sep = '/'), width = width/3.8*1.8, height = height)
      print(umap_label_pic)
      dev.off()
      ggsave(umap_label_pic, filename = paste(out_path, 'umap_label.png', sep = '/'), width = width/3.8*1.8, height = height, dpi = 300)
      }

  } 
  
  
  #UMI counts picture
  if (UMI_stat == TRUE){
    UMI_stat_df <-  object$nCount_Spatial %>% 
      as.data.frame() %>% 
      rownames_to_column(var = 'Barcode') %>% 
      rename(UMI = '.') %>% 
      left_join(barcode_pos, by = 'Barcode')
    # UMI viol
    UMI_viol <- viol(UMI_stat_df, 'UMI_Count')
    # UMI heatmap
    UMI_heatmap <- Space_heatmap(UMI_stat_df, 'nUMI')
    # UMI viol & heatmap
    UMI_viol_heatmap <- UMI_viol + UMI_heatmap + plot_layout(widths = c(1.8,2))
    
    if(output_polt == TRUE){
      # UMI viol & heatmap
      pdf(file = paste(out_path, 'UMI_viol_heatmap.pdf', sep = '/'), width = width, height = height)
      print(UMI_viol_heatmap)
      dev.off()
      ggsave(UMI_viol_heatmap, filename = paste(out_path, 'UMI_viol_heatmap.png', sep = '/'), width = width, height = height, dpi = 300) 
      
      #UMI viol
      pdf(file = paste(out_path, 'UMI_viol.pdf', sep = '/'), width = width/3.8*1.8, height = height)
      print(UMI_viol)
      dev.off()
      ggsave(UMI_viol, filename = paste(out_path, 'UMI_viol.png', sep = '/'), width = width/3.8*1.8, height = height, dpi = 300)
      
      #UMI_heatmap
      pdf(file = paste(out_path, 'UMI_heatmap.pdf', sep = '/'), width = width/3.8*2, height = height)
      print(UMI_heatmap)
      dev.off()
      ggsave(UMI_heatmap, filename = paste(out_path, 'UMI_heatmap.png',sep = '/'), width = width/3.8*2, height = height, dpi = 300)
    }
    
  }
  
  
  #Gene Counts picture
  if (nFeature_stat == TRUE){
    #nFeature file
    nFeature_stat_df <- object$nFeature_Spatial %>% 
      as.data.frame() %>% 
      rownames_to_column(var = 'Barcode') %>% 
      rename(nFeature_Spatial = '.') %>% 
      left_join(barcode_pos, by = 'Barcode')
    # UMI viol
    nFeature_viol <- viol(nFeature_stat_df, 'nFeature_Count')
    # UMI heatmap
    nFeature_heatmap <- Space_heatmap(nFeature_stat_df, 'nFeature')
    # UMI viol & heatmap
    nFeature_viol_heatmap <- nFeature_viol + nFeature_heatmap + plot_layout(widths = c(1.8,2))
    if(output_polt == TRUE){
      # UMI viol & heatmap
      pdf(file = paste(out_path, 'nFeature_viol_heatmap.pdf', sep = '/'), width = width, height = height)
      print(nFeature_viol_heatmap)
      dev.off()
      ggsave(nFeature_viol_heatmap, filename = paste(out_path, 'nFeature_viol_heatmap.png',sep = '/'), width = width, height = height, dpi = 300) 
      
      #UMI viol
      pdf(file = paste(out_path, 'nFeature_viol.pdf', sep = '/'), width = width/3.8*1.8, height = height)
      print(nFeature_viol)
      dev.off()
      ggsave(nFeature_viol, filename = paste(out_path, 'nFeature_viol.png', sep = '/'), width = width/3.8*1.8, height = height, dpi = 300)
      
      #UMI_heatmap
      pdf(file = paste(out_path, 'nFeature_heatmap.pdf', sep ='/'), width = width/3.8*2, height = height)
      print(nFeature_heatmap)
      dev.off()
      ggsave(nFeature_heatmap, filename = paste(out_path, 'nFeature_heatmap.png', sep ='/'), width = width/3.8*2, height = height, dpi = 300)
    }
  }
  
  #single gene picture
  if (Gene_stat == TRUE){
    
    #gene heatmap pic1
    gene_heatmap_func1 <- function(Gene_df, title_name) {
        ggplot(Gene_df, aes(x = pos_w, y = (dim(png)[2] - pos_h)))+
        #background_image(png)+
        geom_point(shape = 16, size = point_size, aes(color = nCount))+
        coord_cartesian(xlim = c(0, dim(png)[2]), ylim = c(0, dim(png)[1]),expand = FALSE)+
        scale_color_gradientn(colours = c('#141414', '#5bb8da', '#bad255', '#e45316', '#841c12'))+
        #scale_color_continuous(low = '#b4cf66', high = '#9b2a14')+
        theme_void()+
        dark_theme_gray()+
        labs(x = NULL, y = NULL, title = title_name)+
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
              axis.text.x = element_blank(), axis.ticks.x = element_blank(),
              plot.title = element_text(size = 14,face = "bold", hjust = 0.5))
      
    }
    #gene heatmap pic2
    gene_heatmap_func2 <- function(Gene_df, title_name) {
      ggplot(Gene_df, aes(x = pos_w, y = (dim(png)[2] - pos_h)))+
        background_image(png)+
        geom_point(shape = 16, size = point_size, aes(color = nCount, alpha = nCount))+
        coord_cartesian(xlim = c(0, dim(png)[2]), ylim = c(0, dim(png)[1]),expand = FALSE)+
        scale_alpha_continuous(range = alpha_continuous)+
        scale_color_gradientn(colours = c('#5c559c', '#77cca3', '#f3fca7', '#fcaa5f', '#be204b'))+
        #scale_color_continuous(low = '#b4cf66', high = '#9b2a14')+
        guides(alpha = 'none')+
        theme_void()+
        labs(x = NULL, y = NULL, title = title_name, alpha = NULL)+
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
              axis.text.x = element_blank(), axis.ticks.x = element_blank(),
              plot.title = element_text(size = 14,face = "bold", hjust = 0.5))
      
    }
      
    
    # single gene file
    if(Custom_gene == TRUE){
      #gene df
      for (i in 1:length(gene_list)) {
        # df
        Gene_df <- object[['Spatial']]@counts %>% as.data.frame() %>% rownames_to_column(var = 'Gene') %>% 
          filter(Gene == gene_list[i]) %>% melt() %>% rename(Barcode = variable, nCount = value) %>% 
          left_join(barcode_pos, by = 'Barcode') %>% dplyr::select(.,-Gene)
        # picture1
        if(dark_background == TRUE){
          gene_heatmap <- gene_heatmap_func1(Gene_df, gene_list[i])
          #save
          pdf(file = paste(out_path, paste0(gene_list[i], '.pdf'), sep = '/'), width = width/3.8*2, height = height)
          print(gene_heatmap)
          dev.off()
          ggsave(gene_heatmap, filename = paste(out_path, paste0(gene_list[i], '.png'), sep = '/'), width = width/3.8*2, height = height, dpi = 300)
        } else {
          gene_heatmap <- gene_heatmap_func2(Gene_df, gene_list[i])
          #save
          pdf(file = paste(out_path, paste0(gene_list[i], '.pdf'), sep = '/'), width = width/3.8*2, height = height)
          print(gene_heatmap)
          dev.off()
          ggsave(gene_heatmap, filename = paste(out_path, paste0(gene_list[i], '.png'), sep = '/'), width = width/3.8*2, height = height, dpi = 300)
        }
      }
  
    } else {
      
      #Markgene
      Markers <- FindAllMarkers(object, only.pos = TRUE, min.pct = min.pct, logfc.threshold = logfc.threshold)
      
      #top gene
      top_gene_df <- Markers %>% group_by(cluster) %>% top_n(n = top_gene, wt = avg_log2FC)
      
      #all top gene pic
      for (i in 1:length(top_gene_df$gene)){
        Gene_df <- object[['Spatial']]@counts %>% as.data.frame() %>% rownames_to_column(var = 'Gene') %>% 
          filter(Gene == top_gene_df$gene[i]) %>% melt() %>% rename(Barcode = variable, nCount = value) %>% 
          left_join(barcode_pos, by = 'Barcode') %>% dplyr::select(.,-Gene)
        
        # heatmap
        gene_heatmap <- gene_heatmap_func2(Gene_df = Gene_df, title_name = top_gene_df$gene[i])
        # save single gene heatmap
        pdf(file = paste(out_path, paste0(top_gene_df$gene[i], '.pdf'), sep = '/'), width = width/3.8*2, height = height)
        print(gene_heatmap)
        dev.off()
        ggsave(gene_heatmap, filename = paste(out_path, paste0(top_gene_df$gene[i], '.png'), sep = '/'), width = width/3.8*2, height = height, dpi = 300)
      }
      
      # Top genes for each subpopulation
      markgene_heatmap <- DoHeatmap(object, features = top_gene_df$gene,size = 2) + NoLegend()
      
      # mark gene vlio
      mark_gene_vilo <- VlnPlot(object, features = top_gene_df$gene, slot = "counts", log = TRUE)
      
      # mark gene tsne
      mark_gene_tsne <- FeaturePlot(object, features = top_gene_df$gene, cols = c("grey", "red"),reduction = "tsne")
      
      # save
      # mark gene file
      write.csv(top_gene_df, file = paste(out_path, paste0('Mark_gene_top', top_gene, '.csv'), sep = '/'), quote = F, row.names = F)
      
      # Top genes for each subpopulation
      pdf(file = paste0(out_path, 'markgene_heatmap.pdf'), width = 8, height = 14*(top_gene/10))
      print(markgene_heatmap)
      dev.off()
      ggsave(markgene_heatmap, filename = paste(out_path, 'markgene_heatmap.png', sep = '/'), width = 8, height = 14, dpi = 300)
      
      
      # mark gene vlio
      pdf(file = paste(out_path, 'markgene_violin.pdf', sep = '/'), width = markpic_width, height = markpic_height)
      print(mark_gene_vilo)
      dev.off()
      ggsave(mark_gene_vilo, filename = paste(out_path, 'markgene_violin.png', sep = '/'), width = markpic_width, height = markpic_height, dpi = 300)
      
      # mark gene tsne
      pdf(file = paste(out_path, 'markgene_tsne.pdf', sep = '/'), width = markpic_width, height = markpic_height)
      print(mark_gene_tsne)
      dev.off()
      ggsave(mark_gene_tsne, filename = paste(out_path,'markgene_tsne.png', sep = '/'), width = markpic_width, height = markpic_width, dpi = 300)
      
    }
    }
  
  #single cluster picture
  if (Single_cluster == TRUE){
    #barcode pos cluster file 
    #cluster file 
    Cluster_df <- object@meta.data %>%
      dplyr::select(.,seurat_clusters) %>%
      rownames_to_column(.,var = "Barcode") %>%
      arrange(.,seurat_clusters) 
    
    #barcode pos & cluster
    barcode_pos_cluster <- left_join(Cluster_df, barcode_pos, by = 'Barcode') %>% 
      rename(cluster =seurat_clusters)
    
    #single cluster picture
    for (i in 0:(length(unique(barcode_pos_cluster$cluster))-1)) {
      p <-  ggplot(barcode_pos_cluster, aes(x = pos_w, y = (dim(png)[2] - pos_h)))+
        background_image(png)+
        geom_point(shape = 16, size = point_size, aes(color = cluster))+
        coord_cartesian(xlim = c(0, dim(png)[2]), ylim = c(0, dim(png)[1]),expand = FALSE)+
        scale_color_manual(values = if_else(unique(barcode_pos_cluster$cluster) == i, 'red', 'grey'))+
        theme_void()+
        #dark_theme_gray()+
        labs(x = NULL, y = NULL, title = paste0('cluster',i))+
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
              axis.text.x = element_blank(), axis.ticks.x = element_blank(),
              legend.position = 'none') + 
        theme(plot.title = element_text(size = 20,face = "bold", hjust = 0.5))
    
      #save pic
      pdf(file = paste(out_path, paste0('cluster', i, '.pdf'), sep = '/'), width = width/3.8*2, height = height)
      print(p)
      dev.off()
      ggsave(p, filename = paste(out_path, paste0('cluster', i, '.png'), sep = '/'), width = width/3.8*2, height = height, dpi = 300)
      
    }
    
  }
  
 return(object) 
}


