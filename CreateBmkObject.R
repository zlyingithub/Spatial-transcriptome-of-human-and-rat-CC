library(Seurat)
library(dplyr)
CreateS1000Object <- function(
  matrix_path,
  png_path,
  spot_radius = NULL,
  min.cells = 5,
  min.features = 100
  ){
  expr <- Seurat::Read10X(matrix_path, cell.column = 1)
  object <- Seurat::CreateSeuratObject(counts = expr,
                               assay = 'Spatial',
                               min.cells=min.cells,
                               min.features=min.features)
  #Image zoom rate
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
  #read png
  png <- png::readPNG(png_path)
  zoom_scale <-  cal_zoom_rate(dim(png)[2], dim(png)[1])
  #read barcode pos file
  ReadBarcodePos <- function(barcode_pos_path){
    barcode_pos <- read.table(gzfile(barcode_pos_path),header = F) %>%
      dplyr::rename(Barcode = V1 , pos_w = V2, pos_h = V3)
    return(barcode_pos)
  }
  #get barcode pos file path
  barcode_pos_path <- paste0(matrix_path,'/barcodes_pos.tsv.gz')
  barcode_pos <- ReadBarcodePos(barcode_pos_path = barcode_pos_path)
  barcode_pos <- barcode_pos %>% dplyr::filter(., Barcode %in% rownames(object@meta.data))
  #make spatial coord file for seurat S4 class
  coord <- data.frame(tissue = 1,
                      row = barcode_pos$pos_h,
                      col = barcode_pos$pos_w,
                      imagerow = barcode_pos$pos_h,
                      imagecol = barcode_pos$pos_w)
  rownames(coord) <- barcode_pos$Barcode
  #spot radius
  spot_radius_lib <- c(0.00063, 0.00179, 0.0027, 0.0039, 0.004, 0.0045, 0.005, NA, NA, NA, NA, NA, 0.0120)
  if(is.null(spot_radius)){
    spot_radius <- spot_radius_lib[as.numeric(gsub('L', '', strsplit(tail(strsplit(matrix_path, '/')[[1]],1), '_')[[1]][1]))]
  }else{
    spot_radius = spot_radius
  }
  #object
  sample1 <-  new(Class = "VisiumV1",
                  image = png,
                  scale.factors = Seurat::scalefactors(zoom_scale, 100, zoom_scale, zoom_scale),
                  coordinates = coord,
                  spot.radius = spot_radius,
                  assay = 'Spatial',
                  key = "sample1_")
  object@images <- list(sample1 = sample1)

  return(object)
}

