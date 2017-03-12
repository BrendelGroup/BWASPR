#' plot_genome_annotation()
#' This function plots a representative figure of genome annotation
#'
#' @param genome_annotation
#'
#'
#' @param mrobj
#'
#' @return A data frame that store the genome annotation data
#'
#' @importFrom
#'
#' @examples
#'   mydatf <- system.file("extdata","Am.dat",package="BWASPR")
#'   myparf <- system.file("extdata","Am.par",package="BWASPR")
#'   myfiles <- setup_BWASPR(datafile=mydatf,parfile=myparf)
#'   AmHE <- mcalls2mkobj(myfiles$datafiles)
#'   genome <- read_genome_info(myfiles$parameters)
#'   genome_annotation_info <- genome_annotation(AmHE, genome, save_output_data = FALSE)
#'   plot_genome_annotation(genome_annotation)
#'
#'
#' @export

plot_genome_annotation <- function(mrobj, genome_annotation_info, if_plot = FALSE, output_file_dir = './BWASPR/R/genome_annotation.pdf', number_of_sites = 200){
  library(ComplexHeatmap)
  library(circlize)

  plot_matrix <- head(genome_annotation_info, number_of_sites)
  label_matrix = as.data.frame(ifelse(plot_matrix == 0, "No", "Yes"))

  #############################################################################################
  # Heatmap
  #############################################################################################
  # issue here, the Rscript could not directly save the pdf as planned.

  sample_list <- getSampleID(mrobj)

if (if_plot == FALSE){
  ha = HeatmapAnnotation(df = data.frame(caste = sample_list))
  ha1 = rowAnnotation(df = label_matrix['gene'], name = 'gene', show_annotation_name = TRUE, width = unit(2, "mm"), col = list(gene = c('No' =  "gray", 'Yes' = "blue")), annotation_name_side = 'top')
  ha2 = rowAnnotation(df = label_matrix['exon'], name = 'exon', show_annotation_name = TRUE, width = unit(2, "mm"), col = list(exon = c('No' =  "gray", 'Yes' = "blue")), annotation_name_side = 'top')
  ha3 = rowAnnotation(df = label_matrix['pcexon'], name = 'pcexon', show_annotation_name = TRUE, width = unit(2, "mm"), col = list(pcexon = c('No' =  "gray", 'Yes' = "blue")), annotation_name_side = 'top')
  ha4 = rowAnnotation(df = label_matrix['CDS'], name = 'CDS', show_annotation_name = TRUE, width = unit(2, "mm"), col = list(CDS = c('No' =  "gray", 'Yes' = "blue")), annotation_name_side = 'top')
  ha5 = rowAnnotation(df = label_matrix['promoter'], name = 'promoter', show_annotation_name = TRUE, width = unit(2, "mm"), col = list(promoter = c('No' =  "gray", 'Yes' = "blue")), annotation_name_side = 'top')

  heatmap <- Heatmap(plot_matrix[, 1:length(sample_list)],
          col = colorRamp2(c(0, 50, 100), c('blue4', 'blue', 'lightskyblue')),
          name = 'methylation level',
          show_row_names = FALSE,
          column_title = 'Representative methylome annotation', row_title = 'methylated sites',
          top_annotation = ha, top_annotation_height = unit(3, "mm")
  )
  draw(heatmap + ha1 + ha2 + ha3 + ha4 + ha5)
}
else{
  pdf(output_file_dir, length(sample_list) * 3, number_of_sites * 0.03)
  draw(heatmap + ha1 + ha2 + ha3 + ha4 + ha5)
  dev.off()
}
}
