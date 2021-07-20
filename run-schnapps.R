read_haplotypes_dlp <- function(paths,
                               cols_to_keep = c("cell_id", "chr", "start", "end", "allele_id", "hap_label", "readcount")){
  
  hapslist <- list()
  for (x in paths){
    hapslist[[x]] <- data.table::fread(x, 
                                       colClasses = list(character = c("chromosome", "cell_id"), 
                                                         integer = c("start", "end", "hap_label", "readcount", "allele_id")))
  }
  
  hapsdata <- data.table::rbindlist(hapslist, use.names = T)
  hapsdata <- dplyr::rename(hapsdata, chr = chromosome)
  hapsdata <- hapsdata[, ..cols_to_keep]
  
  return(hapsdata)
}

read_copynumber_dlp <- function(cnpaths, 
                               metricspaths, 
                               sample_ids = NULL,
                               patient = NULL,
                               filtercells = TRUE, 
                               cols_to_keep = c("cell_id", "chr", "start", "end", "map", "copy", "state"),
                               mappability = 0.99,
                               s_phase_filter = FALSE,
                               quality_filter = 0.75,
                               filter_reads = 0.1e6){
  
  cnlist <- list()
  for (x in cnpaths){
    cnlist[[x]] <- data.table::fread(x, 
                                     colClasses = list(character = c("chr", "cell_id"), 
                                                       integer = c("start", "end", "reads", "state"),
                                                       numeric = c("map", "gc", "copy")))
  }
  
  cndata <- data.table::rbindlist(cnlist, use.names = T)
  cndata <- cndata[, ..cols_to_keep]
  cndata <- cndata[map > mappability]
  
  metricslist <- list()
  i <- 1
  for (x in metricspaths){
    if (!is.null(sample_ids)) {
      metricslist[[x]] <- data.table::fread(x)[, sample_id := sample_ids[i]]
      i <- i + 1
    } else {
      metricslist[[x]] <- data.table::fread(x)
    }
  }
  
  metricsdata <- data.table::rbindlist(metricslist, use.names = T)
  
  if (filtercells == TRUE){
    
    message(paste0("Total number of cell pre-filtering: ", dim(metricsdata)[1]))
    
    cells_to_keep <- metricsdata %>% 
      .[quality > quality_filter] %>% 
      .[is_contaminated == FALSE] %>% 
      .[is_s_phase == FALSE] %>%
      .[total_mapped_reads > filter_reads] %>% 
      .[!stringr::str_detect(experimental_condition, "NTC|NCC|gDNA|GM|CONTROL|control")]
    
    if (s_phase_filter){
      cells_to_keep <- cells_to_keep[is_s_phase == FALSE]
    }
    
    metricsdata <- metricsdata[cell_id %in% cells_to_keep$cell_id]
    message(paste0("Total number of cell post-filtering: ", dim(metricsdata)[1]))
    cndata <- cndata[cell_id %in% cells_to_keep$cell_id]
  }
  
  #cols <- c("cell_id", "sample_id")
  #cndata <- cndata[metricsdata[, ..cols], on = c("cell_id")]
  
  return(list(cn = cndata %>% as.data.frame(), metrics = metricsdata %>% as.data.frame()))
}

callhscn <- function(cn, 
                     qc,
                     haps,
                     ncores = 1){
  
  cn <- as.data.table(cn)
  qc <- as.data.table(qc)
  haps <- as.data.table(haps)
  
  message(paste0("Total number of cells: ", dim(qc)[1]))
  
  cn_ploidy <- cn[, list(ploidy = mean(state), var = var(state), frac_nondiploid = sum(state != 2) / .N), by = "cell_id"]
  
  non_diploidcells <- cn_ploidy[frac_nondiploid > 0.1]
  cn <- cn[cell_id %in% non_diploidcells$cell_id]
  
  message(paste0("Total number of non-diploid cells: ", dim(qc)[1]))
  
  if (length(unique(cn$cell_id)) == 0){
    message("No cells returning NULL")
    return(NULL)
  }
  
  message("Format haplotypes file")
  haps <- haps[cell_id %in% qc$cell_id]
  haps <- format_haplotypes_dlp(haps, cn)
  
  message("Infer HSCN")
  hscn <- callHaplotypeSpecificCN(cn, 
                                  haps, 
                                  likelihood = "auto",
                                  minfrac = 0.7,
                                  cluster_per_chr = TRUE, 
                                  ncores = ncores)
  return(hscn)
}

library(argparse)
library(tidyverse)
library(data.table)
library(schnapps)

parser <- ArgumentParser()

parser$add_argument("--hmmcopyreads", default="character", nargs = "+",
                    help = "hmmcopy reads files")
parser$add_argument("--hmmcopyqc", default=NULL, type="character", nargs = "+",
                    help="hmmcopy QC files")
parser$add_argument("--allelecounts", default=NULL, type="character", nargs = "+",
                    help="hmmcopy QC files")
parser$add_argument("--ncores", default=NULL, type="integer", 
                    help="Number of cores in schnapps inference")
parser$add_argument("--csvfile", default=NULL, type="character", 
                    help="output csvfile")
parser$add_argument("--Rdatafile", default=NULL, type="character", 
                    help="output Rdata file")
parser$add_argument("--heatmap", default=NULL, type="character", 
                    help="heatmap plot")
parser$add_argument("--qcplot", default=NULL, type="character", 
                    help="QC plot")

args <- parser$parse_args()

message("Read in copy number reads and qc")
cndata <- read_copynumber_dlp(cnpaths = args$hmmcopyreads,
                              metricspaths = args$hmmcopyqc)

message("Read in haplotypes data")
haplotypes <- read_haplotypes_dlp(args$allelecounts)

message("Call haplotype specific copy number")
res <- callhscn(cndata$cn, 
                cndata$metrics, 
                haplotypes,
                ncores = args$ncores)

message("Write output files")
saveRDS(file = args$Rdatafile, object = res)
fwrite(x = res$data, file = args$csvfile)


message("Make QC plots")
g1 <- plotBAFperstate(res)
g2 <- plot_variance_state(res)
cowplot::save_plot(args$qcplot,cowplot::plot_grid(g1, g2), 
                   base_width = 10, base_height = 5)

message("Make heatmap")
cl <- umap_clustering(res$data, field = "copy")
h1 <- plotHeatmap(res, 
                  tree = cl$tree, 
                  reorderclusters = T,
                  clusters = cl$clustering,
                  plottree = F,
                  plotcol = "state")
h2 <- plotHeatmap(res, 
                  tree = cl$tree, 
                  reorderclusters = T,
                  clusters = cl$clustering,
                  plottree = F,
                  plotcol = "state_phase", show_clone_label = F, show_library_label = F)


hm1 <- grid::grid.grabExpr(ComplexHeatmap::draw(h1 + h2,  
                                                ht_gap = unit(0.6, "cm"),
                                                column_title = paste0("Number of cells: ", length(unique(res$data$cell_id))), 
                                                column_title_gp = grid::gpar(fontsize = 20),
                                                heatmap_legend_side = "bottom", 
                                                annotation_legend_side = "bottom",
                                                show_heatmap_legend = TRUE), 
                                                width = 40, height = 13/3)

cowplot::save_plot(args$heatmap, 
                   cowplot::plot_grid(hm1), 
                   base_height = 10, 
                   base_width = 20)

message("Finished")

