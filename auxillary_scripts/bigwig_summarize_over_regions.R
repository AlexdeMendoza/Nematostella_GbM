#' Summarize values in bigWig files for a set of genomic regions
#' 
#' Calculates the mean (with or without non-covered bases counting as zeros) or the sum of values in bigWig files 
#' for provided regions. Regions can be provided either as a GRanges object or as a path to a BED file.
#' Columns in the output correspond to input bigWig files and rows correspond to genomic regions. 
#' Uses names of the provided BED file or GRanges as row names so long as they are unique. 
#' Otherwise names regions region_1, region_2, etc. 
#'
#' @param bw_filepaths A vector of filepaths for bigWig files.
#' @param bed_filepath Filepath for a BED file with genomic regions of interest. 
#' One of either bed_filepath or grs must be provided. 
#' @param gr A GRanges object. One of either bed_filepath or gr must be provided. 
#' @param statistic Statistic to compute. One of "mean" (average in bigWig files over just covered bases in regions), 
#' "mean0" (average over bases with non-covered bases counting as zeroes) or "sum".
#' @param column_names A vector of names to use as column names of the output table. 
#' Default is to use the base name of bw_filepaths.
#' @param ncores Number of cores to use. Default is 1
#' @return A data.frame with the results.
#' @export
bigwig_summarize_over_regions = function(bw_filepaths, bed_filepath = NULL, 
  gr = NULL, statistic, column_names = NULL, ncores = 1){
  
  # Check that only one of bed_filepath or gr is provided
  if(!is.null(gr) & !is.null(bed_filepath)){
    stop("Either bed_filepath or gr should be provided, but not both")

  # If gr is provided add names if they are missing or are not unique
  } else if(!is.null(gr) & is.null(bed_filepath)){
    if(is.null(names(gr)) | anyDuplicated(names(gr))){
      message("names(gr) are not unique. Will thus name them region_1, region_2, etc.")
      names(gr) = paste0("region_", seq_along(gr))
    }
    
    # Remove gr metadata and export gr as a temporary BED file
    mcols(gr) = NULL
    bed_filepath = tempfile(pattern = "temp_bed")
    rtracklayer::export.bed(gr, bed_filepath)
  
  # Throw an error if neither bed_filepath or gr are provided
  } else if(is.null(gr) & is.null(bed_filepath)){
    stop("One of bed_filepath or gr must be provided")
  }
  
  # Check that one of the correct choices is provided for statistic
  match.arg(arg = statistic, choices = c("mean", "mean0", "sum"), several.ok = F)
  
  # Check that if column_names is provided, it has the same length as bw_filepaths
  if(!is.null(column_names)){if(length(bw_filepaths) != length(column_names)){stop("Length of bw_filepaths not equal to length of column_names")}}
  
  # Get names of regions from BED file and check that they are unique
  gr = rtracklayer::import.bed(bed_filepath)
  region_names = gr$name
  if(anyDuplicated(region_names)){
    message("Names of regions in BED are not unique. Will thus name them region_1, region_2, etc.")
    region_names = paste0("region_", seq_along(region_names))
    
    # Update region_names and export gr to a temporary BED file
    gr$name = region_names
    bed_filepath = tempfile(pattern = "temp_bed")
    rtracklayer::export.bed(gr, bed_filepath)
  }
  
  # Select the mean, mean0 or sum  columns from the bigWigAverageOverBed output files for downstream use
  statistic_col = ifelse(statistic == "sum", 4, ifelse(statistic == "mean0", 5, 6))
  
  # Create temporary output directory
  temp_dir = tempfile("bigwig_summarize_over_regions_tmp")
  dir.create(temp_dir)
  
  # Create names for the output files that will be produced by bigWigAverageOverBed  
  summary_over_bed_file_names = paste(temp_dir, paste0(tools::file_path_sans_ext(basename(bw_filepaths)), "_summary_over_bed.txt"), sep ="/")
  
  # If no sample names provided, set to the basename of the input bigwig-files
  if(is.null(column_names)){column_names = basename(tools::file_path_sans_ext(bw_filepaths))}
  
  # Make a cluster if more than one core being used
  if(ncores > 1){
    cluster = parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cluster, cores = ncores)
    on.exit(parallel::stopCluster(cluster))
    `%dopar%` = foreach::`%dopar%`
  } else {
    `%dopar%` = foreach::`%do%`
  }
  
  # Use bigWigAverageOverBed to calculate mean value for each bigwig file for each region in the BED file and save the results in output directory
  foreach::foreach(bigwig_number = seq_along(bw_filepaths)) %dopar% {
    system(paste(system.file("extdata", "bigWigAverageOverBed", package = "genomicTools"), bw_filepaths[bigwig_number], bed_filepath, 
      summary_over_bed_file_names[bigwig_number]), ignore.stderr = TRUE)}
  
  # Create a matrix with rows corresponding to genomic regions and columns corresponding to genomic features
  region_values = data.frame(foreach::foreach(result_file = 
   summary_over_bed_file_names, .combine = "cbind") %dopar% {
     data.frame(data.table::fread(result_file, sep = "\t", header = F, select = c(1, statistic_col), nThread = 1), row.names = 1)[region_names, ]
  })
  
  # Remove temporary files
  unlink(temp_dir, recursive = T)
  if(!is.null(gr)){
    unlink(bed_filepath)
  }
  
  # Assign rownames and column names to dataframe
  rownames(region_values) = region_names
  colnames(region_values) = column_names
  
  # Return the 
  return(region_values)
  
}
