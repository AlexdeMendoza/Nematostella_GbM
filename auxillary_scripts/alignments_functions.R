#' Filter a BAM file for alignments to specified regions 
#'
#' @param bam_file
#' @param genomic_regions
#' @param param A scanBamParam object specifying the components of the BAM records to be returned. 
#' Default is the same as scanBam.
#' @param nthreads
#' @return A list returned by running scanBam on the BAM file created by filtering for specified reads
#' @export
filter_bam_for_regions = function(bam_file, genomic_regions, param = NULL, nthreads = 4){
  
  # Check that samtools is available
  if(system("samtools --version", ignore.stdout = T) != 0){
    stop("samtools does not seem to be available")
  }
  
  # Export a temporary BED file with regions in genomic_regions
  S4Vectors::mcols(genomic_regions) = NULL
  temporary_bed = tempfile(pattern = "temp_bed_")
  rtracklayer::export.bed(genomic_regions, temporary_bed)
  
  # Set param to default settings if it is not set
  if(is.null(param)){
    param = ScanBamParam(what=scanBamWhat())
  }
  
  # Filter input BAM for reads
  tempbam = tempfile()
  samtools_command = paste("samtools view", bam_file, "-b -L", temporary_bed, "-@", nthreads, "-o", tempbam)
  system(command = samtools_command)
  on.exit(file.remove(tempbam))
  
  # Index BAM file
  system(paste("samtools index", "-@", nthreads, tempbam))
  on.exit(file.remove(tempbam, paste0(tempbam, ".bai")))
  
  # Return list with the results
  return(Rsamtools::scanBam(file = tempbam, param = param))
  
}
