#' Create a plot of the modifications for reads in a specified locus in a modBAM file using methylartist locus
#'
#' @param modBAM
#' @param region
#' @param genome_fasta
#' @param modification_motif
#' @param modifications
#' @param reads_subset
#' @param outfile
#' @param additional_options
#' @export
methylartist_locus_plot = function(modBAM, region, genome_fasta, modification_motif = "CG", 
  modifications = "m", reads_subset = NULL, outfile, additional_options = ""){
  
  # Check that methylartist is available
  if(system("methylartist --version", ignore.stdout = T, ignore.stderr = T) != 0){
    stop("methylartist does not seem to be available")
  } 
  
  # Check that region is a single genomic region and is either a GRanges or can be converted to one
  if(length(GenomicRanges::GRanges(region)) != 1){
    stop("region should be a single genomic region")
  }
  
  # If reads_subset provided, filter modBAM file for these reads
  if(!is.null(reads_subset)){
    
    # Check that samtools is available
    if(system("samtools --version", ignore.stdout = T, ignore.stderr = T) != 0){
      stop("samtools must be available if using reads_subset")
    }
    
    # Create a temporary file with read names 
    read_names_tempfile = tempfile()
    on.exit(file.remove(read_names_tempfile))
    cat(reads_subset, file = read_names_tempfile, sep = "\n")
    
    # Create a temporary sorted BAM file with specified reads and sort and index it
    tempbam = tempfile()
    on.exit(file.remove(tempbam))
    samtools_command = paste("samtools view", modBAM, "-b -N", read_names_tempfile, "| samtools sort >", tempbam)
    system(command = samtools_command)
    system(paste("samtools index", tempbam))
    
    # Update modBAM to point to the temporary BAM file with filtered reads
    modBAM = tempbam
    
  }
  
  # Create plot with methylartist locus
  methylartist_command = paste("methylartist locus", "-i", region, "-b", modBAM, "-r", genome_fasta, 
    "-n", modification_motif, "-m", modifications, "-o", outfile, additional_options)
  system(command = methylartist_command)
  
}

#' Create a plot of the modifications for reads in a specified locus in a modBAM file using methylartist region
#'
#' @param modBAM
#' @param region
#' @param genome_fasta
#' @param modification_motif
#' @param modifications
#' @param reads_subset
#' @param outfile
#' @param additional_options
#' @export
methylartist_region_plot = function(modBAM, region, genome_fasta, modification_motif = "CG", 
  modifications = "m", reads_subset = NULL, outfile, additional_options = ""){
  
  # Check that methylartist is available
  if(system("methylartist --version", ignore.stdout = T, ignore.stderr = T) != 0){
    stop("methylartist does not seem to be available")
  } 
  
  # Check that region is a single genomic region and is either a GRanges or can be converted to one
  if(length(GenomicRanges::GRanges(region)) != 1){
    stop("region should be a single genomic region")
  }
  
  # If reads_subset provided, filter modBAM file for these reads
  if(!is.null(reads_subset)){
    
    # Check that samtools is available
    if(system("samtools --version", ignore.stdout = T, ignore.stderr = T) != 0){
      stop("samtools must be available if using reads_subset")
    }
    
    # Create a temporary file with read names 
    read_names_tempfile = tempfile()
    on.exit(file.remove(read_names_tempfile))
    cat(reads_subset, file = read_names_tempfile, sep = "\n")
    
    # Create a temporary sorted BAM file with specified reads and sort and index it
    tempbam = tempfile()
    on.exit(file.remove(tempbam))
    samtools_command = paste("samtools view", modBAM, "-b -N", read_names_tempfile, "| samtools sort >", tempbam)
    system(command = samtools_command)
    system(paste("samtools index", tempbam))
    
    # Update modBAM to point to the temporary BAM file with filtered reads
    modBAM = tempbam
    
  }
  
  # Create plot with methylartist locus
  methylartist_command = paste("methylartist region", "-i", region, "-b", modBAM, "-r", genome_fasta, 
    "-n", modification_motif, "-m", modifications, "-o", outfile, additional_options)
  system(command = methylartist_command)
  
}