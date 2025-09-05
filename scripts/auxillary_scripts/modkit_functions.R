#' Create a GRanges object from the output of modkit extract calls
#'
#' @param mod_calls_table A data.frame with lines from a table output by modkit extract calls.
#' @param remove_mcols A logical value indicating whether to return a GRanges without any metadata columns added. 
#' Can be useful if mod_calls_table is very large. Default is FALSE.
#' @return A GRanges object with the location of modification calls in the input table.
#' @export
make_granges_from_mod_calls = function(mod_calls_table, remove_mcols = FALSE){
  
  # Add a column with the end position of modifications 
  mod_calls_table$end_position = mod_calls_table$ref_position + 1
  
  # Create GRanges with modification positions and return
  return(GenomicRanges::makeGRangesFromDataFrame(mod_calls_table, seqnames.field = "chrom", 
    start.field = "ref_position", end.field = "end_position", strand.field = "ref_strand", 
    keep.extra.columns = !remove_mcols, starts.in.df.are.0based = T))
  
}

#' Extract modified base calls from a modBAM for portions of reads overlapping specified genomic regions
#'
#' @param modBAM Path to a modBAM file.
#' @param genomic_regions A GRanges with genomic regions from which to extract modified bases. 
#' @param add_genomic_region_index A logical value indicating whether to add column to the output with the index of 
#' the region in genomic_regions that the mod call overlaps. Ranges in genomic_regions should not overlap if add_genomic_region_index is set. Default is FALSE.
#' @param pass_only A logical value indicating whether to only output modifications passing the call threshold. 
#' Default is TRUE.
#' @param motif An optional sequence motif to extract modifications at. Default is CG. Set to NULL to extract modifications at all positions.
#' @param motif_offset The 0-based offset in the motif for the base to extract modifications at. Default is 0. 
#' @param nthreads Number of threads to use with modkit extract calls. Default is 4.
#' @param reference_fasta Optional path to a FASTA file for the relevant reference genome for the modBAM file for 
#' modkit to extract reference context information from.
#' @return A data.table with the modification calls for the portions of reads overlapping genomic_regions
#' @export
extract_mod_calls = function(modBAM, genomic_regions, add_genomic_region_index = FALSE, 
  pass_only = TRUE, motif = "CG", motif_offset = 0, nthreads = 4, reference_fasta = NULL){
  
  ### Perform checks ###
  
  # Check that modkit is available
  if(system("modkit --version", ignore.stdout = T, ignore.stderr = T) != 0){
    stop("modkit does not seem to be available")
  } 

  # Check that genomic_regions is a GRanges object and that it does not containing overlapping ranges
  if(!is(genomic_regions, "GRanges")){stop(paste("genomic_regions must be a GRanges object"))}
  if(!isDisjoint(genomic_regions) & add_genomic_region_index){
    stop(paste("genomic_regions cannot contain overlapping ranges when add_genomic_region_index is set. Use reduce() to merge overlapping regions."))
  }
  
  # Give a warning if no reference provided
  if(is.null(reference_fasta)){warning("Reference context information will be missing as no reference fasta provided")}
  
  ### Run modkit  ###
  
  # Check is pass_only is set
  if(pass_only){
    pass_string = "--pass-only"
  } else{
    pass_string = ""
  }
  
  # Set motif if it is provided
  if(!is.null(motif)){
    motif_string = paste("--motif", motif, motif_offset)
  } else {
    motif_string = ""
  }
  
  # Export genomic_regions as a BED file
  tempbed = tempfile()
  rtracklayer::export.bed(genomic_regions, tempbed)
  bed_string = paste("--include-bed", tempbed)
  
  # Run modkit with specified number of cores on input modBAM
  temp_modkit_file = tempfile()
  modkit_command = paste("modkit extract calls", "-t", nthreads, "--reference", reference_fasta, 
    pass_string, motif_string, bed_string, modBAM, temp_modkit_file)
  message("Extracting base modification calls with modkit")
  system(command = modkit_command)
  
  # Read temp_modkit_file
  mod_calls_df = data.table::fread(temp_modkit_file, nThread = nthreads, data.table = FALSE)
  
  # Add column indicating which region each mod call overlaps if add_genomic_region_index is set
  if(add_genomic_region_index){
    mod_calls_gr = make_granges_from_mod_calls(mod_calls_df, remove_mcols = T)
    mod_calls_df$genomic_region_index = data.frame(findOverlaps(mod_calls_gr, genomic_regions))$subjectHits
  }
  
  # Return table
  return(mod_calls_df)
  
}

#' Filter a table produced by modkit extract calls for modifications meeting certain conditions
#'
#' @param modkit_table A data.frame with the output of modkit extract calls.
#' @param mapped A logical value indicating whether to filter for mapped positions. Default is TRUE.
#' @param mod_call_prob_threshold Threshold value to filter modification calls (call_code values besides "-") using call_prob. Default value is the 
#' minimum passing call_prob threshold from modkit_table. 
#' @param unmod_call_prob_threshold Threshold value to filter unmodified base calls (call_code of -) using call_prob. Default value is the 
#' minimum passing call_prob threshold from modkit_table. 
#' @param base_call_threshold Threshold value to filter using base_qual. Default value is 0.
#' @param sequence_context An optional character vector giving one or more sequence contexts in the reads to 
#' extract modifications for e.g. "CG". All provided sequence contexts must be of the same length.
#' @param mod_base_offset The offset in query_kmer to the base of interest. 
#' Defaults to the middle of query_kmer e.g. the 3rd position if query_kmer is a 5-mer.
#' @param sequence_context_matches_ref A logical value indicating to filter for modifications where the sequence context of 
#' the reads matches the reference sequence e.g. a CG site in the reads aligns to a CG site in the reference. 
#' Matches against the reverse complement of the reference for reads which align to the - strand. Default value is FALSE.
#' @return A data.frame with the table resulting from applying the specified filters to modkit_table.
#' @export
filter_modkit_calls = function(modkit_table, mapped = TRUE, mod_call_prob_threshold = NULL, unmod_call_prob_threshold = NULL, 
  base_qual_threshold = 20, sequence_context = "CG", mod_base_offset = NULL, sequence_context_matches_ref = TRUE){
  
  # Check that all sequence contexts have the same length if sequence_context provided
  if(!is.null(sequence_context)){
    if(length(unique(nchar(sequence_context))) > 1){
      stop("All provided sequence contexts should be of the same length")
    }
  }
  
  # Identify the minimum passing call_prob threshold from modkit_table if none provided
  if(is.null(mod_call_prob_threshold)){
    mod_call_prob_threshold = min(dplyr::filter(modkit_table, !fail)$call_prob)
  }
  if(is.null(unmod_call_prob_threshold)){
    unmod_call_prob_threshold = min(dplyr::filter(modkit_table, !fail)$call_prob)
  }
  
  # Filter for modification calls where call_prob exceeds the mod_call_prob_threshold and
  # unmodified base calls where call_prob exceeds the unmod_call_prob_threshold
  modkit_table = rbind(
    dplyr::filter(modkit_table, call_code == "-" & call_prob > mod_call_prob_threshold),
    dplyr::filter(modkit_table, call_code != "-" & call_prob > unmod_call_prob_threshold)
  )
  
  # Filter for mapped positions if specified
  if(mapped){modkit_table = dplyr::filter(modkit_table, ref_position > 0, chrom != ".")}
  
  # Filter using base_qual_threshold
  modkit_table = dplyr::filter(modkit_table, base_qual > base_qual_threshold)
  
  # If mod_base_offset is not provided, it is assumed to be the central base in the Kmer
  kmer_length = nchar(modkit_table$query_kmer)[1]
  if(is.null(mod_base_offset)){mod_base_offset = (kmer_length+1)/2}
  
  # If sequence_context provided, filter for modifications occurring in this context. 
  if(!is.null(sequence_context)){
    modkit_table = dplyr::mutate(modkit_table, 
      query_context = substr(query_kmer, mod_base_offset, mod_base_offset + nchar(sequence_context[1]) - 1))
    modkit_table = dplyr::filter(modkit_table, query_context %in% sequence_context)
    
    # If sequence_context_matches_ref is set, filter for modifications where the query sequence context matches the 
    if(sequence_context_matches_ref){
      
      # Convert ref_kmer to a DNAStringSet and get the reverse complement for sequences on the - strand
      ref_kmer_dss = Biostrings::DNAStringSet(modkit_table$ref_kmer)
      ref_kmer_dss[modkit_table$ref_strand == "-"] = 
        Biostrings::reverseComplement(ref_kmer_dss[modkit_table$ref_strand == "-"])
      ref_kmer_dss = as.character(ref_kmer_dss)
      
      # Filter for modifications where the query sequence context matches the reference sequence context
      modkit_table = dplyr::mutate(modkit_table, 
        ref_context = substr(ref_kmer_dss, mod_base_offset, mod_base_offset + nchar(sequence_context[1]) - 1))
      modkit_table = dplyr::filter(modkit_table, query_context == ref_context)
    }
  }
  
  # Return the final table
  return(modkit_table)
}

#' Summarize modification calls per read
#'
#' @param mod_calls_table A data.frame with lines from a table output by modkit extract calls.
#' @return A data.frame with the per-read proportion of modification calls for each type of modification detected. 
#' @export
summarize_mod_calls_per_read = function(mod_calls_table){
  
  # Group by read_id, keeping any existing groups
  mod_calls_table = dplyr::group_by(mod_calls_table, read_id, .add = TRUE)
  
  # Count occurrences of modification calls for each group
  mod_calls_table_summary = dplyr::count(mod_calls_table, call_code)
  
  # Get proportion of modification calls for each group and remove count column
  mod_calls_table_summary = dplyr::mutate(mod_calls_table_summary, proportion = n/sum(n), n = NULL)
  
  # Convert mod_calls_table_summary to wide format and replace any missing modification values with 0
  mod_calls_table_summary = tidyr::pivot_wider(mod_calls_table_summary, 
    names_from = call_code, values_from = proportion, values_fill = 0)
  return(setNames(data.frame(mod_calls_table_summary), names(mod_calls_table_summary)))
  
}
