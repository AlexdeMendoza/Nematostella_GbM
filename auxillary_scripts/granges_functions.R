#' Create relative bins and count the overlaps of relative ranges with relative bins
#'
#' @param relative_ranges A GRanges object returned by relative_ranges
#' @param bin_start Start of bins
#' @param bin_end End of bins
#' @param bin_step Size of bins
#' @param category A category associated with relative_ranges
#' @return A data.frame with the the number of ranges of each category in each bin
#' @export
bin_relative_ranges = function(relative_ranges, bin_start, bin_end, bin_step, category){
  
  bin_ranges = IRanges(
    start = seq(bin_start, bin_end - bin_step, bin_step), 
    end = seq(bin_start + bin_step, bin_end, bin_step) - 1
    )
  
  if(!is.null(category)){
    ranges_list = 
    split(relative_ranges, category)
  } else {ranges_list = list(relative_ranges)}
  
  overlap_df = data.frame(lapply(ranges_list, function(x)
    setNames(countOverlaps(bin_ranges, ranges(x)), start(bin_ranges) + width(bin_ranges)/2)))
  
  overlap_df = tibble::rownames_to_column(overlap_df, "bin_center")
  overlap_df$bin_center = as.numeric(overlap_df$bin_center)
  
  return(overlap_df)

}
