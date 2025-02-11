suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(parallel))
##
# Parse options from the command line.
#

parseOptions <- function()
{
    defaultMarker <- "<REQUIRED>"
    # Define command line arguments
    option_list <- list(
    make_option(c("--mutationsFiles"), type = "character",
                help = "Space-separated list of RDS files for filtered mutations file"),        
    make_option(c("--locusError"), type = "character",
                help = "Space-separated list of RDS files for locus error rates"),
    make_option(c("--cosmicError"), type = "character",
                help = "Space-separated list of RDS files for cosmic error rates"),
    make_option(c("--noCosmicError"), type = "character",
                help = "Space-separated list of RDS files for non-cosmic error rates")
    )

    # Parse options
    opt_parser <- OptionParser(option_list = option_list)
    opt <- parse_args(opt_parser)
    # Utility function to split a space-separated string into a character vector of file paths
    split_files <- function(file_str) {
    file_str <- strsplit(file_str, " ")[[1]]
    if (length(file_str) == 0) {
        stop("Empty file list provided.")
    }
    else{
	return(file_str)
    }
    }
    # Get the file lists for each error type
    return(list(
	mutationsFiles = split_files(opt$mutationsFiles),
	locusError = split_files(opt$locusError),
	cosmicError = split_files(opt$cosmicError),
	noCosmicError = split_files(opt$noCosmicError)
	)
    )
}

merge_rds_files <- function(file_list, output_filename) {
  message("Merging ", length(file_list), " files into ", output_filename)
  data_list <- lapply(file_list, function(file) {
    message("Reading file: ", file)
    readRDS(file)
  })
  # Assuming that the objects can be merged with rbind (e.g. data frames or matrices)
  merged_data <- do.call(rbind, data_list)
  message("Writing merged data to ", output_filename)
  saveRDS(merged_data, file = output_filename)
}

merge_rds_lists <- function(file_list, output_filename) {
    message("Merging ", length(file_list), " files into ", output_filename)
    lst <- lapply(file_list, function(file) {
    message("Reading file: ", file)
    readRDS(file)
    })
    # Check that each file has the same names:
    commonNames <- names(lst[[1]])
    if (!all(sapply(lst, function(x) all(names(x) == commonNames)))) {
    stop("Not all files have the same list elements!")
    }
    # For each element name, merge the corresponding items from all files.
    merged_data <- lapply(commonNames, function(nm) {
    do.call(rbind, lapply(lst, function(x) x[[nm]]))
    })
    # Name the merged list
    names(merged_data) <- commonNames
    # Save the merged list if desired
    saveRDS(merged_data, file = output_filename)
}

##
# The main script, wrapped as a function.
#

main <- function(scriptArgs)
{
  print(scriptArgs$mutationsFiles)
  # Create a list with each merge job's parameters
  merge_jobs <- list(
    list(files = scriptArgs$mutationsFiles, output = "mutation_table.filtered.rds"),
    list(files = scriptArgs$locusError,     output = "locus_error_rates.off_target.rds")
  )
    # Merge each error type into a single RDS file
    # merge_rds_files(merge_jobs[[1]]$files,merge_jobs[[1]]$output)
    results <- mclapply(merge_jobs, function(job) {
    tryCatch({
        message("Merging files into ", job$output)
        merge_rds_files(job$files, job$output)
    }, error = function(e) {
        message("Error processing job with output ", job$output, ": ", e$message)
        return(NULL)
    })
    }, mc.cores = 4)
    # merge list files
    merge_rds_lists(scriptArgs$cosmicError, "error_rates.off_target.cosmic.rds")
    merge_rds_lists(scriptArgs$noCosmicError, "error_rates.off_target.no_cosmic.rds")

}

# Launch it.

invisible(main(parseOptions()))
