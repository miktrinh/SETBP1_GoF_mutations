##=============================================================================##
##              Versioning and Output Management                              ##
##=============================================================================##

#' Create versioned output directory
#' @param base_dir Base output directory
#' @param date_format Date format string
#' @param run_suffix Optional suffix for same-day runs (e.g., "_v2", "_test")
#' @param create_current_symlink Create symlink to current version (Unix only)
#' @return Full path to versioned directory
create_versioned_output_dir <- function(
  base_dir,
  date_format = "%Y_%m_%d",
  run_suffix = NULL,
  create_current_symlink = TRUE
) {
  
  # Create date-stamped folder name
  date_stamp <- format(Sys.Date(), date_format)
  
  if (!is.null(run_suffix)) {
    version_dir <- file.path(base_dir, paste0(date_stamp, run_suffix))
  } else {
    version_dir <- file.path(base_dir, date_stamp)
  }
  
  # Create directory
  if (!dir.exists(version_dir)) {
    dir.create(version_dir, recursive = TRUE)
    message(sprintf("Created versioned directory: %s", version_dir))
  } else {
    warning(sprintf("Directory already exists: %s\nFiles may be overwritten!", version_dir))
  }
  
  # Create symlink to current version (Unix-like systems only)
  if (create_current_symlink && .Platform$OS.type == "unix") {
    current_link <- file.path(base_dir, "current")
    
    # Remove old symlink if exists
    if (file.exists(current_link)) {
      unlink(current_link)
    }
    
    # Create new symlink (relative path)
    tryCatch({
      file.symlink(basename(version_dir), current_link)
      message(sprintf("Created symlink: %s -> %s", current_link, basename(version_dir)))
    }, error = function(e) {
      warning("Could not create symlink: ", e$message)
    })
  }
  
  return(version_dir)
}