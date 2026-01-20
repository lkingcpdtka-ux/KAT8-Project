#!/usr/bin/env Rscript
################################################################################
## KAT8 Knockdown RNA-seq Analysis Pipeline - Master Script (CELLS)
##
## This script runs all 3 parts of the cell analysis in sequence:
##   1. part1_cells_main_analysis.R  - DESeq2 differential expression analysis
##   2. part2_cells_ora.R            - Over-representation analysis (ORA)
##   3. part3_cells_fgsea.R          - Gene set enrichment analysis (fgsea)
##
## USAGE:
##   Rscript run_all_cells.R
##
##   OR from R console:
##   source("run_all_cells.R")
##
## REQUIREMENTS:
##   - All input data files must be in the expected locations
##   - Required R packages must be installed (see individual scripts)
##   - Sufficient memory for DESeq2 analysis (recommended: 8+ GB RAM)
##
## OUTPUT:
##   - All outputs will be saved in the standard output directories
##   - Log files will be saved to output/logs/
##   - Session info will be saved for each part for reproducibility
##
## Author: KAT8 Project Team
## Last Updated: 2026-01-20
################################################################################

## Set up logging
cat("\n")
cat("================================================================================\n")
cat("=== KAT8 KNOCKDOWN RNA-SEQ ANALYSIS PIPELINE (CELLS) ===\n")
cat("================================================================================\n")
cat("Start time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
cat("\n")

## Track timing for each part
pipeline_start_time <- Sys.time()
timing_log <- data.frame(
  Part = character(),
  Start_Time = character(),
  End_Time = character(),
  Duration_Minutes = numeric(),
  Status = character(),
  stringsAsFactors = FALSE
)

## Helper function to run a script with error handling
run_part <- function(script_name, part_number, description) {
  cat("\n")
  cat("================================================================================\n")
  cat("=== PART ", part_number, ": ", toupper(description), " ===\n", sep = "")
  cat("================================================================================\n")
  cat("Script: ", script_name, "\n", sep = "")
  cat("Start time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
  cat("\n")

  part_start <- Sys.time()

  ## Check if script exists
  if (!file.exists(script_name)) {
    cat("[ERROR] Script not found: ", script_name, "\n", sep = "")
    stop("Cannot continue pipeline - missing script file")
  }

  ## Run the script
  tryCatch({
    source(script_name, echo = FALSE)
    part_end <- Sys.time()
    duration <- as.numeric(difftime(part_end, part_start, units = "mins"))

    cat("\n")
    cat("[SUCCESS] Part ", part_number, " completed successfully\n", sep = "")
    cat("Duration: ", round(duration, 2), " minutes\n", sep = "")
    cat("End time: ", format(part_end, "%Y-%m-%d %H:%M:%S"), "\n", sep = "")

    ## Log timing
    timing_log <<- rbind(
      timing_log,
      data.frame(
        Part = paste0("Part ", part_number, ": ", description),
        Start_Time = format(part_start, "%Y-%m-%d %H:%M:%S"),
        End_Time = format(part_end, "%Y-%m-%d %H:%M:%S"),
        Duration_Minutes = round(duration, 2),
        Status = "SUCCESS",
        stringsAsFactors = FALSE
      )
    )

  }, error = function(e) {
    part_end <- Sys.time()
    duration <- as.numeric(difftime(part_end, part_start, units = "mins"))

    cat("\n")
    cat("================================================================================\n")
    cat("[ERROR] Part ", part_number, " FAILED\n", sep = "")
    cat("================================================================================\n")
    cat("Error message: ", conditionMessage(e), "\n", sep = "")
    cat("Duration before failure: ", round(duration, 2), " minutes\n", sep = "")
    cat("\n")

    ## Log timing
    timing_log <<- rbind(
      timing_log,
      data.frame(
        Part = paste0("Part ", part_number, ": ", description),
        Start_Time = format(part_start, "%Y-%m-%d %H:%M:%S"),
        End_Time = format(part_end, "%Y-%m-%d %H:%M:%S"),
        Duration_Minutes = round(duration, 2),
        Status = "FAILED",
        stringsAsFactors = FALSE
      )
    )

    ## Stop pipeline
    stop("Pipeline halted due to error in Part ", part_number)
  })
}

################################################################################
## RUN PIPELINE
################################################################################

## Part 1: DESeq2 differential expression analysis
run_part("part1_cells_main_analysis.R", 1, "DESeq2 Analysis (Cells)")

## Part 2: Over-representation analysis
run_part("part2_cells_ora.R", 2, "ORA Analysis (Cells)")

## Part 3: Gene set enrichment analysis
run_part("part3_cells_fgsea.R", 3, "FGSEA Analysis (Cells)")

################################################################################
## PIPELINE COMPLETE
################################################################################

pipeline_end_time <- Sys.time()
total_duration <- as.numeric(difftime(pipeline_end_time, pipeline_start_time, units = "mins"))

cat("\n")
cat("================================================================================\n")
cat("=== PIPELINE COMPLETE ===\n")
cat("================================================================================\n")
cat("Total duration: ", round(total_duration, 2), " minutes (",
    round(total_duration / 60, 2), " hours)\n", sep = "")
cat("End time: ", format(pipeline_end_time, "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
cat("\n")

## Print timing summary
cat("--- Timing Summary ---\n")
print(timing_log, row.names = FALSE)
cat("\n")

## Save timing log
timing_file <- file.path("output", "logs",
                        paste0("pipeline_timing_cells_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"))

## Create logs directory if it doesn't exist
if (!dir.exists(file.path("output", "logs"))) {
  dir.create(file.path("output", "logs"), recursive = TRUE)
}

write.csv(timing_log, file = timing_file, row.names = FALSE)
cat("[INFO] Timing log saved to: ", timing_file, "\n", sep = "")

cat("\n")
cat("All cell analysis complete! Check the output directory for results.\n")
cat("\n")
