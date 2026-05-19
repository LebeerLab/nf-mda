// # INPUT
params.count_matrix = "data/mda.example.count_data.tsv"
params.formulas = "data/mda.example.formulas.txt"
params.metadata = "data/mda.example.meta_data.tsv"

// # OUTPUT

params.outdir = "results"
// ## MDA METHODS TO RUN
params.skip_lmclr = false



// IMPORT METHODS
include { LMCLR } from './modules/lmclr'


process VALIDATE_INPUT {

    input:
    path(counts)
    path(metadata)
    path(formulas)

    output:
    path("metadata_cleaned.tsv")
    
    script:
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(tidyverse))

    metadata <- read_table("$metadata")
    fmls <- read_lines("$formulas")
    
    mdata_fields <- lapply(fmls, as.formula) %>% 
    lapply(all.vars) %>% 
    flatten() %>% 
    unique()
    
    missing <- setdiff(mdata_fields, colnames(metadata))

    if (any(str_detect(unlist(mdata_fields), " "))) {
        stop("Please refrain from using spaces in variable names")
    }

    if (length(missing) > 0) { 
        stop(paste("Missing fields in metadata:", missing))
    } 

    counts <- read.delim("$counts")
    
    missing_samples <- setdiff(counts[['SampleID']], metadata[['SampleID']])
    if (length(missing_samples) > 0) { 
        stop(paste("Missing samples in metadata:", missing_samples))
    } 

    metadata %>%
    mutate(across(all_of(unlist(mdata_fields)), 
        function(x) str_replace_all(x, " ", ""))) %>%
    mutate(across(all_of(unlist(mdata_fields)), 
        function(x) str_replace_all(x, "[[:punct:]]", ""))) %>%
    write_tsv("metadata_cleaned.tsv")


    """

}


workflow {

  counts = Channel.fromPath(params.count_matrix, checkIfExists: true)
  metadata = Channel.fromPath(params.metadata, checkIfExists: true)
  formulas = Channel.fromPath(params.formulas, checkIfExists: true)

  VALIDATE_INPUT(counts, metadata,formulas)
  metadata_clean = VALIDATE_INPUT.out
  
  ch_results = Channel.empty()
  ch_results_full = Channel.empty()

  if (!params.skip_lmclr) {

    LMCLR(counts, metadata_clean, formulas)
    ch_results = ch_results.mix(
        LMCLR.out.res.first()
    )
    ch_results_full = ch_results_full.mix(
        LMCLR.out.res_full.first()
    )

  }

  ch_results
    .collectFile(name: 'results.tsv', storeDir: params.outdir, keepHeader: true)
  ch_results_full
    .collectFile(name: 'results_full.tsv', storeDir: params.outdir, keepHeader: true)

}
