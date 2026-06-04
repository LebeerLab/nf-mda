// # INPUT
params.count_matrix = "data/mda.example.count_data.tsv"
params.formulas = "data/mda.example.formulas.txt"
params.metadata = "data/mda.example.meta_data.tsv"

// # OUTPUT

params.outdir = "results"
// ## MDA METHODS TO RUN
params.skip_lmclr = false
params.skip_maaslin3 = false
params.skip_maaslin2 = true
params.skip_aldex3 = false
params.skip_aldex2 = true
params.skip_limma = false
params.skip_corncob = true
params.skip_deseq2 = false
params.skip_ancomb2 = false
params.skip_zicoseq = false

// IMPORT METHODS
include { LMCLR } from './modules/lmclr'
include { MAASLIN3 } from './modules/maaslin3'
include { MAASLIN2 } from './modules/maaslin2'
include { DESEQ2 } from './modules/deseq2'
include { ALDEX3 } from './modules/aldex3'
include { ALDEX2 } from './modules/aldex2'
include { LIMMA } from './modules/limma'
include { CORNCOB } from './modules/corncob'
include { ZICOSEQ } from './modules/zicoseq'
include { ANCOMBC2 } from './modules/ancombc2'

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

    metadata_cleaned <- metadata %>%
    mutate(across(all_of(unlist(mdata_fields)), 
        function(x) str_replace_all(x, " ", ""))) %>%
    mutate(across(all_of(unlist(mdata_fields)), 
        function(x) str_replace_all(x, "[[:punct:]]", ""))) 

    if (!identical(metadata_cleaned[['SampleID']],counts[['SampleID']])) {
      # sort according to counts
      metadata_cleaned[['SampleID']] <- factor(metadata_cleaned[['SampleID']], levels=counts[['SampleID']])
      metadata_cleaned <- metadata_cleaned %>% arrange(SampleID)
    }

    metadata_cleaned %>%
    write_tsv("metadata_cleaned.tsv")

    if (!identical(metadata_cleaned[['SampleID']], factor(counts[['SampleID']]))) {
      stop(paste("Something is wrong with the order of samples in the metadata or count table", all.equal(metadata_cleaned[['SampleID']], counts[['SampleID']])))
    }


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

  if (!params.skip_maaslin3) {

    MAASLIN3(counts, metadata_clean, formulas)
    ch_results = ch_results.mix(
        MAASLIN3.out.res.first()
    )
    ch_results_full = ch_results_full.mix(
        MAASLIN3.out.res_full.first()
    )

  }

  if (!params.skip_maaslin2) {

    MAASLIN2(counts, metadata_clean, formulas)
    ch_results = ch_results.mix(
        MAASLIN2.out.res.first()
    )
    ch_results_full = ch_results_full.mix(
        MAASLIN2.out.res_full.first()
    )

  }


  if (!params.skip_limma) {

    LIMMA(counts, metadata_clean, formulas)
    ch_results = ch_results.mix(
        LIMMA.out.res.first()
    )
    ch_results_full = ch_results_full.mix(
        LIMMA.out.res_full.first()
    )

  }

  if (!params.skip_corncob) {

    CORNCOB(counts, metadata_clean, formulas)
    ch_results = ch_results.mix(
        CORNCOB.out.res.first()
    )
    ch_results_full = ch_results_full.mix(
        CORNCOB.out.res_full.first()
    )

  }

  if (!params.skip_deseq2) {

    DESEQ2(counts, metadata_clean, formulas)
    ch_results = ch_results.mix(
        DESEQ2.out.res.first()
    )
    ch_results_full = ch_results_full.mix(
        DESEQ2.out.res_full.first()
    )

  }

  if (!params.skip_aldex3) {

    ALDEX3(counts, metadata_clean, formulas)
    ch_results = ch_results.mix(
        ALDEX3.out.res.first()
    )
    ch_results_full = ch_results_full.mix(
        ALDEX3.out.res_full.first()
    )

  }

  if (!params.skip_aldex2) {

    ALDEX2(counts, metadata_clean, formulas)
    ch_results = ch_results.mix(
        ALDEX2.out.res.first()
    )
    ch_results_full = ch_results_full.mix(
        ALDEX2.out.res_full.first()
    )

  }

  if (!params.skip_ancombc2) {

    ANCOMBC2(counts, metadata_clean, formulas)
    ch_results = ch_results.mix(
        ANCOMBC2.out.res.first()
    )
    ch_results_full = ch_results_full.mix(
        ANCOMBC2.out.res_full.first()
    )

  }


  if (!params.skip_zicoseq) {

    ZICOSEQ(counts, metadata_clean, formulas)
    ch_results = ch_results.mix(
        ZICOSEQ.out.res.first()
    )
    ch_results_full = ch_results_full.mix(
        ZICOSEQ.out.res_full.first()
    )

  }


  ch_results
    .collectFile(name: 'results.tsv', storeDir: params.outdir, keepHeader: true)
  ch_results_full
    .collectFile(name: 'results_full.tsv', storeDir: params.outdir, keepHeader: true)

}
