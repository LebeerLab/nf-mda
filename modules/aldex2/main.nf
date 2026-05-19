process ALDEX2 {

  container "thiesgehrmann/aldex2"

  input:
  path(counts)
  path(metadata)
  path(formulas)

  output:
  path('res.tsv'), emit: res
  path('res_full.tsv'), emit: res_full
  script:
  """
  #!/usr/bin/env Rscript
  suppressPackageStartupMessages(library(tidyverse))
  library(multidiffabundance)

  counts <- read_delim("$counts") %>% column_to_rownames("SampleID")
  metadata <- read_tsv("$metadata") %>% column_to_rownames("SampleID")
  formulas <- read_lines("$formulas")

  D <- mda.create(
                counts,
                metadata,  
                formulas,                   
                usecache=FALSE,
                recache=FALSE)
  res <- mda.aldex2(D)
  res[['res.full']] %>% write_tsv("res_full.tsv")
  res[['res']] %>% write_tsv("res.tsv")
  
  """
}
