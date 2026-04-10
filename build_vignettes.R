setwd('/Users/achazhardenberg/Library/CloudStorage/Dropbox/Repos/because')
devtools::load_all('.', quiet=TRUE)
dir.create('inst/doc', recursive=TRUE, showWarnings=FALSE)

vigs <- list.files('vignettes', pattern='\\.Rmd$', full.names=TRUE)
for (v in vigs) {
  name <- tools::file_path_sans_ext(basename(v))
  cat('Building:', name, '...\n')
  tryCatch({
    rmarkdown::render(v,
      output_dir = 'inst/doc',
      output_format = 'rmarkdown::html_vignette',
      quiet = TRUE
    )
    knitr::purl(v, output = file.path('inst/doc', paste0(name, '.R')), quiet=TRUE)
    cat('  OK\n')
  }, error = function(e) cat('  FAILED:', conditionMessage(e), '\n'))
}
cat('Done.\n')
