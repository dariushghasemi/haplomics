#!/bin/bash

# Install necessary R packages
#R -e "install.packages(c('knitr', 'rmarkdown'), repos='http://cran.us.r-project.org')"


# Check if location argument is provided
if [ -z "$1" ]; then
  echo "Usage: $0 <locus>"
  #exit 1
fi

LOCUS=$1

# Render the RMarkdown file
R -e "rmarkdown::render('04-0_report.Rmd', output_file = paste0('report_', '$LOCUS', '.html'), params = list(LOCUS = '$LOCUS'))"