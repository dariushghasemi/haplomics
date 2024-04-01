install.packages("tidyverse")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("ggtext")
install.packages("config")
install.packages("haplo.stats")
install.packages("pheatmap")
install.packages("digest")
install.packages("stringr")
install.packages("DT")
install.packages("rmarkdown")

is.installed <- function(package_name){
    is.element(package_name, installed.packages()[,1])
}

# check if the package is installed
#if (!is.installed("ggtext")){
#   install.packages("ggtext");
#}
