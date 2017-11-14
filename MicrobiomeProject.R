#Microbiome Project

library(tidyverse) # Load the core tidyverse packages: ggplot2, tibble, 
# tidyr, readr, purrr, and dplyr

setwd("~/Documents/GitHub/DataMining_MicrobiomeProject")

microbiome <- read.csv("MicrobiomeWithMetadata.csv", encoding = 'utf-8', stringsAsFactors = FALSE)
