readFlowsightStats_CC <- function(file_name) {

###############################################
#*********************************************#
# Author: Eric Prince                         #
# Date: 2017-07-17                            #
#*********************************************#
###############################################
  
# The purpose of this package is to read in the txt file generated using the Flowsight Ideas software
# for the cell cycle assay.  A table will be returned from the function call containing all raw statistics
# to be analyzed downstream.
  

# Attach dependencies
  library(tidyverse)
  library(magrittr)
  
# Bring in data file, skipping the first 3 lines because they are not of importance to the final table.
in.file <- read.table(file_name, 
                      sep = "\t",
                      skip = 3,
                      fill = TRUE,
                      stringsAsFactors = FALSE)

# Take top row as column names
names(in.file) <- in.file[1,]

# Organize the table for output
out_table <- in.file %>%
                  select(-`NA`) %>% # delete column 'NA'
                  filter(!grepl("This file was created", File)) %>% # delete trailing comment on stats table
                  slice(2:nrow(.)) %>% # start from second row because column heads and the first row are duplicates
                  mutate(File = gsub(".daf", "", File), # convert columns to numeric
                         G0G1 = as.numeric(G0G1),
                         S = as.numeric(S),
                         G2M = as.numeric(G2M))

return(out_table)
}