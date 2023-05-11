## Alice Herneisen
## Last updated: 26 February 2023
## Perform thermal profiling analysis, generate melting curves, save output files

# packages
library(tidyverse)
library(TPP)

## set directories: input files (kept in inputdata folder)
setwd("~/Dropbox (MIT)/LOURIDO/COLLABORATIONS/Aarti_Krishnan/inputfiles")
file1 = "Sema_TPP_analysis_1_master_proteins_only_Proteins.txt"

## Load the protein data from Proteome Discoverer - "raw" data frames
DF1 <- read.table(file1, 
                  header = T, 
                  fill = T, 
                  stringsAsFactors = F)

## Keep most relevant columns of Proteome Discoverer data frame and add columns of 1's
## corresponding to the lowest temperature
## columns to keep from Proteome Discoverer 
DMSO_1 <- DF1 %>% dplyr::select(Description, Accession, Number.of.PSMs, Number.of.Unique.Peptides, 
                             "Abundance.Ratio.44.DMSO..40.DMSO", 
                             "Abundance.Ratio.47.DMSO..40.DMSO",
                             "Abundance.Ratio.51.DMSO..40.DMSO",
                             "Abundance.Ratio.55.DMSO..40.DMSO",
                             "Abundance.Ratio.58.DMSO..40.DMSO",
                             "Abundance.Ratio.60.DMSO..40.DMSO",
                             "Abundance.Ratio.63.DMSO..40.DMSO",
                             "Abundance.Ratio.66.DMSO..40.DMSO",
                             "Abundance.Ratio.70.DMSO..40.DMSO") %>% 
  rename(gene_name = Accession) %>%
  rename(gene_description = Description) %>% 
  add_column(Abundance_Ratio_1 = as.numeric(1), .before = 5)

SEMA_1 <- DF1 %>% dplyr::select(Description, Accession, Number.of.PSMs, Number.of.Unique.Peptides, 
                             "Abundance.Ratio.44.sema..40.sema", 
                             "Abundance.Ratio.47.sema..40.sema",
                             "Abundance.Ratio.51.sema..40.sema",
                             "Abundance.Ratio.55.sema..40.sema",
                             "Abundance.Ratio.58.sema..40.sema",
                             "Abundance.Ratio.60.sema..40.sema",
                             "Abundance.Ratio.63.sema..40.sema",
                             "Abundance.Ratio.66.sema..40.sema",
                             "Abundance.Ratio.70.sema..40.sema") %>% 
  rename(gene_name = Accession) %>%
  rename(gene_description = Description) %>% 
  add_column(Abundance_Ratio_1 = as.numeric(1), .before = 5)


## Create list of TMT labels for renaming the columns (see below)

TMTrange <- c("126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C","131")

## create data frames corresponding to input relative frequencies of the vehicle and treatment
## temperature range lysate experiments. These data frames are input to the TPP package.


names(SEMA_1)[5:14] <- paste("rel_fc_", TMTrange, sep="")
names(DMSO_1)[5:14] <- paste("rel_fc_", TMTrange, sep="")


####################  CURVE FITTING ####################

exptconfig <- read.csv("Ecoli_SEMA_config_table.csv", 
                       header = T, 
                       stringsAsFactors = F, 
                       check.names=F)
dflist <- list(DMSO_1, 
               SEMA_1)
names(dflist) <- c("DMSO_1", "SEMA_1")

resultPath <- "~/Dropbox (MIT)/LOURIDO/COLLABORATIONS/Aarti_Krishnan/exported_plots"
fits <- analyzeTPPTR(configTable = exptconfig,
                     data = dflist, 
                     resultPath = resultPath, 
                     normalize = TRUE, 
                     nCores = 1, 
                     verbose = TRUE, 
                     xlsxExport = FALSE)
                          
                        
## save the output data frames
write.table(fits, file = "DMSO_SEMA_tr_fits_R1.txt", 
                                      quote = FALSE, 
                                      sep = "\t",
                                      row.names = FALSE)

p <- ggplot(data = fits, aes(x = meltPoint_DMSO_1, 
                             y = meltPoint_SEMA_1,  
                             color = -log10(pVal_adj_SEMA_1_vs_DMSO_1),
                             text = paste0(Protein_ID, "\n",
                                           gene_description_DMSO_1))) +
                            geom_point() +
                            theme(aspect.ratio = 1)
p
ggplotly(p, tooltip = c("text"))
                          
    
p <- ggplot(data = fits, aes(x = rank(p_adj_NPARC, ties.method = "first"), 
                             y = -log10(p_adj_NPARC),
                             text = paste0(Protein_ID, "\n",
                                           gene_description_DMSO_1))) +
  geom_point(color = "snow3") +
  geom_hline(yintercept = -log10(0.01)) +
  theme(aspect.ratio = 1)
p
ggplotly(p, tooltip = c("text"))                      
