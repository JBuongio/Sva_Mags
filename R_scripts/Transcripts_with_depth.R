# April 27th, 2023
# Aim: to understand if there are depth trends associated with MAG-specific transcription patterns in terms of number of genes turned on.

library(dplyr)
library(reshape)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(ggpubr)
library(plyr)

getwd() #"/Users/joy.buongiorno/Documents/Research/Svalbard_MAG_analysis/2021_R_Analysis/For_GitHub"
master_dat<-read.csv("Master_gene_file_TPM.csv")


#remove metagenome values
colnames(master_dat)
master_dat_trans<-master_dat[-c(5:7)]

#melt so can add site and depth information
master_dat_trans_m<-melt(master_dat_trans)
max(master_dat_trans$VKAB134_TRANSCRIPTS_MAGS_NO_RIBO)

#read in file that has variable name (sample name) alongside site and depth
variables<-read.csv("variable_to_site_depth.csv")
master_dat_trans_m_sites<-left_join(master_dat_trans_m, variables, by ="variable")

#remove rows that do not contain TPM information
master_dat_trans_m_sites_clean<-master_dat_trans_m_sites %>% drop_na()
#ensure we removed rows:
nrow(master_dat_trans_m_sites) - nrow(master_dat_trans_m_sites_clean)

#remove rows with 0 TPM
master_dat_trans_m_sites_clean_nozeros<-filter(master_dat_trans_m_sites_clean, value > 0)
#ensure we removed rows:
nrow(master_dat_trans_m_sites_clean) - nrow(master_dat_trans_m_sites_clean_nozeros)

#count
counts_all<-master_dat_trans_m_sites_clean_nozeros %>%
  dplyr::count(Prokka.Prodigal, Organism.Name, site, depth, name= "Count")

#summarize
counts_MAGs<-counts_all %>%
  dplyr::count(Organism.Name, site, depth, name = "Count")
unique(counts_MAGs$Organism.Name)

#make plots

d <- counts_MAGs

make_plot <- function(d, save_plot=FALSE, print_plot=FALSE, filename=NULL, ...) {
  p <- ggplot(d, aes(x=depth, y=Count, fill=site)) + 
    geom_point(pch=21, aes(fill=site)) + 
    geom_line(aes(color=site)) +
    scale_x_reverse() +
    coord_flip() 
  
  # Do you want to print hte plot to the screen?
  if(print_plot) {
    print(p)
  }
  
  # Do you want to save the plot?
  if(save_plot) {
    if(is.null(filename)) { # create a filename for the plot automatically, if one hasn't been specified, and add .png
      filename <- paste0(d$Organism.Name[1], ".png")
    }
    ggsave(filename, p, ...)
  }
  
  p
}

# Test this function on one MAG that I pull out manually
test_set <- d[d$Organism.Name == unique(d$Organism.Name)[1], ]
test_plot <- make_plot(test_set)
print(test_plot)

# Use dlply to make a list of each data set, and save them
plot_list <- dlply(d, c("Organism.Name"), make_plot, save_plot=TRUE, print_plot=FALSE, height=4, width=3, units="in", dpi=300)

# Later on you can do other things, like print them
l_ply(plot_list, print)
  

####Plotting cross plots of MAG mean coverage vs. number of genes transcribed

counts_MAGs$site.depth<-paste(counts_MAGs$site, counts_MAGs$depth, sep =".")

mean_cov<-read.csv("Mean_coverage_all.csv")
mean_cov_t<-melt(mean_cov)

counts_MAGs_sample<-left_join(counts_MAGs[,-c(2:3)], variables, by = "site.depth")
counts_MAGs_sample$MAG.sample<-paste(counts_MAGs_sample$Organism.Name,counts_MAGs_sample$variable, sep=".")
mean_cov_t$MAG.sample<-paste(mean_cov_t$Organism.Name, mean_cov_t$variable, sep=".")

mean_cov_with_trans_counts<-left_join(counts_MAGs_sample, mean_cov_t[,-c(1)], by = "MAG.sample")


k<- mean_cov_with_trans_counts



make_plot <- function(k, save_plot=FALSE, print_plot=FALSE, filename=NULL, ...) {
  p <- ggplot(k, aes(x=value, y=Count)) + 
    geom_point(aes(fill=depth, shape=site), size=3) + 
    scale_shape_manual(values=c(21,22,23,24)) +
    geom_smooth(method='lm') +
    stat_cor() + 
    scale_fill_gradient(high='black', low= 'yellow') +
    xlab("Mean coverage") +
    ylab("Transcript Count")
    #geom_line(aes(color=site)) +
    #scale_x_reverse() +
    #coord_flip() 
  
  # Do you want to print hte plot to the screen?
  if(print_plot) {
    print(p)
  }
  
  # Do you want to save the plot?
  if(save_plot) {
    if(is.null(filename)) { # create a filename for the plot automatically, if one hasn't been specified, and add .png
      filename <- paste0(k$Organism.Name[1], ".png")
    }
    ggsave(filename, p, ...)
  }
  
  p
}

# Test this function on one MAG that I pull out manually
test_set <- k[k$Organism.Name == unique(k$Organism.Name)[1], ]
test_plot <- make_plot(test_set)
print(test_plot)

# Use dlply to make a list of each data set, and save them
plot_list <- dlply(k, c("Organism.Name"), make_plot, save_plot=TRUE, print_plot=FALSE, height=4, width=3, units="in", dpi=300)

# Later on you can do other things, like print them
l_ply(plot_list, print)

###Want to see how sparse some MAGs are in terms of transcription



