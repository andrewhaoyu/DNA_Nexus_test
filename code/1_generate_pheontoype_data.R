#Goal: Generate phenotypes data using 1000 Genomes data
#Randomly select 1% SNPs as causal SNPs
#Generate the effect sizes following normal distribution
#The total heritability is 0.4


#Install R packages
install.packages("data.table")
#install.packages("tidyverse")
#library(tidyverse)
library(data.table)
library(dplyr)

#copy the data from 1000 Genoems data from DCEG2 to local instance
system("cp -r /mnt/project/1000G_5Pop /home/rstudio/dnanexus/")


#Combine the SNPs from all ancestries together
data_list = list()
eth_vec = c("AFR", "AMR", "EAS", "EUR", "SAS")
i = 1
bim = fread(paste0("./dnanexus/1000G_5Pop/",
                   eth_vec[i],".bim"))




for(i in 1:5){
  bim = fread(paste0("./dnanexus/1000G_5Pop/",
                     eth_vec[i],".bim"))
}
