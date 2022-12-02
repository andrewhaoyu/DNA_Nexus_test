#Goal: Generate phenotypes data using 1000 Genomes data
#Randomly select 1% SNPs as causal SNPs
#Generate the effect sizes following normal distribution
#The total heritability is 0.4


#Install R packages
install.packages("data.table")
install.packages("dplyr")
#install.packages("tidyverse")
#library(tidyverse)
library(data.table)
library(dplyr)

#copy the data from 1000 Genoems data from DCEG2 to local instance
#system("cp -r /mnt/project/1000G_5Pop /home/rstudio/dnanexus/")


#Combine the SNPs from all ancestries together
data_list = list()
eth_vec = c("AFR", "AMR", "EAS", "EUR", "SAS")


for(i in 1:5){
  bim = fread(paste0("../dnanexus/1000G_5Pop/",
                     eth_vec[i],".bim"))
  data_list[[i]] = bim
}
#combine the bim file together
data_com = rbindlist(data_list)
#find the unique SNPs
data_com = unique(data_com)

#randomly select 1% causal SNPs
n_snp = nrow(data_com)
n_cau_snp = as.integer(0.01*n_snp)
cau_id = sample(c(1:n_snp), n_cau_snp, replace = F)
data_com_cau = data_com[cau_id, ]
h2 = 0.4
#for now, I am putting the effect-size for all ancestries to be the same
beta = rnorm(n_cau_snp, 0, sd = sqrt(h2/n_cau_snp))
#We only need a subset of SNPs to generate the phenotype
#First use PLINK1.9 to extract the causal SNPs for each ancestry
cau_snp_id = data_com_cau$V2
system("mkdir temp_result")
write.table(cau_snp_id, 
            file = "./temp_result/cau_snp_id", 
            row.names = F,
            col.names = F, 
            quote = F)
soft_dir = "../dnanexus/Software/"
data_dir = "../dnanexus/1000G_5Pop/"
temp_dir = "./temp_result/"
for(i in 1:5){
  plink_command = paste0(soft_dir, "plink2 ",
                         "--bfile ",data_dir, eth_vec[i]," ",
                         "--extract ",temp_dir,"cau_snp_id ",
                         "--make-bed ",
                         "--out ",temp_dir,"cau_",eth_vec[i])
  system(plink_command)
  
}
#use GCTA to generate phenotypes data
#prepare the effect size table
effect_size = data.frame(snp = cau_snp_id, effect_size = beta)
write.table(effect_size, 
            file = paste0(temp_dir, "effect_size"),
            row.names = F,
            col.names = F,
            quote = F)
for(i in 1:5){
  gcta_command = paste0(soft_dir, "gcta ",
                        "--bfile ",temp_dir,"cau_",eth_vec[i]," ",
                        "--simu-qt ",
                        "--simu-causal-loci ",paste0(temp_dir, "effect_size "),
                        "--simu-hsq ",h2," ",
                        "--simu-rep 1 ",
                        "--out ",paste0(temp_dir, eth_vec[i])
  )
  system(gcta_command)
}

system("mkdir 1000G_5Pop_simulate_phenotypes")
out_dir = "1000G_5Pop_simulate_phenotypes/"
#move the result from temp to 1000G_5Pop_simulate_phenotypes
for(i in 1:5){
  mv_command = paste0("mv ",temp_dir,eth_vec[i],".phen ",
                      out_dir)
  system(mv_command)
}
#simulate the covariates data
for(i in 1:5){
  fam = read.table(paste0(data_dir, eth_vec[i],".fam"))
  
  covar = data.frame(V1 = fam$V1,
                     V2 = fam$V2,
                     covar = rnorm(nrow(fam), 0, 1))
  write.table(covar, file = paste0(out_dir, eth_vec[i],".covar"))  
}
#upload the data from out_dir to DNANexus project in ttyd
