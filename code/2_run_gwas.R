#Goal:run GWAS with plink2 with the available data in each consortium
soft_dir = paste0("../Software/")
data_dir = paste0("../Data/")
system("mkdir temp_result")
out_dir = "temp_result/"
eth_vec = c("AFR", "AMR", "EAS", "EUR", "SAS")
i = 2
#run GWAS with plink2
plink2_command = paste0(soft_dir, "plink2 ",
                        "--bfile ",data_dir,eth_vec[i]," ",
                        "--pheno ",data_dir,eth_vec[i],".phen ",
                        "--covar ",data_dir,eth_vec[i],".covar ",
                        "--linear ",
                        "--out ",out_dir,"AMR")
system(plink2_command)
