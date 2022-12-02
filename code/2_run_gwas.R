#Goal:run GWAS with plink2 with the available data in each consortium
soft_dir = paste0("../Software/")
data_dir = paste0("../Data/")
system("mkdir temp_result")
out_dir = "temp_result/"
#run GWAS with plink2
plink2_command = paste0(soft_dir, "plink2 ",
                        "--bfile ",data_dir,"AMR ",
                        "--pheno ",data_dir,"AMR.phen ",
                        "--covar ",data_dir,"AMR.covar ",
                        "--linear ",
                        "--out ",out_dir,"AMR")
system(plink2_command)
