#   
#   Author : Dolors Pelegrí-Sisó
#
#   The Very Large Database of Lipids (VLDL) (VLDL)
#
#
#
#   Study Record Detail :
#       https://clinicaltrials.gov/ct2/show/results/NCT01698489
#
#
# Data extracted from different studys :
#   - Deficient serum 25-hydroxy vitamin D is associated with an atherogenic lipid profile: The Very Large Database of Lipids (VLDL-3) Study
#       https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4762185/
#   - Non–High-Density Lipoprotein Cholesterol, Guideline Targets, and Population Percentiles for Secondary Prevention in 1.3 Million Adults: The VLDL-2 Study (Very Large Database of Lipids)
#       https://www.sciencedirect.com/science/article/pii/S0735109713030751?via%3Dihub
#   - Vitamin D deficiency and non-lipid biomarkers of cardiovascular risk
#       https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5510501/

# Samples : 1,310,440
#           Adults > 18 years with triglyceride levels <400 mg/dl (Estimated  using the Friedewald formula)
#                                                       >= 400 mg/dl => Excluded
#
#           Patients were 59 ± 15 years of age (mean ± SD)
#           52% were women
#
#           LDL-C cutpoints of 70, 100, 130, 160, and 190 mg/dl
#           non–HDL-C values of 93, 125, 157, 190, and 223 mg/dl,
#
#
##


library(readxl)
library(tibble)
library(dplyr)

library(BigDataStatMeth)

params_chol <- read_excel("generate_cholesterol.xlsx",
                          sheet = "PMC4762185_Transformed_R")

nSimSamples <- 1000000

params_chol <- column_to_rownames(params_chol, var = "Values")
get_pcent_sampes <- as.numeric(params_chol["n",c(3,6,9)]) / sum(as.numeric(params_chol["n",c(3,6,9)]))



lw25D <- tibble(TCholesterol = rnorm(nSimSamples*get_pcent_sampes[1] , mean = params_chol["TC","mean_lw20"], sd=params_chol["TC","sd_lw20"]),
                Age = rnorm(nSimSamples*get_pcent_sampes[1], mean = params_chol["Age","mean_lw20"], sd=params_chol["Age","sd_lw20"]),
                Insulin = rnorm(nSimSamples*get_pcent_sampes[1], mean = params_chol["Insulin","mean_lw20"], sd=params_chol["Insulin","sd_lw20"]),
                Triglycerides = rnorm(nSimSamples*get_pcent_sampes[1], mean = params_chol["Triglycerides","mean_lw20"], sd=params_chol["Triglycerides","sd_lw20"]),
                HDL_C  = rnorm(nSimSamples*get_pcent_sampes[1], mean = params_chol["HDL-C","mean_lw20"], sd=params_chol["HDL-C","sd_lw20"]),
                LDL_C  = rnorm(nSimSamples*get_pcent_sampes[1], mean = params_chol["LDL-Cd","mean_lw20"], sd=params_chol["LDL-Cd","sd_lw20"]),
                Sex  = rbinom( nSimSamples*get_pcent_sampes[1], 1, ( as.numeric( params_chol["Sex","25OHD_lw20"])/100) )
)
lw25D$D_25OHD = "<20"


bt203025D <- tibble(TCholesterol = rnorm(nSimSamples*get_pcent_sampes[2], mean = params_chol["TC","mean_2030"], sd=params_chol["TC","sd_2030"]),
                    Age = rnorm(nSimSamples*get_pcent_sampes[2], mean = params_chol["Age","mean_2030"], sd=params_chol["Age","sd_2030"]),
                    Insulin = rnorm(nSimSamples*get_pcent_sampes[2], mean = params_chol["Insulin","mean_2030"], sd=params_chol["Insulin","sd_2030"]),
                    Triglycerides = rnorm(nSimSamples*get_pcent_sampes[2], mean = params_chol["Triglycerides","mean_2030"], sd=params_chol["Triglycerides","sd_2030"]),
                    HDL_C  = rnorm(nSimSamples*get_pcent_sampes[2], mean = params_chol["HDL-C","mean_2030"], sd=params_chol["HDL-C","sd_2030"]),
                    LDL_C  = rnorm(nSimSamples*get_pcent_sampes[2], mean = params_chol["LDL-Cd","mean_2030"], sd=params_chol["LDL-Cd","sd_2030"]),
                    Sex  = rbinom( nSimSamples*get_pcent_sampes[2], 1, ( as.numeric( params_chol["Sex","25OHD_2030"])/100) )
)
bt203025D$D_25OHD = ">=20-30"


gt25D <- tibble(TCholesterol = rnorm(nSimSamples*get_pcent_sampes[3], mean = params_chol["TC","mean_gt30"], sd=params_chol["TC","sd_gt30"]),
                Age = rnorm(nSimSamples*get_pcent_sampes[3], mean = params_chol["Age","mean_gt30"], sd=params_chol["Age","sd_gt30"]),
                Insulin = rnorm(nSimSamples*get_pcent_sampes[3], mean = params_chol["Insulin","mean_gt30"], sd=params_chol["Insulin","sd_gt30"]),
                Triglycerides = rnorm(nSimSamples*get_pcent_sampes[3], mean = params_chol["Triglycerides","mean_gt30"], sd=params_chol["Triglycerides","sd_gt30"]),
                HDL_C  = rnorm(nSimSamples*get_pcent_sampes[3], mean = params_chol["HDL-C","mean_gt30"], sd=params_chol["HDL-C","sd_gt30"]),
                LDL_C  = rnorm(nSimSamples*get_pcent_sampes[3], mean = params_chol["LDL-Cd","mean_gt30"], sd=params_chol["LDL-Cd","sd_gt30"]),
                Sex  = rbinom( nSimSamples*get_pcent_sampes[3], 1, ( as.numeric( params_chol["Sex","25OHD_gt30"])/100) )
)
gt25D$D_25OHD = ">30"

all25OHd <- bind_rows(lw25D, bt203025D, gt25D) %>%
    mutate_if( is.character,as.factor)

dim(all25OHd)
head(all25OHd)


## ############################### ##
##    Create text file with data   ##
## ############################### ##

Y <- as.matrix(all25OHd$TCholesterol)
model <- model.matrix( ~ Age + Insulin  + Triglycerides + HDL_C + LDL_C, data = all25OHd)

# Save data as tab delimited text
write.table(model, "data/model.csv", quote = FALSE, sep = "\t", dec = ".", row.names = FALSE)
write.table(Y, "data/Y.csv", quote = FALSE, sep = "\t", dec = ".", row.names = FALSE)





