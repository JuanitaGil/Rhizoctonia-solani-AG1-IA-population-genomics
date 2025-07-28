## Install and load required libraries

library(viridis)
library(ggplot2)
#if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install("LEA")
library(LEA)
library(tidyr)



## sNMF requires a geno file. Convert filtered, clone-corrected vcf to geno file using the library LEA
geno_cc <- vcf2geno("/path/to/pop_cc.vcf", output.file = "pop.geno", force = TRUE)

#create project object snmf (genofile, K=the K or range of Ks to test, repetitions=number of repetitions, entropy= TRUE to calculate the cross entropy, project="new" for new project)
# Firs look at the overall cross entropy. No need to run repetitions
project_pop.snmf <- snmf(geno_cc,
                        K = 1:10,
                        #repetitions = 10,
                        entropy = TRUE,
                        project = "new", CPU = 8)

#To see the K with the less entropy, and therefore the best one:
plot(project_cc.snmf, col = "#21918c", pch = 19, cex = 1.2)

#Now run snmf again with 10 reps per K
project_cc_10reps.snmf <- snmf(geno_cc,
                          K = 1:10,
                          repetitions = 10,
                          entropy = TRUE,
                          project = "new", CPU = 8)

#To see the K with the less entropy, and therefore the best one:
plot(project_cc_10reps.snmf, col = "#21918c", pch = 19, cex = 1.2)

# get the cross-entropy for all runs for K = 1 to K =10 and then get best (lowest) entropy
# This is to select the best run per K
# Initialize vectors to store cross-entropy and best indices
#ces <- numeric(10)
#best_ces <- numeric(10)

# Loop through K values and compute cross-entropy
#for (K in 1:10) {
#  ces[K] <- cross.entropy(project_cc_10reps.snmf, K = K)
#  best_ces[K] <- which.min(ces[K])
#}

ces_1 <- cross.entropy(project_cc_10reps.snmf, K = 1)
best_ces_1 <- which.min(ces_1)
ces_2 <- cross.entropy(project_cc_10reps.snmf, K = 2)
best_ces_2 <- which.min(ces_2)
ces_3 <- cross.entropy(project_cc_10reps.snmf, K = 3)
best_ces_3 <- which.min(ces_3)
ces_4 <- cross.entropy(project_cc_10reps.snmf, K = 4)
best_ces_4 <- which.min(ces_4)
ces_5 <- cross.entropy(project_cc_10reps.snmf, K = 5)
best_ces_5 <- which.min(ces_5)
ces_6 <- cross.entropy(project_cc_10reps.snmf, K = 6)
best_ces_6 <- which.min(ces_6)
ces_7 <- cross.entropy(project_cc_10reps.snmf, K = 7)
best_ces_7 <- which.min(ces_7)
ces_8 <- cross.entropy(project_cc_10reps.snmf, K = 8)
best_ces_8 <- which.min(ces_8)
ces_9 <- cross.entropy(project_cc_10reps.snmf, K = 9)
best_ces_9 <- which.min(ces_9)
ces_10  <- cross.entropy(project_cc_10reps.snmf, K = 10)
best_ces_10 <- which.min(ces_10)
