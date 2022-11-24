# codes used to run the corHMM function on each dataset version ----

require(corHMM)

# i in the file names varies from 1 to 6 and stands for the dataset version

phy <- readRDS("Phylogenetic tree of dataset vi.RDS")

dat <- read.csv("Data frame of dataset vi.csv")




# models of independent evolution of drought adaptation ----

Dro_dat <- dat[, c("Plant_species", "Drought_adaptation")]

# with ER model structure, j hidden rates categories, and replicated start x

indepDroERj_x <- corHMM(phy, Dro_dat, rate.cat = j, model = ER) # j varies from 1 to 3 and x from 1 to 5

saveRDS(indepDroERj_x, "indepDroERj_x.RDS")

# with SYM model structure, j hidden rates categories, and replicated start x

indepDroSYMj_x <- corHMM(phy, Dro_dat, rate.cat = j, model = SYM) # j varies from 1 to 3 and x from 1 to 5

saveRDS(indepDroSYMj_x, "indepDroSYMj_x.RDS")

# with ARD model structure, j hidden rates categories, and replicated start x

indepDroARDj_x <- corHMM(phy, Dro_dat, rate.cat = j) # j varies from 1 to 3 and x from 1 to 5

saveRDS(indepDroARDj_x, "indepDroARDj_x.RDS")




# models of independent evolution of mycorrhizal strategy ----

Myc_dat <- dat[, c("Plant_species", "Mycorrhizal_strategy")]

# with ER model structure, j hidden rates categories, and replicated start x

indepMycERj_x <- corHMM(phy, Myc_dat, rate.cat = j, model = ER) # j varies from 1 to 3 and x from 1 to 5

saveRDS(indepMycERj_x, "indepMycERj_x.RDS")

# with SYM model structure, j hidden rates categories, and replicated start x

indepMycSYMj_x <- corHMM(phy, Myc_dat, rate.cat = j, model = SYM) # j varies from 1 to 3 and x from 1 to 5

saveRDS(indepMycSYMj_x, "indepMycSYMj_x.RDS")

# with ARD model structure, j hidden rates categories, and replicated start x

indepMycARDj_x <- corHMM(phy, Myc_dat, rate.cat = j) # j varies from 1 to 3 and x from 1 to 5

saveRDS(indepMycARDj_x, "indepMycARDj_x.RDS")




# models of dependent evolution ----

# with ER model structure, j hidden rates categories, and replicated start x

depERj_x <- corHMM(phy, dat, rate.cat = j, model = ER) # j varies from 1 to 3 and x from 1 to 5

saveRDS(depERj_x, "depERj_x.RDS")

# with SYM model structure, j hidden rates categories, and replicated start x

depSYMj_x <- corHMM(phy, dat, rate.cat=j, model = SYM) # j varies from 1 to 3 and x from 1 to 5

saveRDS(depSYMj_x, "depSYMj_x.RDS")

# with ARD model structure, j hidden rates categories, and replicated start x

depARDj_x <- corHMM(phy, dat, rate.cat = j) # j varies from 1 to 3 and x from 1 to 5

saveRDS(depARDj_x, "depARDj_x.RDS")


# code to run the ComputeCI function and obtain confidence intervals ----

require(corHMM)

#Dataset_v1

depARD2_3 <- readRDS("depARD2_3.RDS")

CI_dataset_v1 <- ComputeCI(depARD2_3, desired.delta = 2, 5000)

saveRDS(CI_dataset_v1, "CI_dataset_v1.rds")

write.csv(CI_dataset_v1$all_ranges, "CI_dataset_v1.csv")

#Dataset_v2

depARD2_2 <- readRDS("depARD2_2.RDS")

CI_dataset_v2 <- ComputeCI(depARD2_2, desired.delta = 2, 5000)

saveRDS(CI_dataset_v2, "CI_dataset_v2.rds")

write.csv(CI_dataset_v2$all_ranges, "CI_dataset_v2.csv")

#Dataset_v3

depARD2_2 <- readRDS("depARD2_2.RDS")

CI_dataset_v3 <- ComputeCI(depARD2_2, desired.delta = 2, 5000)

saveRDS(CI_dataset_v3, "CI_dataset_v3.rds")

write.csv(CI_dataset_v3$all_ranges, "CI_dataset_v3.csv")

#Dataset_v4

depARD2_3 <- readRDS("depARD2_3.RDS")

CI_dataset_v4 <- ComputeCI(depARD2_3, desired.delta = 2, 5000)

saveRDS(CI_dataset_v4, "CI_dataset_v4.rds")

write.csv(CI_dataset_v4$all_ranges, "CI_dataset_v4.csv")

#Dataset_v5

depSYM2_3 <- readRDS("depSYM2_3.RDS")

CI_dataset_v5 <- ComputeCI(depSYM2_3, desired.delta = 2, 5000)

saveRDS(CI_dataset_v5, "CI_dataset_v5.rds")

write.csv(CI_dataset_v5$all_ranges, "CI_dataset_v5.csv")

#Dataset_v6

depSYM2_1 <- readRDS("depSYM2_1.RDS")

CI_dataset_v6 <- ComputeCI(depSYM2_1, desired.delta = 2, 5000)

saveRDS(CI_dataset_v6, "CI_dataset_v6.rds")

write.csv(CI_dataset_v6$all_ranges, "CI_dataset_v6.csv")

