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



