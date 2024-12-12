pacman::p_load(tidyverse,
               purrr,
               ape,
               rncl,
               rotl,
               brms,
               metafor,
               emmeans)
set.seed(234)

data <- read_csv(file = "data/processed_data.csv")

unique(data$species) # 25 species

# Match species to the OTL database
taxa <- unique(data$species)
resolved_names <- tnrs_match_names(taxa, context_name = "Animals")
resolved_names # All good, no approximate matches

# Plot tree
tree <- tol_induced_subtree(ott_ids = resolved_names$ott_id)
tree$tip.label <- strip_ott_ids(tree$tip.label) # Remove ott IDs for presentation
plot(tree, no.margin=TRUE, cex=1)

tree <- multi2di(tree, random = TRUE) # Resolve polytomy at random, but it matches classification from Cornetti et al. 2019. Molecular Phylogenetics and Evolution.
phylo_tree <- compute.brlen(tree, method = "Grafen", power = 1) # Compute branch lengths

plot(phylo_tree)

# Add an additional column for phylogeny
data <- mutate(data, phylogeny = gsub(" ", "_", species),
               obs = row_number())

# Filter to only warm temperatures
data_warm <- filter(data, warm_cold == "warm")

# Drop tips that are not in the warm dataset
tips_to_drop <- setdiff(phylo_tree$tip.label, data_warm$phylogeny)
phylo_tree_warm <- drop.tip(phylo_tree, tips_to_drop)

# Compute phylogenetic correlation matrix
phylo_matrix_warm <- vcv(phylo_tree_warm, cor = T)  # The vcv function returns a variance-covariance matrix

# Convert tibble to data frame
data_warm <- as.data.frame(data_warm)

# Calculate VCV matrix
VCV_lnRR_warm <- vcalc(vi = var_lnRR,
                        cluster = shared_control_ID, 
                        rho = 0.5,
                        obs = obs,
                        data = data_warm) 

# Model specification
formula <- bf(lnRR ~ 0 + trait_type:constant_increasing + trait_type:constant_increasing:assay_temp_diff + # Separate effects by trait type
                (trait_type-1|ref) + # Correlation between traits among studies
                (1|species) + # Species-level random effect
                (1|gr(phylogeny, cov = phylo_matrix_warm)) + # Phylogenetic relatedness
                (1|experiment_ID) + # Experiment-level random effect
                (1|obs) + # Observation-level random effect
                fcor(VCV_lnRR_warm)) # Variance covariance matix of correlated sampling variances

# Define priors
prior = c(
  prior(constant(1), class = "sigma")
) # Because the residual variance-covariance structure is specified in fcor, there is no need to estimate sigma so it is left as a constant

start <- Sys.time()

# Fit model
lnRR_model_constant_increasing <- brm(formula, 
                                       family = gaussian(),
                                       data = data_warm, 
                                       data2 = list(phylo_matrix_warm = phylo_matrix_warm,
                                                    VCV_lnRR_warm = VCV_lnRR_warm),
                                       prior = prior,
                                       control = list(adapt_delta = 0.99, max_treedepth = 15),
                                       iter = 4000, 
                                       warmup = 2000,
                                       chains = 4, 
                                       cores = 4, 
                                       seed = 123) 

end <- Sys.time()
end-start

saveRDS(lnRR_model_constant_increasing, file = "RData/lnRR_model_constant_increasing.rds")