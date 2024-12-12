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

# Compute phylogenetic correlation matrix
phylo_matrix <- vcv(phylo_tree, cor = T)  # The vcv function returns a variance-covariance matrix

# Add an additional column for phylogeny
data <- mutate(data, phylogeny = gsub(" ", "_", species),
               obs = row_number())

# Convert tibble to data frame
data <- as.data.frame(data)

# Calculate VCV matrix
VCV_lnCVR <- vcalc(vi = var_lnCVR,
                   cluster = shared_control_ID, 
                   rho = 0.5,
                   obs = obs,
                   data = data) 

# Model specification
formula <- bf(lnCVR ~ trait_type -1 + # Separate effects by trait type
                (trait_type-1|ref) + # Correlation between traits among studies
                (1|species) + # Species-level random effect
                (1|gr(phylogeny, cov = phylo_matrix)) + # Phylogenetic relatedness
                (1|experiment_ID) + # Experiment-level random effect
                (1|obs) + # Observation-level random effect
                fcor(VCV_lnCVR)) # Variance covariance matix of correlated sampling variances

# Define priors
prior = c(
  prior(constant(1), class = "sigma")
) # Because the residual variance-covariance structure is specified in fcor, there is no need to estimate sigma so it is left as a constant

start <- Sys.time()

# Fit model
lnCVR_model_overall <- brm(formula, 
                             family = gaussian(),
                             data = data, 
                             data2 = list(phylo_matrix = phylo_matrix,
                                          VCV_lnCVR = VCV_lnCVR),
                             prior = prior,
                             control = list(adapt_delta = 0.99, max_treedepth = 15),
                             iter = 4000, 
                             warmup = 2000,
                             chains = 4, 
                             cores = 4, 
                             seed = 123) 

end <- Sys.time()
end-start

saveRDS(lnCVR_model_overall, file = "RData/lnCVR_model_overall.rds")
