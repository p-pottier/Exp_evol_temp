pacman::p_load(tidyverse,
               purrr,
               ape,
               rncl,
               rotl,
               brms,
               metafor,
               emmeans,
               future.apply)
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

# Get a list of study IDs for the LOO
study_ids <- unique(data$ref)
num_studies <- length(study_ids)

# Initialize a list to store the results
loo_results <- vector("list", num_studies)
names(loo_results) <- study_ids

# Set up parallel processing (16 models at a time)
plan(multisession, workers = 16)  

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


# Function to fit the model after iteratively removing one study
leave_one_out <- function(study_id) {
  # Subset the data
  data_subset <- subset(data, ref != study_id)
  
  # Recalculate VCV matrix
  VCV_lnCVR_subset <- vcalc(
    vi = var_lnCVR,
    cluster = shared_control_ID, 
    rho = 0.5,
    obs = obs,
    data = data_subset
  )
  
  # Fit the model
  model_subset <- brm(
    formula, 
    family = gaussian(),
    data = data_subset, 
    data2 = list(
      phylo_matrix = phylo_matrix,
      VCV_lnCVR = VCV_lnCVR_subset
    ),
    prior = prior,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    iter = 4000, 
    warmup = 2000,
    chains = 2, # Only two chains to reduce computational demands
    cores = 1,   # Set cores to 1 inside the function
    seed = 123
  )
  
  # Extract the fixed effects
  fixed_effects <- fixef(model_subset)
  
  # Return the fixed effects
  return(list(study_id = study_id, fixed_effects = fixed_effects))
}

start <- Sys.time()

# Run the function in parallel over the studies
loo_results <- future_lapply(study_ids, leave_one_out)

end <- Sys.time()
end-start

saveRDS(loo_results, file = "RData/lnCVR_LOO_results.rds")
