#HYPERGEOMETRIC TEST FOR TABLE 2 

  #LEVEL 0 SOMATIC LC DRIVER GENE ENRICHMENT (LUNG GRN) 
  #Define the parameters of the hypergeometric distribution
n_sample_successes <- 1    # number of successes in the sample
n_sample <- 64          # sample size
n_pop_success <- 43       # number of successes in the population
n_pop <- 15855            # population size

# Perform the one-sided hypergeometric test
pvalue <- phyper(n_sample_successes - 1, n_pop_success, n_pop - n_pop_success, n_sample, lower.tail = FALSE)

# Print the p-value
print(pvalue)

  #LEVEL 1 SOMATIC LC DRIVER GENE ENRICHMENT (LUNG PPIN) - STRING 
  #Define the parameters of the hypergeometric distribution
n_sample_successes <- 4    # number of successes in the sample
n_sample <- 140          # sample size
n_pop_success <- 43       # number of successes in the population
n_pop <- 11507            # population size

# Perform the one-sided hypergeometric test
pvalue <- phyper(n_sample_successes - 1, n_pop_success, n_pop - n_pop_success, n_sample, lower.tail = FALSE)

# Print the p-value
print(pvalue)

  #LEVEL 2 SOMATIC LC DRIVER GENE ENRICHMENT (LUNG PPIN) - STRING 
  #Define the parameters of the hypergeometric distribution
n_sample_successes <- 6    # number of successes in the sample
n_sample <- 347          # sample size
n_pop_success <- 43       # number of successes in the population
n_pop <- 11507            # population size

# Perform the one-sided hypergeometric test
pvalue <- phyper(n_sample_successes - 1, n_pop_success, n_pop - n_pop_success, n_sample, lower.tail = FALSE)

# Print the p-value
print(pvalue)

  #LEVEL 1 SOMATIC LC DRIVER GENE ENRICHMENT (REST-OF-BODY PPIN) - STRING
  #Define the parameters of the hypergeometric distribution
n_sample_successes <- 2    # number of successes in the sample
n_sample <- 66          # sample size
n_pop_success <- 19       # number of successes in the population
n_pop <- 8056            # population size

# Perform the one-sided hypergeometric test
pvalue <- phyper(n_sample_successes - 1, n_pop_success, n_pop - n_pop_success, n_sample, lower.tail = FALSE)

# Print the p-value
print(pvalue)

  #LEVEL 2 SOMATIC LC DRIVER GENE ENRICHMENT (REST-OF-BODY PPIN) - STRING
  #Define the parameters of the hypergeometric distribution
n_sample_successes <- 2    # number of successes in the sample
n_sample <- 189          # sample size
n_pop_success <- 19       # number of successes in the population
n_pop <- 8056            # population size

# Perform the one-sided hypergeometric test
pvalue <- phyper(n_sample_successes - 1, n_pop_success, n_pop - n_pop_success, n_sample, lower.tail = FALSE)

# Print the p-value
print(pvalue)

#HYPERGEOMETRIC TEST FOR SUPPLEMENTARY TABLE 2 

  #LEVEL 1 SOMATIC LC DRIVER GENE ENRICHMENT (LUNG PPIN) - STRING and PROPER
  #Define the parameters of the hypergeometric distribution
n_sample_successes <- 5    # number of successes in the sample
n_sample <- 380          # sample size
n_pop_success <- 43       # number of successes in the population
n_pop <- 11648            # population size

# Perform the one-sided hypergeometric test
pvalue <- phyper(n_sample_successes - 1, n_pop_success, n_pop - n_pop_success, n_sample, lower.tail = FALSE)

# Print the p-value
print(pvalue)

  #LEVEL 2 SOMATIC LC DRIVER GENE ENRICHMENT (LUNG PPIN) - STRING and PROPER
  #Define the parameters of the hypergeometric distribution
n_sample_successes <- 27    # number of successes in the sample
n_sample <- 5335          # sample size
n_pop_success <- 43       # number of successes in the population
n_pop <- 11648            # population size

# Perform the one-sided hypergeometric test
pvalue <- phyper(n_sample_successes - 1, n_pop_success, n_pop - n_pop_success, n_sample, lower.tail = FALSE)

# Print the p-value
print(pvalue)

  #LEVEL 1 SOMATIC LC DRIVER GENE ENRICHMENT (LUNG PPIN) - PROPER
#Define the parameters of the hypergeometric distribution
n_sample_successes <- 2    # number of successes in the sample
n_sample <- 248          # sample size
n_pop_success <- 40       # number of successes in the population
n_pop <- 5947            # population size

# Perform the one-sided hypergeometric test
pvalue <- phyper(n_sample_successes - 1, n_pop_success, n_pop - n_pop_success, n_sample, lower.tail = FALSE)

# Print the p-value
print(pvalue)

  #LEVEL 2 SOMATIC LC DRIVER GENE ENRICHMENT (LUNG PPIN) - PROPER
#Define the parameters of the hypergeometric distribution
n_sample_successes <- 28    # number of successes in the sample
n_sample <- 5120          # sample size
n_pop_success <- 40       # number of successes in the population
n_pop <- 5947            # population size

# Perform the one-sided hypergeometric test
pvalue <- phyper(n_sample_successes - 1, n_pop_success, n_pop - n_pop_success, n_sample, lower.tail = FALSE)

# Print the p-value
print(pvalue)

