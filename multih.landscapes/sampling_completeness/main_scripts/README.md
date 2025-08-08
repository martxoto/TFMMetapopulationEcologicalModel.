# Exploring sampling completeness in sampled networks

The following scripts explore the level of sampling completeness of pollination interactions within triads and monads:

* `habitat_sc.R`: This script draws rarefaction curves/sampling coverage curves, and explores how we can subsample dataset to the same level of sampling completeness.

* `habitat_sc_same_plcomsize.R`: This script does a similar analysis to the above, and specifically checks whether preserving the size of the plant community within the null monads changes their sampling completeness and the number of interaction events to subsample for the second null model.

The following cover all types of interactions covered by the sampling protocole :

* `site_sc.R` : This script calculates the interaction sampling coverage within each site. All analyses rely on transposing the concepts of species richness estimation to interaction richness.

* `int-evenness-vs-samp-cover.R` : This script explores the relationship between interaction eveneness and interaction sampling coverage.

All of these analyses rely on transposing the concepts of species richness estimation to interaction richness. These scripts require the R package `iNEXT` by [Hsieh et al. (2016)](https://besjournals.onlinelibrary.wiley.com/doi/pdf/10.1111/2041-210X.12613) to run.
