`shared_habitat_variation.R` gathers two analyses:
* It examines variance partitioning to estimate how much variation in the plant communities among triads could be accounted for by study site versus habitat type within triads.
* It test for differences in beta-diversity between monads, dyads and triads.

This script requires two input files describing floral abundance in each site and habitat-components in different formats:
* `floral_units_per_site_hab.csv`: a dataframe which lines provide the number of floral units for each plant species sampled across all sites (each column is named with the scientific name of a plant species) in each habitat-component (column `Habitat`) of a given site (column `Site`). Each line is encoded with the corresponding site and habitat (column `Code`), and the type of landscape is given by column `Type`.
* `floral_units_per_site.csv`: a dataframe listing plant species (column `Plant_Species`) and the number of open floral units (column `FU_Open`) in each site. The landscape type (column `Type`) indicates whether the site is a monad (`mon`), a dyad (`di`) or a tryad (`tri`). Hence, each line corresponds to the number of open floral units for a given plant species in a given site.
