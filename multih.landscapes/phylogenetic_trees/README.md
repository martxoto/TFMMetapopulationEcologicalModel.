This folders gathers scripts aiming at generating phylogenetic trees for the plant communities observed across sites:

* `phylogeny_mm.R` creates a phylogenetic tree for the plant species found across sites thanks to Daphne (a dated phylogeny of European flora by Durka & Michalski, 2012), implementing the solutions described in the Methods for observed plant taxa that are missing in Daphne.
* `phylogeny_func.R` creates the phylogenetic tree for a given plant community, and calculates the phylogenetic diversity - measured as the sum of branch lengths. The tree is a chop of the Daphne tree. The scripts also illustrates the use of the function `PD()`. Note that this function is isolated and used in the analysis on potential emergent properties within triads (see `emergent_properties`).

