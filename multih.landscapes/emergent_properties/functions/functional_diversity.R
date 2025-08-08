# Here are gathered a set of functions used for the functional diversity analysis.

required_libraries <- c("bipartite", "ape", "FD", "vegan", "usedist")

for (lib in required_libraries){
  if (!require(lib, character.only = TRUE)){
    install.packages(lib)
  }
  library(lib, character.only = TRUE)
}


# freq2prop() changes interaction frequencies into proportions of interaction events for a given taxon.
freq2prop <- function(x){
  newx <- x / sum(x)
  return(newx)
}

# bespokedist() returns a distance matrix of diets for taxa within the upper guild.
bespokedist <- function(web, datatype, distmeth){
  # verify input variables
  if (!is.matrix(web)){
    stop("web must be a matrix which non-zero entries correspond to interaction between species in rows and columns.")
  }
  if (!any(datatype %in% c("bin", "raw", "sumto1"))){
    stop("datatype is a character string describing the type of entries to web. It is either 
         equal to 'bin' (for binary entries), or 'raw' (for no transformation of web entries),
         or 'sumto1' (to convert interaction frequencies into proportions).")
  }
  if (!any(distmeth %in% c("bray", "jaccard"))){
    stop("distmeth is a character string for the type of dissimilarity index used to calculate the distance matrix.
         It is either 'bray' (for Bray-Curtis) or 'jaccard' (for Jaccard index).")
  }
  
  # transform data if required
  if (datatype == "bin"){
    web[web != 0] <- 1
  }
  if(datatype == "sumto1"){
    web <- apply(web, 2, freq2prop) # for each taxon from the upper guild, turn interaction frequencies into probabilities
  }
  # calculate the distance matrix
  # distmat <- as.matrix(vegdist(t(web), method = distmeth, diag = TRUE, upper = FALSE))
  distmat <- vegdist(t(web), method = distmeth, diag = TRUE, upper = FALSE)
  
  return(distmat)
}

# func_div_multisites() returns various indices of functional diversity across multiple networks
func_div_multsites <- function(df_web, dist_mat, ndim){
  # verify input variables
  if (!is.data.frame(df_web)){
    stop("df_web must be a dataframe describing interactions across sites.")
  }
  else {
    if ((!"Site" %in% colnames(df_web)) | (!"Lower_Taxon" %in% colnames(df_web)) | (!"Upper_Taxon" %in% colnames(df_web))){
      stop("df_web must have at least three columns named 'Site', 'Lower_Taxon', and 'Upper_Taxon'.")
    }
  }
  ug_taxa <- sort(unique(df_web$Upper_Taxon)); N_ug_taxa <- length(ug_taxa) # the set of taxa within the upper guild
  if (!any(is(dist_mat) == "dist")){
    stop("dist_mat must be a matrix describing the dissimilarity of diet among taxa within the upper guild.")
  }
  if ((!is.numeric(ndim)) | (any(ndim %% 1 != 0))){
    stop("ndim must be an integer value.")
  }
  if ((!all(ug_taxa %in% labels(dist_mat))) | (!all(labels(dist_mat) %in% ug_taxa))){
    stop("Taxa within the upper guild in web (column names) must be the same as those named in the distance matrix ('dist_mat').")
  }

  site_list <- sort(unique(df_web$Site)); N_sites <- length(site_list)
  abund_ug <- matrix(0, N_sites, N_ug_taxa, dimnames = list(site_list, ug_taxa))
  for (site in site_list){
    ug_taxa_site <- sort(unique(df_web$Upper_Taxon[df_web$Site == site]))
    abund_ug_site <- table(df_web$Upper_Taxon[df_web$Site == site])
    abund_ug[site, ug_taxa_site] <- abund_ug_site[match(ug_taxa_site, names(abund_ug_site))]
  }
  occ_ug <- abund_ug; occ_ug[occ_ug != 0] <- 1
  
  FD_df <- data.frame(Site = rep(site_list, length(ndim)), NoDim = rep(ndim, N_sites), FRic_occ = NA, FEve_occ = NA, FDiv_occ = NA, FDis_occ = NA, FEve_abund = NA, FDiv_abund = NA, FDis_abund = NA)
  for (nb in ndim){
    FD_occ <- dbFD(dist_mat, occ_ug, corr = "cailliez", stand.FRic = TRUE, m = nb, messages = FALSE)
    FD_abund <- dbFD(dist_mat, abund_ug, corr = "cailliez", stand.FRic = TRUE, m = nb, messages = FALSE)
    FD_df[(FD_df$NoDim == nb), ] <- data.frame(Site = site_list, NoDim = nb, FRic_occ = FD_occ$FRic, FEve_occ = FD_occ$FEve, FDiv_occ = FD_occ$FDiv, FDis_occ = FD_occ$FDis,
                                             FEve_abund = FD_abund$FEve, FDiv_abund = FD_abund$FDiv, FDis_abund = FD_abund$FDis)
  }
  
  return(FD_df)
}