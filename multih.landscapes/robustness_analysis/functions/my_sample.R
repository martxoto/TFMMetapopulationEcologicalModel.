# my_sample() is a modification of R function sample()
# it implement the case of sampling one species in a 1-element vector

my_sample <- function(set, n, repl, pi = NULL){
  # set is a vector of elements to sample
  # n is the number of elements to sample in set
  # repl is a logical, =TRUE if same element of set can be sampled several times
  # pi is a vector of probability to pick a given element from set
  if ((n > length(set)) & (repl == FALSE)){
    stop("Cannot sample more elements than present in set if repl = FALSE")
  }
  if (n<0){
    stop("n should be a postive integer.")
  }
  if ((length(set) != length(pi)) & (is.null(pi) == FALSE)){
    stop("pi and set should have the same length.")
  }
  res_sample <- numeric(0)
  
  if (length(set) != 1){
    res_sample <- sample(set, n, replace = repl, prob = pi)
  }
  else {
    res_sample <- rep(set, n)
  }
  return(res_sample)
}