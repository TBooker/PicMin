#' @title orderStatsPValues
#' @param p_list the vector of p-values from your genome scans
#' @importFrom "stats" "pbeta"
orderStatsPValues <- function(p_list){
  ## This function returns a list of the p-values for each of marginal p-values
  # sort the list of $p$-values
  p_sort <- sort(p_list)[2:length(p_list)]
  # calculate the number of species minus 1
  n = length(p_sort)
  # get a vector of the 'a' parameters for each of the marginal distributions
  the_as = 1:length(p_sort)
  # get a vector of the 'b' parameters for each of the marginal distributions
  the_bs = n+1-the_as
  # calculate the p-values for each of the marginals and return
  return(pbeta(p_sort, the_as, the_bs))
}


#' @title Generate data under the null hypothesis
#' @param pList the vector of p-values from your genome scans
#' @param correlationMatrix the correlation matrix under the null hypothesis
#' @param numReps the number of replicate draws to perform when building the empirical distributing for calculating the Tippett p-value
#' @importFrom "poolr" "tippett"
PicMin <- function(pList, correlationMatrix){
  # Calculate the p-value for the order statistics
  ord_stats_p_values <- orderStatsPValues(pList)
  # Apply the Tippett/Dunn-Sidak Correction
  p_value <- tippett(ord_stats_p_values, adjust = "empirical",
                     R = correlationMatrix,
                     side = 1,
                     size = c(1000, 10000, 100000, 5000000),
                     threshold = c(.10, .01, 0.001))$p
  return(list(p=p_value,
              config_est=which.min(ord_stats_p_values)+1))
}



#' @title Generate data under the null hypothesis
#' @param adaptation_screen The threshold used to determine adaptation
#' @param a the 'a' parameter of a beta distribution of p-values for the false null
#' @param b the 'b' parameter of a beta distribution of p-values for the false null
#' @param n the number of species in the test
#' @param genes the number of genes in the genome use to calculate empirical p-values
#' @importFrom "stats" "rbeta"
GenerateNullData <- function(adaptation_screen, a, b, n, genes){
  temp <- c( rbeta(1,a,b),replicate(n-1, sample(genes,1)/genes) )
  while (sum(temp<adaptation_screen)==0){
    temp <- c( rbeta(1,a,b),replicate(n-1, sample(genes,1)/genes) )
  }
  return(temp)
}



#' @title Calculate empirical p-values from a vector of summary statistics
#' @param vector_of_values A set of summary statistics that you want to convert to empirical p-values
#' @param large_i_small_p Do you want large values to have small p-values (e.g. Fst)?
EmpiricalPs <- function( vector_of_values, large_i_small_p = FALSE ){
  if  (large_i_small_p==FALSE){
    1-(rank(vector_of_values))/length(vector_of_values)
  }
  else{
    rank(vector_of_values)/length(vector_of_values)
  }
}
