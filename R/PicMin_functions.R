#' @title orderStatsPValues_fix
#' @param p_list the vector of p-values from your genome scans
#' @importFrom "stats" "pbeta"
orderStatsPValues <- function(p_list){
  ## This function returns a list of the p-values for each of marginal p-values
  # sort the list of $p$-values
  p_sort <- sort(p_list)[2:length(p_list)]
  # calculate the number of species minus 1
  n = length(p_sort)
  # get a vector of the 'a' parameters for each of the marginal distributions
  the_as = 2:(length(p_sort) + 1)
  # get a vector of the 'b' parameters for each of the marginal distributions
  the_bs = n+2-the_as
  # calculate the p-values for each of the marginals and return
  return(pbeta(p_sort, the_as, the_bs))
}





#' @title PicMin
#' @param pList the vector of p-values from your genome scans
#' @param correlationMatrix the correlation matrix under the null hypothesis
#' @param numReps the number of replicate draws to perform when building the empirical distributing for calculating the Tippett p-value
#' @importFrom "poolr" "tippett"
PicMin <- function(pList, correlationMatrix, numReps = 100000){
  # Calculate the p-value for the order statistics
  ord_stats_p_values <- orderStatsPValues_fix(pList)
  # Apply the Tippett/Dunn-Sidak Correction
  p_value <- tippett(ord_stats_p_values, adjust = "empirical",
                     R = correlationMatrix,
                     side = 1,
                     size = c(100,1000,10000,100000,1000000),
                     threshold = c(0.1,0.01,0.001,0.0001, 0.00001)
                     )$p
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
GenerateNullData <- function(adaptation_screen, n, a = 0.3, b = 5, genes = 20000, shuffle  = TRUE){
  temp <- c( order( c(rbeta(1,0.3,5), runif(genes-1)))[1]/genes,
             replicate(n-1, sample(genes,1)/genes) )
  while (sum(temp<adaptation_screen)==0){
    temp <- c( order( c(rbeta(1,0.3,5), runif(genes-1)))[1]/genes,
               replicate(n-1, sample(genes,1)/genes) )
  }
  if (shuffle==TRUE){
    return(sample(temp))
  }
  else{
    return(temp)

  }
}


#' @title Calculate empirical p-values from a vector of summary statistics
#' @param vector_of_values A set of summary statistics that you want to convert to empirical p-values
#' @param large_i_small_p Do you want large values to have small p-values (e.g. Fst)?
EmpiricalPs <- function( vector_of_values, large_i_small_p = FALSE ){
  if  (large_i_small_p==TRUE){
    rank(vector_of_values * -1,na.last = "keep")/sum(is.na (vector_of_values) == F)
  }
  else{
    rank(vector_of_values,na.last = "keep")/sum(is.na (vector_of_values) == F)
  }
}


#' @title Extract locus names from a dataframe of the kind used in the vignette  (written for the vignette, so may not be widely applicable)
#' @param lineage_df A dataframe containing the data you want to extract names from. Must have the columns "scaff" and "start".
get_names <- function( lineage_df ){
  return( paste( lineage_df$scaff, lineage_df$start,
                 sep = "_") )
}



#' @title Reduce a dataframe to contain only the empirical p-values for a particular lineage (written for the vignette, so may not be widely applicable)
#' @param lineage_df A dataframe containing the data you want to extract names from.
#' @param lineage_name A unique identifier to use as a column heading
min_lin <- function( lineage_df , lineage_name){
  tmp <- data.frame( emp_p = lineage_df$emp_p,
                     window = lineage_df$name)
  names(tmp) <- c(lineage_name,
                  "window")
  return( tmp )
}
