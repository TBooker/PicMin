rm(list = ls())
library(poolr) # A dependency of PicMin
library(ggplot2) # For making pretty plots
library(cowplot) # For arranging those plots
library(devtools) # Needed to install the PicMin package from Tom's GitHub
#install_github("TBooker/PicMin")
library(tidyverse)
library(PicMin)



# Sam's sugestion...
orderStatsPValues_fix <- function(p_list){
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



FixMin <- function(pList, correlationMatrix, numReps = 100000){
  # Calculate the p-value for the order statistics
  ord_stats_p_values <- orderStatsPValues_fix(pList)
  # Apply the Tippett/Dunn-Sidak Correction
  p_value <- tippett(ord_stats_p_values, adjust = "empirical",
                     R = correlationMatrix,
                     side = 1,
                     size = numReps)$p
  return(list(p=p_value,
              config_est=which.min(ord_stats_p_values)+1))
}


### Perform PicMin on WZA simulation results

n_lineages = 7
# Read WZA_results:

bc_0.003 <- read.csv("~/UBC/GEA/pMax/WZA_analysis/s0.003/BC_Map.wza.csv")
cl_0.003 <- read.csv("~/UBC/GEA/pMax/WZA_analysis/s0.003/cline.wza.csv")
cl_0.003$rep <- cl_0.003$rep + 30
tr_0.003 <- read.csv("~/UBC/GEA/pMax/WZA_analysis/s0.003/trunc.wza.csv")
tr_0.003$rep <- tr_0.003$rep + 60

bc_0.0136 <- read.csv("~/UBC/GEA/pMax/WZA_analysis/s0.0136/BC_Map.wza.csv")
bc_0.0136$rep <- bc_0.0136$rep + 90
cl_0.0136 <- read.csv("~/UBC/GEA/pMax/WZA_analysis/s0.0136/cline.wza.csv")
cl_0.0136$rep <- cl_0.0136$rep + 120
tr_0.0136 <- read.csv("~/UBC/GEA/pMax/WZA_analysis/s0.0136/trunc.wza.csv")
tr_0.0136$rep <- tr_0.0136$rep + 150

neutral <- read.csv("~/UBC/GEA/pMax/WZA_analysis/2D_SteppingStone.wza.csv")
neutral$rep <- neutral$rep + 180


wza_sims_minNeu <- rbind( bc_0.003, 
                   cl_0.003, 
                   tr_0.003,
                   bc_0.0136,
                   cl_0.0136,
                   tr_0.0136)
wza_sims <- rbind(neutral, 
                  wza_sims_minNeu[,names(wza_sims_minNeu)%in%names(neutral)])

reps <- sample( unique( wza_sims$rep ) )

reps1 <- sample( unique( wza_sims$rep ) )
reps2 <- sample( unique( wza_sims$rep ) )
reps3 <- sample( unique( wza_sims$rep ) )
reps4 <- sample( unique( wza_sims$rep ) )
reps5 <- sample( unique( wza_sims$rep ) )
reps6 <- sample( unique( wza_sims$rep ) )

all_reps <- c(reps, reps1, reps2, reps3, reps4, reps5, reps6)

#lineage_sets <- split(reps, ceiling(seq_along(reps)/n_lineages))
lineage_sets <- split(all_reps, ceiling(seq_along(all_reps)/n_lineages))

#lineage_set_i <- lineage_sets[[1]]

make_mini_wza_df <- function(wza_df_for_a_rep, lineage_name){
  tmp <- data.frame( emp_p = wza_df_for_a_rep$empirical_p, 
                     orthogroup = wza_df_for_a_rep$gene)
  names(tmp) <- c(lineage_name, 
                  "orthogroup")
  return( tmp )
}

make_mini_LA_df <- function(wza_df_for_a_rep, lineage_name){
  tmp <- data.frame( LA = wza_df_for_a_rep$LA, 
                     orthogroup = wza_df_for_a_rep$gene)
  names(tmp) <- c(lineage_name, 
                  "orthogroup")
  return( tmp )
}


results_list = list()



for (  r in 1:length(lineage_sets)){

     lineage_set_r <- lineage_sets[[r]]

    
    wza_data_frames <- list()
    LA_data_frames <- list()
    for (p in 1:n_lineages){
      lil_wza <- wza_sims[ wza_sims$rep == lineage_set_r[p], ]
      lil_wza$empirical_p <- PicMin:::EmpiricalPs(vector_of_values = lil_wza$empR_Z, large_i_small_p=TRUE)
      wza_data_frames[[p]] = make_mini_wza_df(lil_wza, paste("species_",p, sep = ""))
      LA_data_frames[[p]] = make_mini_LA_df(lil_wza, paste("species_",p, sep = ""))
      }
    
    # Merge all data frames in list - WZA
    all_lins_wza <- wza_data_frames %>% reduce(full_join, by='orthogroup')
    all_lins_p <- all_lins_wza[ , !(names(all_lins_wza) %in% c("orthogroup"))]
    rownames(all_lins_p) <- all_lins_wza$orthogroup
    
    # Merge all data frames in list - LA
    all_lins_LA <- LA_data_frames %>% reduce(full_join, by='orthogroup')
    all_lins_LA_dropped <- all_lins_LA[ , !(names(all_lins_LA) %in% c("orthogroup"))]
    rownames(all_lins_LA_dropped) <- all_lins_LA$orthogroup
    
    
    # Choose the adaptation threshold - with 4 lineages alpha_a = 0.05 is suitable
    alpha_a <- 0.999
    
    # How many lineages are you analyzing?
    nLins <- 7
    
    # How many samples do you want to take to build the empirical distribution of p-values in poolr?
    # This number needs to be quite large to get accurate p-values. I use 100 for testing and demonstration, but when analyzing the real data it needs to be LARGE - at least 1e6, if not more. When it is large, it takes a long time to run the analysis, but it's important that it's as big as can be. 
    num_reps=10000
    
    # Initialize an output variable 
    count = 0
    
    # Create a container list to store the results in 
    results = list()
    
    count = count + 1
      
      # Run 10,000 replicate simulations of this situation and build the correlation matrix
    emp_p_null_dat <- t(replicate(40000, PicMin:::GenerateNullData(alpha_a, 
                                                                     0.5, 
                                                                     3, 
                                                                    nLins, 
                                                                     10000)))
      
      # Calculate the order statistics p-values for each simulation
    #    emp_p_null_dat_unscaled <- t(apply(emp_p_null_dat ,1, PicMin:::orderStatsPValues))
    emp_p_null_dat_unscaled <- t(apply(emp_p_null_dat ,1, orderStatsPValues_fix))
         
      # Use those p-values to construct the correlation matrix
    null_pMax_cor_unscaled <- cor( emp_p_null_dat_unscaled )
      
      
      # Screen out gene with no evidence for adaptation
    lins_p_screened <- all_lins_p[ apply(all_lins_p<alpha_a,1,function(x) sum(na.omit(x)))!=0, ]
      # Grab out the genes that have the speficied (i.e. the loop variable) amount of missing data
    lins_p_n_screened <-  as.matrix(lins_p_screened[rowSums(is.na(lins_p_screened)) == nLins-nLins,])
      
      # If there is no data in this category, move to the next loop variable
    if (dim(lins_p_n_screened)[1] ==0){
        next
      }
      # Initialise a container for the p-values
    res_p <- rep(-1, 
                   nrow(lins_p_n_screened))
      # Initialise a container for the n_est values
    res_n <- rep(-1, 
                   nrow(lins_p_n_screened))
      
      # Loop over the items in the dataframe (sorry that I'm looping - I'm as disgusted as you!)
    for (i in seq(nrow(lins_p_n_screened)) ){
      print(i)
      test_result <- FixMin(na.omit(lins_p_n_screened[i,]), 
                                       null_pMax_cor_unscaled, 
                                       numReps = num_reps)
      res_p[i] <- test_result$p
      res_n[i] <- test_result$config_est
    }
      
      # Make a dataframe with the results  
    results = data.frame(numLin = 7,
                                    p = res_p,
                                    q = p.adjust(res_p, method = "fdr"),
                                    n_est = res_n,
                                    locus = row.names(lins_p_n_screened) )
      
    my_LA <- all_lins_LA_dropped[row.names(all_lins_LA_dropped)%in%results$locus,]
    numLA <- apply( my_LA, 1, function(x) sum(x>0.005)) 
    numLinked <- apply( my_LA, 1, function(x) sum(x==-99)) 
    results$numLA <- numLA
    results$numLinked <- numLinked
    results$index <- 1:nrow(results)
    results$replicate <- r
    results_list[[r]] = results
}


ggplot(data = results,
       aes(x = index,
           y = -log10(q),
           col = as.factor(numLA),
           shape = as.factor(numLinked)))+
  geom_point()+
  geom_hline(aes(yintercept = -log10(0.05)))

all_results <- do.call(rbind, results_list )


#
write.csv(all_results, "~/UBC/GEA/pMax/WZA_analysis/all_results_r10000.7sets.FixMin099.csv",
          row.names = F)


ggplot(data = all_results,
       aes(x = index,
           y = -log10(q),
           col = as.factor(numLA),
           shape = as.factor(numLinked)))+
  geom_point()+
  facet_wrap(~replicate)+
  geom_hline(aes(yintercept = -log10(0.05)))


ggplot(data = all_results,
       aes(x = index,
           y = -log10(q),
           col = as.factor(n_est == numLA),
           shape = as.factor(numLinked)))+
  geom_point()+
  facet_wrap(~replicate)+
  geom_hline(aes(yintercept = -log10(0.05)))
