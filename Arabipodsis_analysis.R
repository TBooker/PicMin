rm(list = ls())

## You'll need to install poolr to use the PicMin function
library(poolr)
library(ggplot2)
library(cowplot)

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
PicMin <- function(pList, correlationMatrix, numReps = 100000){
  # Calculate the p-value for the order statistics
  ord_stats_p_values <- orderStatsPValues(pList)
  # Apply the Tippett/Dunn-Sidak Correction
  p_value <- tippett(ord_stats_p_values, adjust = "empirical", 
                     R = correlationMatrix, 
                     side = 1, 
#                     size = numReps,
                     size = c(1000, 10000, 1000000), threshold = c(.10, .01))$p
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

## Perform Picmin on the Arabidopsis data from Bohutinska et al 

## quick function for calculating empirical p-values

empirical_ps <- function( vector_of_values ){
  1-rank(vector_of_values)/length(vector_of_values)
}


## quick function for getting names

get_names <- function( lineage_df ){
  return( paste( lineage_df$scaff, lineage_df$start,
         sep = "_") )
}


## quick function for minimising the dataframes

min_lin <- function( lineage_df , lineage_name){
  tmp <- data.frame( emp_p = lineage_df$emp_p, 
              window = lineage_df$name)
  names(tmp) <- c(lineage_name, 
                  "window")
  return( tmp )
}


## Read in the data, calculate empirical p-values and use a common naming scheme for all lineages

lin_1 <- read.csv("~/UBC/GEA/pMax/Arabidopsis/BPM/BALTIS_WS1000_MS1_BPM.txt",
                  sep = "\t")
lin_1$emp_p <- empirical_ps(lin_1$FstH)
lin_1$name <- get_names( lin_1 )
lin_1_m <- min_lin( lin_1, "BALTIS" )

lin_2 <- read.csv("~/UBC/GEA/pMax/Arabidopsis/BPM/HCADRG_WS1000_MS1_BPM.txt",
                  sep = "\t")
lin_2$emp_p <- empirical_ps(lin_2$FstH)
lin_2$name <- get_names( lin_2 )
lin_2_m <- min_lin( lin_2, "HCADRG" )

lin_3 <- read.csv("~/UBC/GEA/pMax/Arabidopsis/BPM/INECAR_WS1000_MS1_BPM.txt",
                  sep = "\t")
lin_3$emp_p <- empirical_ps(lin_3$FstH)
lin_3$name <- get_names( lin_3 )
lin_3_m <- min_lin( lin_3, "INECAR" )

lin_4 <- read.csv("~/UBC/GEA/pMax/Arabidopsis/BPM/OBIGUN_WS1000_MS1_BPM.txt",
                  sep = "\t")
lin_4$emp_p <- empirical_ps(lin_4$FstH)
lin_4$name <- get_names( lin_4 )
lin_4_m <- min_lin( lin_4, "OBIGUN" )

lin_5 <- read.csv("~/UBC/GEA/pMax/Arabidopsis/BPM/TKOHRA_WS1000_MS1_BPM.txt",
                  sep = "\t")
lin_5$emp_p <- empirical_ps(lin_5$FstH)
lin_5$name <- get_names( lin_5 )
lin_5_m <- min_lin( lin_5, "TKOHRA" )

lin_6 <- read.csv("~/UBC/GEA/pMax/Arabidopsis/BPM/WILKAS_WS1000_MS1_BPM.txt",
                  sep = "\t")
lin_6$emp_p <- empirical_ps(lin_6$FstH)
lin_6$name <- get_names( lin_6 )
lin_6_m <- min_lin( lin_6, "WILKAS" )

lin_7 <- read.csv("~/UBC/GEA/pMax/Arabidopsis/BPM/ZEPSUB_WS1000_MS1_BPM.txt",
                  sep = "\t")
lin_7$emp_p <- empirical_ps(lin_7$FstH)
lin_7$name <- get_names( lin_7 )
lin_7_m <- min_lin( lin_7, "ZEPSUB" )

## Merge dataframes

library(tidyverse)

#put all data frames into list
df_list <- list(lin_1_m,
                lin_2_m,
                lin_3_m,
                lin_4_m,
                lin_5_m,
                lin_6_m,
                lin_7_m)

#merge all data frames in list
all_lins <- df_list %>% reduce(full_join, by='window')
all_lins_p <- all_lins[ , !(names(all_lins) %in% c("window"))]
rownames(all_lins_p) <- all_lins$window
c(all_lins_p[rownames(all_lins_p)=="scaffold_7_8621000",])


alpha_a <- 0.05
nLins <- 7

count = 0
results = list()
for (n in c(4,5,6,7)){
  count = count + 1
  # Run 10,000 replicate simulations of this situation and build the correlation matrix
  emp_p_null_dat <- t(replicate(40000, GenerateNullData(alpha_a, 0.5, 3, n, 10000)))
  # Calculate the order statistics p-values for each simulation
  emp_p_null_dat_unscaled <- t(apply(emp_p_null_dat ,1, orderStatsPValues))
  # Use those p-values to construct the correlation matrix
  null_pMax_cor_unscaled <- cor( emp_p_null_dat_unscaled )

  
  # Screen out gene with no evidence for adaptation
  lins_p_screened <- all_lins_p[ apply(all_lins_p<alpha_a,1,function(x) sum(na.omit(x)))!=0, ]
  lins_p_n_screened <-  as.matrix(lins_p_screened[rowSums(is.na(lins_p_screened)) == nLins-n,])

  if (dim(lins_p_n_screened)[1] ==0){
    next
  }
  res_p <- rep(-1, 
               nrow(lins_p_n_screened))
  res_n <- rep(-1, 
               nrow(lins_p_n_screened))
  
  for (i in seq(nrow(lins_p_n_screened)) ){
    test_result <- PicMin(na.omit(lins_p_n_screened[i,]), null_pMax_cor_unscaled, numReps = 10000)
    res_p[i] <- test_result$p
    res_n[i] <- test_result$config_est
  }
  results[[count]] = data.frame(numLin = n ,
                                p = res_p,
                                 q = p.adjust(res_p, method = "fdr"),
                                 n_est = res_n,
                                 locus = row.names(lins_p_n_screened) )

  }


picMin_results <- do.call(rbind, results) 

head(picMin_results)

picMin_results$pooled_q <- p.adjust(picMin_results$p, method = "fdr")


picMin_results <- cbind( picMin_results, 
       read.csv(text=picMin_results$locus, header=FALSE,
         sep = "_",
         col.names=c('redundan','scaffold','start'))
)

library(ggplot2)

str(picMin_results)

col_pal <- c("white", "#8ec641", "#897696", "#e93826", "#13a4f5", "#f89b56")

picMin_results$scaffold <- factor(picMin_results$scaffold,
                                  labels = paste("Scaffold",1:8))
  

SuppFigure4 <- ggplot(data = picMin_results, 
       aes(x = start/1e6,
           y = -log10(pooled_q),
           fill = factor(n_est)))+
  geom_point(shape = 21,
             size = 4)+
  geom_hline(aes(yintercept = -log10(0.05)),
             lty=2)+
  facet_wrap(~scaffold,
             ncol = 4,
             scales = "free_x")+
  scale_fill_manual(expression(italic(n[est])),values = col_pal)+
  scale_y_continuous(expression(-log[10]*"("*italic("q")*"-value)"))+
  scale_x_continuous("Position in Scaffold (Mbp)")+
  theme_half_open() +
  theme(strip.background = element_blank())+
  background_grid()# always place this after the theme
  
#write.csv(picMin_results,
#          file = "~/UBC/GEA/pMax/Arabidopsis/PicMinResults.csv", row.names = F)

ggsave(SuppFigure4,
       file = "~/UBC/GEA/pMax/Arabidopsis/arabidopsisPicMin.png",
       width = 9.50,
       height = 5.00)


picMin_results<- read.csv("~/UBC/GEA/pMax/Arabidopsis/PicMinResults.csv")

picMin_results[picMin_results$pooled_q < 0.2,]

hits <- picMin_results[picMin_results$p<0.01,] 
dim(hits)
hits[hits$n_est==2,]
#write.csv(hits,
#          file = "~/UBC/GEA/pMax/Arabidopsis/arabidopsisHits.csv", row.names = F)


### Subset


#put all data frames into list
df_list <- list(lin_1_m,
#                lin_2_m,
                lin_3_m,
#                lin_4_m,
                lin_5_m,
                lin_6_m,
                lin_7_m)

#merge all data frames in list
all_lins <- df_list %>% reduce(full_join, by='window')
all_lins_p <- all_lins[ , !(names(all_lins) %in% c("window"))]
rownames(all_lins_p) <- all_lins$window


alpha_a <- 0.05
nLins <- 5

count = 0
results = list()
for (n in c(3,4,5)){
  count = count + 1
  # Run 10,000 replicate simulations of this situation and build the correlation matrix
  emp_p_null_dat <- t(replicate(40000, GenerateNullData(alpha_a, 0.5, 3, n, 10000)))
  # Calculate the order statistics p-values for each simulation
  emp_p_null_dat_unscaled <- t(apply(emp_p_null_dat ,1, orderStatsPValues))
  # Use those p-values to construct the correlation matrix
  null_pMax_cor_unscaled <- cor( emp_p_null_dat_unscaled )
  
  
  # Screen out gene with no evidence for adaptation
  lins_p_screened <- all_lins_p[ apply(all_lins_p<alpha_a,1,function(x) sum(na.omit(x)))!=0, ]
  lins_p_n_screened <-  as.matrix(lins_p_screened[rowSums(is.na(lins_p_screened)) == nLins-n,])
  
  if (dim(lins_p_n_screened)[1] ==0){
    next
  }
  res_p <- rep(-1, 
               nrow(lins_p_n_screened))
  res_n <- rep(-1, 
               nrow(lins_p_n_screened))
  
  for (i in seq(nrow(lins_p_n_screened)) ){
    test_result <- PicMin(na.omit(lins_p_n_screened[i,]), null_pMax_cor_unscaled, numReps = 10000)
    res_p[i] <- test_result$p
    res_n[i] <- test_result$config_est
  }
  results[[count]] = data.frame(numLin = n ,
                                p = res_p,
                                q = p.adjust(res_p, method = "fdr"),
                                n_est = res_n,
                                locus = row.names(lins_p_n_screened) )
  
}


picMin_results <- do.call(rbind, results) 

head(picMin_results)

picMin_results$pooled_q <- p.adjust(picMin_results$p, method = "fdr")


picMin_results <- cbind( picMin_results, 
                         read.csv(text=picMin_results$locus, header=FALSE,
                                  sep = "_",
                                  col.names=c('redundan','scaffold','start'))
)

library(ggplot2)

str(picMin_results)

col_pal <- c("#fcbe57", "#8ec641", "#897696", "#e93826", "#13a4f5", "#f89b56")

picMin_results$scaffold <- factor(picMin_results$scaffold,
                                  labels = paste("Scaffold",1:8))

picMin_results <- read.csv("~/UBC/GEA/pMax/Arabidopsis/PicMinResults_subset.csv")

ggplot(data = picMin_results, 
       aes(x = start/1e6,
           y = -log10(pooled_q),
           fill = factor(n_est)))+
  geom_point(shape = 21,
             size = 4)+
  geom_hline(aes(yintercept = -log10(0.05)),
             lty=2)+
  facet_wrap(~scaffold,
             ncol = 4,
             scales = "free_x")+
  scale_fill_manual("Number of\nLineages\n",values = col_pal)+
  scale_y_continuous(expression(-log[10]*"(q-value)"))+
  scale_x_continuous("Position in Scaffold (Mbp)")+
  theme_half_open() +
  theme(strip.background = element_blank())+
  background_grid()# always place this after the theme

write.csv(picMin_results,
          file = "~/UBC/GEA/pMax/Arabidopsis/PicMinResults_subset.csv", row.names = F)


  hits <- picMin_results[picMin_results$pooled_q<0.05,] 

write.csv(hits,
          file = "~/UBC/GEA/pMax/Arabidopsis/arabidopsisHits_subset.csv", row.names = F)



all_lins_p$binary_classification <- apply(all_lins_p<0.05, 1, function(x) sum(na.omit(x)) )

hits_all_lins <- all_lins_p[row.names(all_lins_p)%in%hits$locus,]
hits_all_lins$locus <- row.names(hits_all_lins)

merged_hits <- list(hits,hits_all_lins) %>% reduce(full_join, by='locus')

plot( -log10(merged_hits$p), merged_hits$binary_classification/merged_hits$n_est )

