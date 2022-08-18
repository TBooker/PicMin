rm(list=ls())

colorBlindBlack8  <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

library(rSymPy)
library(poolr)
library(reshape2)
library(ggplot2)
library(foreach)
library(doParallel)
library(ggpubr)


#' @title Beta probability
#' @param p p-value
#' @param i rank
#' @param n The number of inputs
#' @importFrom "stats" "pbeta" "qbeta"
F_i = function(p, i, n)
{
  a = i
  b = n-i+1
  res = pbeta(q = p, shape1 = a, shape2 = b, lower.tail = T)
  return(res)
}


#' @title ordmeta
#' @description Minimum Marginal P-value in joint order distribution
#' @param p A vector of p-values
#' @param is.onetail Logical. If set TRUE, p-values are combined without considering the direction of effect, and vice versa. Default: TRUE.
#' @param eff.sign A vector of signs of effect sizes. It works when is.onetail = FALSE
#' @importFrom "rSymPy" "Var" "sympy"
#' @return p : Combined p-value
#' @return optimal_rank : Optimal rank where minimum marginal p-value exists.
#' @return eff.p.idx : Index of effective p-values
#' @return MMP : Minimum marginal p-value
#' @return overall.eff.direction : The direction of combined effects.
#' @examples \donttest{ordmeta(p=c(0.01, 0.02, 0.8, 0.25), is.onetail=FALSE, eff.sign = c(1,1,1,-1))}
#' @export

ordmeta = function(p, is.onetail = TRUE, eff.sign=NULL)
{
  direc = eff.sign
  if(is.null(p)){stop("Input p-values are required.")}
  if(!is.onetail & is.null(eff.sign)){stop("Input the direction of effects.")}
  idx_na = which(is.na(p))
  if(length(idx_na)>0){p = p[-idx_na]; eff.sign = eff.sign[-idx_na]}
  ordmeta = function(p2)
  {
    ord = order(p2, decreasing = F)
    pord = sort(p2, decreasing = F)

    # get alpha = MIN(F_(i)(x)) {i={1..n}}
    N = length(p2)
    alpha = 1.01 # an arbitrary number larger than 1
    for(i in 1:N)
    {
      alpha_temp = F_i(pord[i], i, N)
      if(alpha_temp < alpha){idx_minimum = i; alpha = alpha_temp}
    }
    # symbolic integral
    for(i in 1:N)
    {
      x = Var("x")
      y = Var("y")
      if(i==1)
      {
        templete = paste(i,"*integrate(1, (x, lob, y))")
        lob = qbeta(p = alpha,shape1 = i, shape2 = N-i+1, lower.tail = T)
        templete = gsub("lob", lob, templete)
      }else if(i>1 & i<N){
        integ = gsub(pattern = "y", replacement = "x", x = integ)
        templete = paste(i, "*integrate(",integ,", (x, lob, y))")
        lob = qbeta(p = alpha,shape1 = i, shape2 = N-i+1, lower.tail = T)
        templete = gsub("lob", lob, templete)
      }else if(i==N)
      {
        integ = gsub(pattern = "y", replacement = "x", x=integ)
        templete = paste(i, "*integrate(",integ,", (x, lob, 1))")
        lob = qbeta(p = alpha,shape1 = i, shape2 = N-i+1, lower.tail = T)
        templete = gsub("lob", lob, templete)
      }
      #print(templete)
      integ = sympy(templete)
    }
    res = 1-as.numeric(integ)
    return(list(p=res, optimal_rank = idx_minimum, eff.p.idx = ord[1:idx_minimum], MMP = alpha))
  }
  if(is.onetail)
  {
    RES = ordmeta(p2 = p)
    return(RES)
  }else{
    p1 = p2 = p
    idx_pos = which(eff.sign >= 0)
    idx_neg = which(eff.sign < 0)
    p1[idx_pos] = p[idx_pos]/2
    p1[idx_neg] = 1-p[idx_neg]/2
    p2[idx_pos] = 1-p[idx_pos]/2
    p2[idx_neg] = p[idx_neg]/2

    RES1 = ordmeta(p2 = p1)
    RES2 = ordmeta(p2 = p2)
    if(RES1$p<=RES2$p){
      RES = RES1; RES$overall.eff.direction = "+"
    }else{
      RES = RES2; RES$overall.eff.direction = "-"
    }
    RES$p = RES$p * 2
    if(RES$p > 1.0){RES$p = 1.0}
    return(RES)
  }
}


p_value_adjustment <- function(p_vec){
  p_adj = 1 - (1 - min(p_vec))^length(p_vec)
  return(p_adj)
}

scale_p_values <- function(p_list, screen = 0.05){
    # sort the list of $p$-values
    p_sort_unscaled=sort(p_list)[2:length(p_list)]
    p_sort <- (p_sort_unscaled - screen)/(1-screen)
    return(p_sort)
}


order_stats_p_values_scaled <- function(p_list){
  # sort the list of $p$-values
  p_sort_unscaled <- sort(p_list)[2:length(p_list)]
  p_sort <- (p_sort_unscaled - min(p_list))/(1-min(p_list))
  #p_sort <- (p_sort_unscaled - min(p_list))/(1-min(p_list))
  n = length(p_sort)
  out_p = c()
  for (j in 1:length(p_sort)){
    # Grab the relevant $p$-value from the list
    p_j=p_sort[j]
    # Calculate the $p$-values of each of the order statistics (excluding the first gene)
    out_p = append(out_p, pbeta(p_j,j,n+1-j) )
  }
  return(out_p)
}


order_stats_p_values_unscaled <- function(p_list){
  # sort the list of $p$-values
  p_sort <- sort(p_list)[2:length(p_list)]
  #p_sort <- (p_sort_unscaled - min(p_list))/(1-min(p_list))
  #p_sort <- (p_sort_unscaled - min(p_list))/(1-min(p_list))
  n = length(p_sort)
  out_p = c()
  for (j in 1:length(p_sort)){
    # Grab the relevant $p$-value from the list
    p_j=p_sort[j]
    # Calculate the $p$-values of each of the order statistics (excluding the first gene)
    out_p = append(out_p, pbeta(p_j,j,n+1-j) )
  }
  return(out_p)
}


x<-function(){
## Let's try a Bayesian approach

# Bayes theorem:
# P(A|B) = P(B|A)P(A)/P(B)

# In our case, P(A|B) is what is the probability of the first order statistic given that
# we've screened on alpha_adap
#
# P(A) - the probability of the kth order stat
# P(B) - the probability of a gene being screened using alpha_adap
# P(B|A) - the probability that a gene was screened given that it was the kth order stat

# Here's a test case:

p_values <- c(0.001, 0.1, 0.12, 0.3, 0.4, 0.5, 0.8)

scaled_p_values <- scale_p_values(p_values)

# For the first order stat
alpha_adap <- 0.01

p_a <-  pbeta(p_values[2],1,6+1-1)
p_b <- 1-alpha_adap
p_ba <- 1-pbeta(alpha_adap,1,6+1-1)

p_ab <- (p_ba*p_a)/p_b

pbeta(0.1,1,6+1-1)
}




generateNullData <- function(adaptation_screen, a, b, n){
  temp <- c( rbeta(1,a,b),replicate(n-1, sample(10000,1)/10000) )
  while (sum(temp<adaptation_screen)==0){
    temp <- c( rbeta(1,a,b),replicate(n-1, sample(10000,1)/10000) )
  }
  return(temp)
}

generateAltData <- function(adaptation_screen, max, k, n){

  if (n-k ==0){
    temp <-  c(replicate(k, sample(max,1)/10000))
  }
  else if(k==0){
    temp<- c(replicate(n, sample(10000,1)/10000) )
  }
  else{
    temp <-  c(replicate(k, sample(max,1)/10000)
               ,c(replicate(n-k, sample(10000,1)/10000)) )
  }

  while (sum(temp<adaptation_screen)==0){
      if (n-k ==0){
        temp <-  c(replicate(k, sample(max,1)/10000))
      }
      else if(k==0){
        temp <- c(replicate(n, sample(10000,1)/10000) )
      }
      else{
        temp <-  c(replicate(k, sample(max,1)/10000)
                   ,c(replicate(n-k, sample(10000,1)/10000)) )
      }
  }
  return(temp)
}



runPowerAnalysis <- function(reps,  alpha_1, n, power_vec, speciesVec){
  emp_p_null_dat <- t(replicate(10000, generateNullData(alpha_1, 1, 1, n)))
  emp_p_null_dat_unscaled <- t(apply(emp_p_null_dat ,1, order_stats_p_values_unscaled))
  null_pMax_cor_unscaled <- cor( emp_p_null_dat_unscaled )

  result_count= 0
  output_DFs = list()

  for (GEA_power in c(0.05, 1:9/10) ){
    print(GEA_power)
    topHits = 500/GEA_power
    for (nSpecies in speciesVec){
      print(nSpecies)
      result_count = result_count +1
      test_data <- t(replicate(reps,
                               generateAltData(alpha_1,
                                               topHits,
                                               nSpecies,
                                               n)))

      test_data <- t(apply(test_data, 1, sort))
      test_data_order_stats_unscaled <- t(apply(test_data ,1, order_stats_p_values_unscaled))
      test_data_order_stats_unscaled

      binom_results = rep(-1, reps)
      binom_adj_results = rep(-1, reps)
      tippett_results = rep(-1, reps)
#      ordmeta_all_results = rep(-1, reps)
#      ordmeta_2plus_results = rep(-1, reps)


      for (i in 1:reps){
              print(i)
        #      print(test_data[i,])
        unscaled_p_value_vector = test_data_order_stats_unscaled[i,]

        temp_binom <-  binom.test( sum(test_data[i,]<alpha_1), n, alpha_1, , alternative='greater')$p.value
        temp_binom_adj <-  binom.test( sum(test_data[i,]<alpha_1)-1, n-1, alpha_1, , alternative='greater')$p.value
        temp_tippett <- tippett(unscaled_p_value_vector, adjust = "empirical", R = null_pMax_cor_unscaled, side = 1, size = 10000)$p
#        temp_ordmeta_all <- ordmeta(test_data[i,])$p
 #       temp_ordmeta_2plus <- ordmeta(sort(test_data[i,])[2:length(test_data[i,])])$p

        print(c(temp_tippett) )

        binom_results[i] = temp_binom
        binom_adj_results[i] = temp_binom_adj
        tippett_results[i] = temp_tippett
#        ordmeta_all_results[i] = temp_ordmeta_all
#        ordmeta_2plus_results[i] = temp_ordmeta_2plus

      }
      output_DFs[[result_count]] = data.frame(GEA_Power=GEA_power,
                                              nSpeciesConvergent=nSpecies,
                                              power_1_binomial=sum(binom_results<power_vec[1])/reps,
                                              power_1_binomial_adj=sum(binom_adj_results<power_vec[1])/reps,
                                              power_1_tippett=sum(tippett_results<power_vec[1])/reps,
                                           #   power_1_ordmeta_all=sum(ordmeta_all_results<power_vec[1])/reps,
                                           #   power_1_ordmeta_2plus=sum(ordmeta_2plus_results<power_vec[1])/reps,

                                              power_2_binomial=sum(binom_results<power_vec[2])/reps,
                                              power_2_binomial_adj=sum(binom_adj_results<power_vec[2])/reps,
                                              power_2_tippett=sum(tippett_results<power_vec[2])/reps)
                                          #    power_2_ordmeta_all=sum(ordmeta_all_results<power_vec[2])/reps,
                                          #    power_2_ordmeta_2plus=sum(ordmeta_2plus_results<power_vec[2])/reps)
    }

  }


  output_temp <- do.call(rbind, output_DFs)
  return(output_temp)

}


#################
### 7 species

Reps = 1000
Alpha_1 = 0.05
N = 7
Power_vec=c(0.01, 0.003757)
SpeciesVec = c(1:7)

#output_df_n7 <- runPowerAnalysis( Reps,  Alpha_1, N, Power_vec, SpeciesVec )
#write.csv(output_df_n7,"~/UBC/GEA/pMax/Analysis/7_species_alpha_0.05_adap.csv")
output_df_n7 <- read.csv("~/UBC/GEA/pMax/Analysis/7_species_alpha_0.05_adap.csv")

#output_df_n7_alpha_1_01 <- runPowerAnalysis( Reps,  0.01, N, Power_vec, SpeciesVec )
#write.csv(output_df_n7_alpha_1_01,"~/UBC/GEA/pMax/Analysis/7_species_alpha_0.01_adap.csv")
output_df_n7_alpha_1_01 <- read.csv("~/UBC/GEA/pMax/Analysis/7_species_alpha_0.01_adap.csv")

#################
## N = 30

Reps_n30 = 1000
Alpha_1_n30 = 0.01
N_n30 = 30
Power_vec_n30=c(0.01, 0.003318)
SpeciesVec_n30 = c(1,2,5,10,20,30)


#output_df_n30 <- runPowerAnalysis( Reps_n30,  Alpha_1_n30, N_n30, Power_vec_n30, SpeciesVec_n30 )
#write.csv(output_df_n30,"~/UBC/GEA/pMax/Analysis/30_species_alpha_0.01_adap.csv")
output_df_n30 <- read.csv("~/UBC/GEA/pMax/Analysis/30_species_alpha_0.01_adap.csv")
#output_df_n30$nSpeciesConvergent <- rep(c(1,2,5,10,20,30),10)
#################
## N = 30 - for sam

Reps_n30_s = 1000
Alpha_1_n30_s = 0.01
N_n30_s = 30
Power_vec_n30_s=c(0.03615, 0.003318)
SpeciesVec_n30_s = c(1,3,5,7,9,11)


#output_df_n30_s <- runPowerAnalysis( Reps_n30_s,  Alpha_1_n30_s, N_n30_s, Power_vec_n30_s, SpeciesVec_n30_s )
#write.csv(output_df_n30_s,"~/UBC/GEA/pMax/Analysis/30_species_alpha_0.01_adap_forSam.csv")
output_df_n30_s <- read.csv("~/UBC/GEA/pMax/Analysis/30_species_alpha_0.01_adap_forSam.csv")


#################
## N = 3

Reps_n3 = 1000
Alpha_1_n3 = 0.05
N_n3 = 3
Power_vec_n3=c(0.0073)
SpeciesVec_n3 = c(0,1,2,3)


#output_df_n3 <- runPowerAnalysis( Reps_n3,  Alpha_1_n3, N_n3, Power_vec_n3, SpeciesVec_n3 )
#write.csv(output_df_n3,"~/UBC/GEA/pMax/Analysis/3_species_alpha_0.05_adap.csv")
output_df_n3 <- read.csv("~/UBC/GEA/pMax/Analysis/3_species_alpha_0.05_adap.csv")




#################
## N = 2

Reps_n2 = 1000
Alpha_1_n2 = 0.05
N_n2 = 2
Power_vec_n2 = c(0.0073)
SpeciesVec_n2 = c(0,1,2)


output_df_n2 <- runPowerAnalysis( Reps_n2,  Alpha_1_n2, N_n2, Power_vec_n2, SpeciesVec_n2 )
write.csv(output_df_n2,"~/UBC/GEA/pMax/Analysis/2_species_alpha_0.05_adap.csv")
output_df_n2 <- read.csv("~/UBC/GEA/pMax/Analysis/2_species_alpha_0.05_adap.csv")

###########
# Plot the data
#
# Start with n=7
# With Binomial Test
output_df_n7$nSpeciesConvergent <- factor(output_df_n7$nSpeciesConvergent,
                                          levels = c(1,2,3,4,5,6,7),
                                           labels = c("1/7",
                                                      "2/7",
                                                      "3/7",
                                                      "4/7",
                                                      "5/7",
                                                      "6/7",
                                                      "7/7"))
power_1_plot_n7 <- ggplot(data = output_df_n7)+
  geom_hline(aes(yintercept = Power_vec[2]))+
  geom_line(aes(x = GEA_Power,
                y = power_2_tippett,
                col = as.factor(nSpeciesConvergent),
                lty = "PicMin"),
            lwd=1)+
   geom_line(aes(x = GEA_Power,
                 y = power_2_binomial,
                 col = as.factor(nSpeciesConvergent),
                 lty = "Binomial Test"),
             lwd = 1)+
#  ggtitle("Adaptation Outliers Identified\nUsing 95th Percentile")+
  scale_y_continuous(expression("Probability of Rejecting Null Hypothesis ("*alpha[Repeated]*"=0.01)"),
                     limits = c(0,1))+
  scale_x_continuous("Power of Genome Scan",
                     breaks = 1:9/10)+
  scale_colour_manual("Number of\nLineages\nWhere the\nGene is Causal", values = colorBlindBlack8)+
  scale_linetype_manual("Statistical Test", values = c(2,1))+
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5)
  )

### Now  n=30 with binomial
output_df_n30$nSpeciesConvergent <- factor(output_df_n30$nSpeciesConvergent,
                                           levels = c(1,2,5,10,20,30),
                                           labels =  c("1/30",
                                                      "2/30",
                                                      "5/30",
                                                      "10/30",
                                                      "20/30",
                                                      "30/30"))

power_1_plot_n30_adap_01 <- ggplot(data = output_df_n30)+
  geom_hline(aes(yintercept = Power_vec_n30[2]))+
  geom_line(aes(x = GEA_Power,
                y = power_2_tippett,
                col = as.factor(nSpeciesConvergent),
                lty = "PicMin"),
            lwd=1)+
   geom_line(aes(x = GEA_Power,
                 y = power_2_binomial,
                 col = as.factor(nSpeciesConvergent),
                 lty = "Binomial Test"),
             lwd = 1)+
#  ggtitle("Adaptation Outliers Identified\nUsing 99th Percentile")+
  scale_y_continuous(expression("Probability of Rejecting Null Hypothesis ("*alpha[Repeated]*"=0.00332)"),
                     limits = c(0,1))+
  scale_x_continuous("Power of Genome Scan",
                     breaks = 1:9/10)+
  scale_colour_manual("Number of\nLineages\nWhere the\nGene is Causal", values = colorBlindBlack8)+
  scale_linetype_manual("Statistical Test", values = c(2,1))+
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5)
  )


# Supp Figure 5
S5 <- ggarrange(power_1_plot_n7,
          power_1_plot_n30_adap_01,
          ncol = 1,
          nrow = 2,
          labels = "AUTO")

png("~/UBC/GEA/pMax/writeUp/Plots/SuppFigure5.png",width = 500, height = 700)
print(S5)
dev.off()


# Same plot, but without the Binomial Test

F1_power_1_plot_n7 <- ggplot(data = output_df_n7_alpha_1_01)+
  geom_hline(aes(yintercept = Power_vec[1]))+
  geom_line(aes(x = GEA_Power,
                y = power_2_tippett,
                col = as.factor(nSpeciesConvergent)),
            lwd=1)+
  scale_y_continuous(expression("Power ("*alpha[Repeated]*"=0.01)"),
                     limits = c(0,1))+
  scale_x_continuous("Power of Genome Scan",
                     breaks = 1:9/10)+
  scale_colour_manual("Number of\nLineages\nWhere the\nGene is Causal", values = colorBlindBlack8)+
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5)
  )

### Now  n=30 with binomial
F1_power_1_plot_n30_adap_01 <- ggplot(data = output_df_n30)+
  geom_hline(aes(yintercept = Power_vec_n30[1]))+
  geom_line(aes(x = GEA_Power,
                y = power_2_tippett,
                col = as.factor(nSpeciesConvergent)),
            lwd=1)+
  scale_y_continuous(expression("Power ("*alpha[Repeated]*"=0.01)"),
                     limits = c(0,1))+
  scale_x_continuous("Power of Genome Scan",
                     breaks = 1:9/10)+
  scale_colour_manual("Number of\nLineages\nWhere the\nGene is Causal", values = colorBlindBlack8)+
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5)
  )

# Figure 2
Figure2 <- ggarrange(F1_power_1_plot_n7,
                     F1_power_1_plot_n30_adap_01,
                ncol = 1,
                nrow = 2,
                labels = "AUTO")

pdf("~/UBC/GEA/pMax/writeUp/Plots/Figure2.pdf",width = 5, height = 7)
print(Figure2)
dev.off()


###################
### Num Species Convergent
###
###
###
###
###

n = 7
alpha_1 = 0.05

emp_p_null_dat <- t(replicate(10000, generateNullData(alpha_1, 1, 1, n)))
emp_p_null_dat_unscaled <- t(apply(emp_p_null_dat ,1, order_stats_p_values_unscaled))
null_pMax_cor_unscaled <- cor( emp_p_null_dat_unscaled )

reps = 200
result_count= 1
output_DFs = list()
threshold=  0.05

for (GEA_power in c(0.2, 0.4, 0.6, 0.8) ){
    topHits = 500/GEA_power
        for (nSpecies in c(3,5,7)){
#    for (nSpecies in c(3)){
    for (i in 1:reps){
        pass_test = FALSE

        while (pass_test == FALSE){

                  test_data <- generateAltData(alpha_1,
                        topHits,
                        nSpecies,
                        n)
                  test_data_order_stats_unscaled <- order_stats_p_values_unscaled(test_data)
                  temp_binom <-  binom.test( sum(test_data<alpha_1), n, alpha_1, , alternative='greater')$p.value
#                  temp_binom_adj <-  binom.test( sum(test_data[i,]<alpha_1)-1, n-1, alpha_1, , alternative='greater')$p.value
                  temp_tippett <- tippett(test_data_order_stats_unscaled, adjust = "empirical", R = null_pMax_cor_unscaled, side = 1)$p

                  if ((temp_binom<threshold)&(temp_tippett<threshold)){

                    binom_count = sum(test_data<alpha_1)
                    PicMin_estimate = which.min( test_data_order_stats_unscaled )+1
                    print(test_data)
                    print(c(GEA_power, nSpecies, i, binom_count, PicMin_estimate))
                    print("'")
                    output_DFs[[result_count]] = data.frame(GEA_Power=GEA_power,
                                                            nSpeciesConvergent=nSpecies,
                                                            binomial_estimate=binom_count,
                                                            PicMin_estimate=PicMin_estimate)
                    result_count = result_count +1
                    pass_test = TRUE


                  }
                }

    }
  }
}

#output_temp <- do.call(rbind, output_DFs)
#write.csv(output_temp, "~/UBC/GEA/pMax/Analysis/convergenceConfigurations_n7_a0.05.csv")
output_temp <- read.csv("~/UBC/GEA/pMax/Analysis/convergenceConfigurations_n7_a0.05.csv")

data <- melt(output_temp, id = c("GEA_Power", "nSpeciesConvergent", "X"))

data$variable  <- factor( data$variable,
                          levels = c("binomial_estimate",
                                     "PicMin_estimate"),
                          labels = c(expression("Simple Threshold ("*italic(ep)*"-value < 0.05)"),
                                     expression(italic("PicMin")*" Estimate"))
                          )


data$GEA_Power <- paste("Genome Scan Power = ", data$GEA_Power, sep = "")




safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499",
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

colourPal <- safe_colorblind_palette[c(3,7,11)]

data$indicator <- data$value == data$nSpeciesConvergent

tbl <- with(data, table(nSpeciesConvergent, variable, GEA_Power, value))

plot_data <- as.data.frame(tbl)

plot_data$indicator <- as.numeric(as.character(plot_data$value)) == as.numeric(as.character(plot_data$nSpeciesConvergent))
plot_data$nSpeciesConvergent_lab <- paste(plot_data$nSpeciesConvergent,  " Lineages With Causal Gene", sep = "")

parse.labels <- function(x) parse(text = x)


config_plot <- ggplot(plot_data, aes(value, Freq, fill = variable, col = indicator)) +
  geom_col(position = 'dodge')+
  scale_color_manual(values = c("white","black"))+

  scale_x_discrete("Estimate of the Number of Lineages Where the Gene is Causal",
                     breaks = c(1,2,3,4,5,6,7)) +
  scale_y_continuous(limits = c(0,200))+
  #coord_cartesian(clip = "off") +
  scale_fill_manual("", values = colourPal, labels = parse.labels)+
  guides(alpha = "none", colour = "none")+
  facet_grid(nSpeciesConvergent_lab~GEA_Power)+
  theme_bw()+
  theme(
    axis.title.y = element_blank() ,
    axis.title.x = element_text( hjust = 0.5),
    strip.background = element_blank(),
    strip.text.y = element_text(vjust = 1),
    legend.position = "bottom"
  )

ggsave("~/UBC/GEA/pMax/writeUp/Plots/Figure3.pdf",config_plot,
       width = 8,
       height = 7)

ggsave("~/UBC/GEA/pMax/writeUp/Plots/Figure3.png", config_plot,
       device = png,
       units = "in",
       width = 8, height = 7,
       res = 1000)


###########
#### Config plot n=30


n = 30
alpha_1 = 0.01

emp_p_null_dat <- t(replicate(10000, generateNullData(alpha_1, 1, 1, n)))
emp_p_null_dat_unscaled <- t(apply(emp_p_null_dat ,1, order_stats_p_values_unscaled))
null_pMax_cor_unscaled <- cor( emp_p_null_dat_unscaled )

reps = 200
result_count= 1
output_DFs = list()
threshold=  0.045

for (GEA_power in c(0.2, 0.4, 0.6, 0.8) ){
  topHits = 500/GEA_power
  for (nSpecies in c(3,5,7,9,11)){
    #    for (nSpecies in c(3)){
    for (i in 1:reps){
      pass_test = FALSE

      while (pass_test == FALSE){

        test_data <- generateAltData(alpha_1,
                                     topHits,
                                     nSpecies,
                                     n)
        test_data_order_stats_unscaled <- order_stats_p_values_unscaled(test_data)
        temp_binom <-  binom.test( sum(test_data<alpha_1), n, alpha_1, , alternative='greater')$p.value
        #                  temp_binom_adj <-  binom.test( sum(test_data[i,]<alpha_1)-1, n-1, alpha_1, , alternative='greater')$p.value
        temp_tippett <- tippett(test_data_order_stats_unscaled, adjust = "empirical", R = null_pMax_cor_unscaled, side = 1)$p

        if ((temp_binom<threshold)&(temp_tippett<threshold)){

          binom_count = sum(test_data<alpha_1)
          PicMin_estimate = which.min( test_data_order_stats_unscaled )+1
          print(test_data)
          print(c(GEA_power, nSpecies, i, binom_count, PicMin_estimate))
          print("'")
          output_DFs[[result_count]] = data.frame(GEA_Power=GEA_power,
                                                  nSpeciesConvergent=nSpecies,
                                                  binomial_estimate=binom_count,
                                                  PicMin_estimate=PicMin_estimate)
          result_count = result_count +1
          pass_test = TRUE


        }
      }

    }
  }
}

output_temp <- do.call(rbind, output_DFs)


data <- melt(output_temp, id = c("GEA_Power", "nSpeciesConvergent"))

data$variable  <- factor( data$variable,
                          levels = c("binomial_estimate",
                                     "PicMin_estimate"),
                          labels = c(expression(atop("Simple Threshold", paste("("*italic(ep)*"-value <0.05)"))),
                                     expression(atop(italic("PicMin"), paste("Estimate")))
                                     )
                          )

data$GEA_Power <- paste("Genome Scan Power = ", data$GEA_Power, sep = "")

expression(atop("Simple Threshold", paste("("*italic(ep)*"-value <0.05)")))

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499",
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

colourPal <- safe_colorblind_palette[c(3,7,11)]

data$indicator <- data$value == data$nSpeciesConvergent

tbl <- with(data, table(nSpeciesConvergent, variable, GEA_Power, value))

plot_data <- as.data.frame(tbl)

plot_data$indicator <- as.numeric(as.character(plot_data$value)) == as.numeric(as.character(plot_data$nSpeciesConvergent))
plot_data$nSpeciesConvergent_lab <- paste(plot_data$nSpeciesConvergent,  " Species Using Gene", sep = "")


plot_data$nSpeciesConvergent_lab <- factor(plot_data$nSpeciesConvergent_lab,
                         levels = c("3 Species Using Gene",
                                    "5 Species Using Gene",
                                    "7 Species Using Gene",
                                    "9 Species Using Gene",
                                    "11 Species Using Gene"))

config_plot <- ggplot(plot_data, aes(value, Freq, fill = variable)) +
  geom_col(position = 'dodge')+
  geom_col(data= plot_data[plot_data$indicator == TRUE,], aes(fill = variable),position = 'dodge',col ="black")+
  scale_x_discrete("Estimate of the Number of Lineages Using the Gene",
                   breaks = c(5,10,15,20,25,30)) +
  scale_y_continuous(limits = c(0,200))+
  #coord_cartesian(clip = "off") +
  scale_fill_manual("", values = colourPal)+
  guides(alpha = "none", colour = "none")+
  facet_grid(nSpeciesConvergent_lab~GEA_Power)+
  theme_bw()+
  theme(
    axis.title.y = element_blank() ,
    axis.title.x = element_text( hjust = 0.5),
    strip.background = element_blank(),
    strip.text.y = element_text(vjust = 1),
    legend.position = "bottom"
  )

ggsave("~/UBC/GEA/pMax/Analysis/convergenceConfigurations_n30_a0.05_t0.045_s.pdf",config_plot,
       width = 14,
       height = 10)


######
#
# Full genome comparison
#
#

nGenes <- 10000
nSpecies <- 7
nHits <- 100

# First, simulate p-values for each species
speciesEmpiricalP <- function(nGenes, nHits, a, b){
                              species_p_values <- c( rbeta( nHits, a,b),
                                                   runif(nGenes-nHits) )

                              species_empirical_p_values <- rank(species_p_values)/nGenes
                              return(species_empirical_p_values)
}


sum(rbeta(1000, 0.2,5)<0.05)

GEA_power_list <- list(
                list(a=0.8, b = 5, power=0.3), # ~30% power
#                list(a=0.66, b = 5, power=0.4), # ~40% power
                list(a=0.5, b = 5, power=0.5), # ~50% power
#                list(a=0.4, b = 5, power=0.6), # ~60% power
                list(a=0.3, b = 5, power = 0.7), # ~70% power
#                list(a=0.2, b = 5, power = 0.8), # ~70% power
                list(a=0.1, b = 5, power = 0.9) # ~90% power
)

alpha_a <- 0.05


# Run 10,000 replicate simulations of this situation and build the correlation matrix
emp_p_null_dat <- t(replicate(40000, generateNullData(alpha_a, 0.5, 3, 7)))

# Calculate the order statistics p-values for each simulation
emp_p_null_dat_unscaled <- t(apply(emp_p_null_dat ,1, order_stats_p_values_unscaled))

# Use those p-values to construct the correlation matrix
null_pMax_cor_unscaled <- cor( emp_p_null_dat_unscaled )


output_DF_list <- list()
result_count=0
for (k in c(3,5,7)){
#  for (k in c(7)){
    for (g in GEA_power_list){
    for (rep in 1:10){
        print(c(k,g$power,rep))
        result_count = result_count+1
        if (k < 7){
          dataSet <- cbind( c(replicate(k, speciesEmpiricalP(nGenes, nHits, g$a, g$b)) ),
                        c(replicate(7-k, speciesEmpiricalP(nGenes, 0, 1, 1) ) ))
        }
        else if (k == 7){
          dataSet <- cbind( replicate(k, speciesEmpiricalP(nGenes, nHits, g$a, g$b)) )
        }
                adap_screen_vector <- apply(dataSet<0.05, 1, sum)!=0
        gene_names <- c(1:nGenes)[adap_screen_vector]
        screenedDataSet <- dataSet[adap_screen_vector, ]

        registerDoParallel(numCores)

        output_p_values<- foreach (i=1:nrow(screenedDataSet), .combine=c) %dopar% {
#        output_p_values<- foreach (i=1:100, .combine=c) %dopar% {
          PicMin(screenedDataSet[i,], null_pMax_cor_unscaled)$p
        }

        q_values <- p.adjust( output_p_values, method = "fdr" )

        temp_df <- as.data.frame(screenedDataSet)
        temp_df$gene_id =  gene_names
        temp_df$q <- q_values

        true_positives <- sum((temp_df$q<0.05)&(temp_df$gene_id<=nHits))
        false_positives <- sum((temp_df$q<0.05)&(temp_df$gene_id>nHits))

        output_DF_list[[result_count]] = data.frame(GEA_Power=g$power,
                                                    rep = rep,
                                                nSpeciesConvergent=k,
                                                true_positives=true_positives,
                                                false_positives=false_positives)
        }
  }
}

temp <- do.call(rbind, output_DF_list)

#write.csv(temp, "~/UBC/GEA/pMax/Analysis/genome_analysis_3-7.csv")
temp <- read.csv("~/UBC/GEA/pMax/Analysis/genome_analysis_3-7.csv")


temp$GEA_Power <- as.character(temp$GEA_Power)
temp$nSpeciesConvergent <- as.character(temp$nSpeciesConvergent)

temp_melt <- melt(temp,
                  id = c("GEA_Power", "rep", "nSpeciesConvergent", "X"))


temp_melt$variable <- factor(temp_melt$variable,
                                levels = c("true_positives","false_positives"),
                                labels = c("True Positives", "False Positives"))

genome_fdr_plot <- ggplot(data = temp_melt, aes(x = as.factor(GEA_Power), y = value, fill = nSpeciesConvergent,col = nSpeciesConvergent))+
  geom_boxplot(alpha = 0.5)+
  scale_y_continuous("Number of Genes Identified\n(q<0.05)",
                     limits = c(0,100))+
  scale_x_discrete("Genome Scan Power")+
  facet_wrap(~variable)+
  scale_fill_manual("Number of\nLineages\nWhere the\nGene is Causal"
    ,values = colorBlindBlack8[c(3,5,7)])+
  scale_color_manual("Number of\nLineages\nWhere the\nGene is Causal"
                    ,values = colorBlindBlack8[c(3,5,7)])+
  theme_bw()+
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 12)
)


ggsave("~/UBC/GEA/pMax/writeUp/Plots/Plot4.pdf",
       genome_fdr_plot,
       width = 8,
       height = 4)
ggsave("~/UBC/GEA/pMax/writeUp/Plots/Plot4.png",
       genome_fdr_plot,
       units = "in",
       device = "png",
       width = 8,
       height = 4)




############
### Alpha_a threshold
###


runPowerAnalysis_alpha1 <- function(reps,  GEA_power_vec, n, power_vec, speciesVec){
  output_DFs = list()
  result_count= 0

  for (alpha_1 in c(10,50,100,500,1000)/10000 ){

    emp_p_null_dat <- t(replicate(10000, generateNullData(alpha_1, 1, 1, n)))
    emp_p_null_dat_unscaled <- t(apply(emp_p_null_dat ,1, order_stats_p_values_unscaled))
    null_pMax_cor_unscaled <- cor( emp_p_null_dat_unscaled )

    for (GEA_power in GEA_power_vec){

      print(GEA_power)
      topHits = 500/GEA_power
      for (nSpecies in speciesVec){
        print(nSpecies)
        result_count = result_count +1
        test_data <- t(replicate(reps,
                                 generateAltData(alpha_1,
                                                 topHits,
                                                 nSpecies,
                                                 n)))

        test_data <- t(apply(test_data, 1, sort))
        test_data_order_stats_unscaled <- t(apply(test_data ,1, order_stats_p_values_unscaled))
        test_data_order_stats_unscaled

  #      binom_results = rep(-1, reps)
  #      binom_adj_results = rep(-1, reps)
        tippett_results = rep(-1, reps)
        #      ordmeta_all_results = rep(-1, reps)
        #      ordmeta_2plus_results = rep(-1, reps)


        for (i in 1:reps){
          print(i)
          #      print(test_data[i,])
          unscaled_p_value_vector = test_data_order_stats_unscaled[i,]

  #        temp_binom <-  binom.test( sum(test_data[i,]<alpha_1), n, alpha_1, , alternative='greater')$p.value
  #        temp_binom_adj <-  binom.test( sum(test_data[i,]<alpha_1)-1, n-1, alpha_1, , alternative='greater')$p.value
          temp_tippett <- tippett(unscaled_p_value_vector, adjust = "empirical", R = null_pMax_cor_unscaled, side = 1, size = 10000)$p
          #        temp_ordmeta_all <- ordmeta(test_data[i,])$p
          #       temp_ordmeta_2plus <- ordmeta(sort(test_data[i,])[2:length(test_data[i,])])$p

          print(c(temp_tippett) )

  #        binom_results[i] = temp_binom
  #        binom_adj_results[i] = temp_binom_adj
          tippett_results[i] = temp_tippett
          #        ordmeta_all_results[i] = temp_ordmeta_all
          #        ordmeta_2plus_results[i] = temp_ordmeta_2plus

        }
        output_DFs[[result_count]] = data.frame(alpha_a=alpha_1,

                                                nSpeciesConvergent=nSpecies,
                                                GEA_Power=GEA_power,
                                      #          power_1_binomial=sum(binom_results<power_vec[1])/reps,
                                      #          power_1_binomial_adj=sum(binom_adj_results<power_vec[1])/reps,
                                                power_1_tippett=sum(tippett_results<power_vec[1])/reps,
                                                #   power_1_ordmeta_all=sum(ordmeta_all_results<power_vec[1])/reps,
                                                #   power_1_ordmeta_2plus=sum(ordmeta_2plus_results<power_vec[1])/reps,

                                      #          power_2_binomial=sum(binom_results<power_vec[2])/reps,
                                      #          power_2_binomial_adj=sum(binom_adj_results<power_vec[2])/reps,
                                                power_2_tippett=sum(tippett_results<power_vec[2])/reps)
        #    power_2_ordmeta_all=sum(ordmeta_all_results<power_vec[2])/reps,
        #    power_2_ordmeta_2plus=sum(ordmeta_2plus_results<power_vec[2])/reps)
        print( output_DFs[[result_count]])
      }



    }
  }
  output_temp <- do.call(rbind, output_DFs)
  return(output_temp)
}

#######
### 7 species

Reps = 1000
GEA_power_vec = c(0.2,0.50, 0.7, 0.9)
N = 7
Power_vec=c(0.05, 0.01)
SpeciesVec = c(1:7)

#output_df_n7_alpha_screen <- runPowerAnalysis_alpha1( Reps,  GEA_power_vec, N, Power_vec, SpeciesVec )
output_df_n7_alpha_screen <- read.csv("~/UBC/GEA/pMax/Analysis/falsePositives_alphaAdapt_n7.csv")
#write.csv(output_df_n7_alpha_screen,
#          "~/UBC/GEA/pMax/Analysis/falsePositives_alphaAdapt_n7.csv")

single_hit <- output_df_n7_alpha_screen[output_df_n7_alpha_screen$nSpeciesConvergent==1,]

single_hit_melt <- melt(single_hit, id = c("alpha_a","nSpeciesConvergent","GEA_Power", "X"))
#single_hit_melt <- melt(single_hit, id = c("alpha_a","nSpeciesConvergent","X"))


## Calculate Clopper-Pearson intervals
single_hit_melt$lower_p <-qbeta(0.05/2,
                                single_hit_melt$value*Reps,
                                Reps - single_hit_melt$value*Reps +1)

single_hit_melt$upper_p <-qbeta(1-0.05/2,
                                single_hit_melt$value*Reps+1,
                                Reps - single_hit_melt$value*Reps)

single_hit_melt$variable <- factor(single_hit_melt$variable,
                                   levels = c("power_1_tippett","power_2_tippett"),
                                   labels = c(expression(alpha[Repeated]*" = 0.05"),
                                              expression(alpha[Repeated]*" = 0.01")))


lineref_df <- data.frame(value = c(0.05,0.01), variable = c("power_1_tippett","power_2_tippett"))
lineref_df$variable <- factor(lineref_df$variable,
                                   levels = c("power_1_tippett","power_2_tippett"),
                                   labels = c(expression(alpha[Repeated]*" = 0.05"),
                                              expression(alpha[Repeated]*" = 0.01")))


##############
###  30 species
###

Reps = 1000
GEA_power_vec = c(0.2, 0.50, 0.7, 0.9)
N = 30
Power_vec=c(0.05, 0.01)
SpeciesVec = c(1)

#output_df_n30_alpha_screen <- runPowerAnalysis_alpha1( Reps,  GEA_power_vec, N, Power_vec, SpeciesVec )
output_df_n30_alpha_screen <- read.csv("~/UBC/GEA/pMax/Analysis/falsePositives_alphaAdapt_n30.csv")
#write.csv(output_df_n30_alpha_screen,
 #         "~/UBC/GEA/pMax/Analysis/falsePositives_alphaAdapt_n30.csv")

single_hit_n30 <- output_df_n30_alpha_screen[output_df_n30_alpha_screen$nSpeciesConvergent==1,]

single_hit_n30_melt <- melt(single_hit_n30, id = c("alpha_a","nSpeciesConvergent","GEA_Power", "X"))
#single_hit_melt <- melt(single_hit, id = c("alpha_a","nSpeciesConvergent","X"))


## Calculate Clopper-Pearson intervals
single_hit_n30_melt$lower_p <-qbeta(0.05/2,
                                single_hit_n30_melt$value*Reps,
                                Reps - single_hit_n30_melt$value*Reps +1)

single_hit_n30_melt$upper_p <-qbeta(1-0.05/2,
                                single_hit_n30_melt$value*Reps+1,
                                Reps - single_hit_n30_melt$value*Reps)

single_hit_n30_melt$variable <- factor(single_hit_n30_melt$variable,
                                   levels = c("power_1_tippett","power_2_tippett"),
                                   labels = c(expression(alpha[Repeated]*" = 0.05"),
                                              expression(alpha[Repeated]*" = 0.01")))

single_hit_n30_melt$n <- "30"
single_hit_melt$n <- "7"
plot_single_hit_data <- rbind(single_hit_melt,single_hit_n30_melt)
plot_single_hit_data$n <- factor( plot_single_hit_data$n ,
                                  levels = c("7","30"),
                                  labels = c(expression(7*" Species"),
                                             expression(30*" Species")))

alpha_adapt_FP_plot <- ggplot(data = plot_single_hit_data)+
               geom_point(aes( x= as.factor(GEA_Power),
                               y= value,
                               col = as.factor(alpha_a)),
                          position=position_dodge(width=0.25))+
  geom_errorbar(aes( x= as.factor(GEA_Power),
                     y= value,
                     col = as.factor(alpha_a),
                     ymin = lower_p,
                     ymax = upper_p),
                position=position_dodge(width=0.25),
                width = 0.01)+
  scale_x_discrete("Genome Scan Power")+
  scale_y_continuous("Proportion of False Positives",
                     limits =c(0,0.15))+
  geom_hline(data = lineref_df, aes(yintercept = value), lty = 2)+
  facet_grid(n~variable,
             scale = "free_y",
             label=label_parsed)+
  scale_color_brewer(expression(alpha[Adapt]),
                     palette="Dark2")+
  theme_bw()+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_text(size = 12)
  )

ggsave("~/UBC/GEA/pMax/Analysis/falsePositives_alphaAdapt.pdf", alpha_adapt_FP_plot,
       width = 7, height = 5)
ggsave("~/UBC/GEA/pMax/Analysis/falsePositives_alphaAdapt.png", device = png,
       alpha_adapt_FP_plot,
       units = "in",
       width = 7, height = 5,
       res = 1000)



multi_n_data <- output_df_n7_alpha_screen[output_df_n7_alpha_screen$nSpeciesConvergent>1,]

multi_n_data_melt <- melt(multi_n_data,
                          id = c("X","alpha_a","nSpeciesConvergent","GEA_Power"))
label_vec <- c(expression(2*" Species Convergent"),
  expression(3*" Species Convergent"),
  expression(4*" Species Convergent"),
  expression(5*" Species Convergent"),
  expression(6*" Species Convergent"),
  expression(7*" Species Convergent"))

multi_n_data_melt$variable <- factor(multi_n_data_melt$variable,
                                       levels = c("power_1_tippett","power_2_tippett"),
                                       labels = c("alpha[Repeated]*' = 0.05'",
                                                  "alpha[Repeated]*' = 0.01'"))

multi_n_data_melt$n <-  factor(multi_n_data_melt$nSpeciesConvergent,
                               levels = 2:7,
                               labels = label_vec)
  #paste(multi_n_data_melt$nSpeciesConvergent, "Species Convergent")

alpha_adapt_TP_plot <- ggplot(data = multi_n_data_melt)+
  geom_line(aes( x= as.factor(GEA_Power),
                  y= value,
                  col = as.factor(alpha_a),
                 group = as.factor(alpha_a)),
             lwd= 1)+
  facet_grid(n~variable,
             labeller = label_parsed)+
  scale_x_discrete("Genome Scan Power")+
  scale_y_continuous("Probability of Rejecting Null Hypothesis")+
  scale_color_brewer(expression(alpha[Adapt]),
                     palette="Dark2")+
  theme_bw()+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_text(size = 12)
  )

ggsave("~/UBC/GEA/pMax/Analysis/truePositives_alphaAdapt.pdf", alpha_adapt_TP_plot,
       width = 7, height = 12)
ggsave("~/UBC/GEA/pMax/Analysis/truePositives_alphaAdapt.png", alpha_adapt_TP_plot,
       device = png,
       units = "in",
       width = 7, height = 12,
       res = 1000)


###################
# p-value distributions under the null




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

replicates = 100
n=7
results_DFs <- list()
count = 0

for (alpha_adapt in c(0.005, 0.01, 0.05, 1.0)){
  null_data_for_mat <- t(replicate(40000, GenerateNullData(alpha_adapt, 0.5, 3, n, 10000)))
  
  emp_p_null_dat_unscaled <- t(apply(null_data_for_mat,
                                     1,
                                     order_stats_p_values_unscaled))
  
  # Use those p-values to construct the correlation matrix
  null_pMax_cor_unscaled <- cor( emp_p_null_dat_unscaled )
  for (maxP in c(600, 1000, 2000)){
    count = count + 1
    print(c(alpha_adapt, maxP))


  null_data <- t( replicate( replicates,
                             generateAltData(alpha_adapt,
                                                  maxP,
                                                 1,
                                                 n) ) )
  null_data_order_stats <- t(apply(null_data_for_mat,
                                   1,
                                   order_stats_p_values_unscaled))

  result_vector = rep(0,
                      replicates)
  print("!")
  for (s in seq_len(replicates)){
    print(s)
    temp_tippett <- tippett(null_data_order_stats[s,],
                            adjust = "empirical",
                            R = null_pMax_cor_unscaled,
                            side = 1,
                            size = 10000)$p

    result_vector[s] = temp_tippett
  }
  results_DFs[[count]] = data.frame(p_vals = result_vector,
                                     alpha_adapt = rep(alpha_adapt, replicates),
                                     maxP = rep(maxP, replicates))
  }
}

histo_data <- do.call( rbind, results_DFs )

histo_data$alpha_adapt <- factor(histo_data$alpha_adapt,
                                 levels = c(0.005,0.01,0.05,1.0),
                                 labels = c(expression(alpha[Adapt]*"=0.005"),
                                            expression(alpha[Adapt]*"=0.01"),
                                            expression(alpha[Adapt]*"=0.05"),
                                            expression(alpha[Adapt]*"=1.0")))

histo_data$maxP <- 1/(histo_data$maxP/500)
histo_data$scan_power <- factor(histo_data$maxP,
                                levels = c(5/6, 0.5, 0.25),
                                labels = c(expression("Genome scan power"*" = 0.8333..."),
                                           expression("Genome scan power"*" = 0.5"),
                                           expression("Genome scan power"*" = 0.25")))

null_hypothesis_histogram <- ggplot(data = histo_data,
       aes(x = p_vals))+
  geom_histogram( binwidth = 0.05, boundary  = 1, fill = "grey", col = "black" )+
  facet_grid(scan_power~alpha_adapt,
             label = label_parsed)+
  scale_x_continuous(expression(italic("p")*"-value"))+
#  scale_y_continuous("Count",
#                     limits = c(0,700))+
  theme_bw()+
    theme(
      strip.text.x = element_text(size = 13),
      strip.background = element_blank()
    )

ggsave("~/UBC/GEA/pMax/Analysis/null_distribution_histogram_n7.pdf",
       width = 9.50,
       height = 6.50,
        units = "in")



replicates = 10000
n=30
results_DFs <- list()
count = 0

for (alpha_adapt in c(0.005, 0.01, 0.05, 1.0)){
  for (maxP in c(600, 1000, 2000)){
    count = count + 1
    null_data_for_mat <- t(replicate(40000, GenerateNullData(alpha_adapt, 0.5, 3, n, 10000)))

    emp_p_null_dat_unscaled <- t(apply(null_data_for_mat,
                                       1,
                                       order_stats_p_values_unscaled))

    # Use those p-values to construct the correlation matrix
    null_pMax_cor_unscaled <- cor( emp_p_null_dat_unscaled )

    null_data <- t( replicate( replicates,
                               generateAltData(alpha_adapt,
                                               maxP,
                                               1,
                                               n) ) )
    null_data_order_stats <- t(apply(null_data_for_mat,
                                     1,
                                     order_stats_p_values_unscaled))

    result_vector = rep(0,
                        replicates)

    for (s in seq_len(replicates)){
      temp_tippett <- tippett(null_data_order_stats[s,],
                              adjust = "empirical",
                              R = null_pMax_cor_unscaled,
                              side = 1,
                              size = 10000)$p

      result_vector[s] = temp_tippett
    }
    results_DFs[[count]] = data.frame(p_vals = result_vector,
                                      alpha_adapt = rep(alpha_adapt, replicates),
                                      maxP = rep(maxP, replicates))
  }
}

histo_data <- do.call( rbind, results_DFs )

histo_data$alpha_adapt <- factor(histo_data$alpha_adapt,
                                 levels = c(0.005,0.01,0.05,1.0),
                                 labels = c(expression(alpha[Adapt]*"=0.005"),
                                            expression(alpha[Adapt]*"=0.01"),
                                            expression(alpha[Adapt]*"=0.05"),
                                            expression(alpha[Adapt]*"=1.0")))

histo_data$maxP <- 1/(histo_data$maxP/500)
histo_data$scan_power <- factor(histo_data$maxP,
                                levels = c(5/6, 0.5, 0.25),
                                labels = c(expression("Genome scan power"*" = 0.8333..."),
                                           expression("Genome scan power"*" = 0.5"),
                                           expression("Genome scan power"*" = 0.25")))



null_hypothesis_histogram <- ggplot(data = histo_data,
                                    aes(x = p_vals))+
  geom_histogram( binwidth = 0.05, boundary  = 1, fill = "grey", col = "black" )+
  facet_grid(scan_power~alpha_adapt,
             label = label_parsed)+
  scale_x_continuous(expression(italic("p")*"-value"))+
  scale_y_continuous("Count",
                     limits = c(0,700))+
  theme_bw()+
  theme(
    strip.text.x = element_text(size = 13),
    strip.background = element_blank()
  )

ggsave("~/UBC/GEA/pMax/Analysis/null_distribution_histogram_n30.pdf",
       width = 9.50,
       height = 6.50,
       units = "in")
