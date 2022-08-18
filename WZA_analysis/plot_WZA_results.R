rm(list = ls())

library(ggplot2)
library(ggpubr)

wza <- read.csv("~/UBC/GEA/pMax/WZA_analysis/all_results_r10000.7sets.FixMin099.csv")

wza[wza$p<0.01,]

results <- list()
FDR_results <- list()


all_positives <- wza[wza$q<0.05,]

all_positives$diff <- abs(all_positives$n_est - all_positives$numLA)

hist(all_positives$diff)

FDR = nrow( all_positives[all_positives$numLA<2,] ) / nrow(all_positives)
conf_int_fdr <- binom.test( nrow( all_positives[all_positives$numLA<2,] ) ,nrow(all_positives), 0.05)$conf.int
FDR_results[[1]] = c(-1, FDR, conf_int_fdr[1], conf_int_fdr[2])

n_wza <- wza[ (wza$numLA<2), ]
power = sum( n_wza$q < 0.05 )/nrow(n_wza)
conf_int <- binom.test( sum( n_wza$q < 0.05 ), nrow(n_wza), 0.05)$conf.int
results[[1]] = c(0, power, conf_int[1], conf_int[2])

for (i in (1:6)){
  n_wza <- wza[ (wza$numLA==i+1), ]
  power = sum( n_wza$q < 0.05 )/nrow(n_wza)
  conf_int <- binom.test( sum( n_wza$q < 0.05 ), nrow(n_wza), 0.05)$conf.int
  results[[i+1]] = c(i+1, power, conf_int[1], conf_int[2])
}

fdr_results_df <- as.data.frame(do.call(rbind,FDR_results))
names( fdr_results_df ) <- c("numSignif",
                            "power",
                            "power_up",
                            "power_down")

power_results <- as.data.frame(do.call(rbind,results))
names( power_results ) <- c("numSignif",
                            "power",
                            "power_up",
                            "power_down")
power_results$facetter <- power_results$numSignif==-1
power_results$facetter <- factor(power_results$facetter,
                                 levels = c(FALSE),
                                 labels = c("Probability of Rejecting Null Hypothesis (q<0.05)"))

power_results$numSignif <- factor(power_results$numSignif,
                                  levels = c(0,2,3,4,5,6,7),
                                  labels = c("0/1\n(False Positive Rate)","2","3","4","5","6","7"))


library(ggplot2)

panel_power <- ggplot(data = power_results,
       aes(x = as.factor(numSignif),
           y = power))+
  geom_point(size = 3)+
  geom_errorbar(aes(ymax = power_up,
                    ymin = power_down),
                width = 0.3)+
  ggtitle("FixMin")+
    scale_y_continuous(limits = c(0,1))+
  ylab("Probability of Rejecting Null Hypothesis (q<0.05)")+
  xlab("Number of Lineages Where the Gene is Causal")+
  theme_bw()+
  scale_color_brewer(expression(alpha[Adapt]),
                     palette="Dark2")+
  theme_bw()+
  theme(
    strip.background = element_blank(),
    axis.title = element_text(size = 12, color = "black"),
    strip.text.x = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
  )



panel_FDR <- ggplot(data = fdr_results_df,
       aes(x = as.factor(numSignif),
           y = power))+
  geom_point(size = 3)+
  geom_errorbar(aes(ymax = power_up,
                    ymin = power_down),
                width = 0.3)+
  ggtitle("FixMin")+
  geom_hline( aes (yintercept = 0.05), lty = 2)+
  scale_y_continuous( limits = c(0,1))+
  ylab("False Dicovery Rate")+
  xlab("")+
  theme_bw()+
  scale_color_brewer(expression(alpha[Adapt]),
                     palette="Dark2")+
  theme_bw()+
  theme(
    strip.background = element_blank(),
    axis.title = element_text(size = 12, color = "black"),
    strip.text.x = element_text(size = 12, color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 10, color = "black"),
  )

combined_figure <- ggarrange(panel_power,
          panel_FDR,
          widths = c(3.5,1),
          align ="h",
          labels = "AUTO")

pdf("~/UBC/GEA/pMax/writeUp/Plots/WZA_results_plot.FixMin.pdf",
    width = 8.74,
    height = 5.16)
print(combined_figure)
dev.off()
