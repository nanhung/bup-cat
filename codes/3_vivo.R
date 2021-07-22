# Load packages ----------------------------------------------------------------
library(data.table)
library(rstan)
library(tidyverse)
library(bayesplot)
library(GGally)
library(ggpubr)
library(scales)
library(cowplot)
library(magrittr)
library(LaplacesDemon)
library(bayestestR)
library(PowerTOST)
library(foreach)
library(doParallel) 
library(PKNCA)
library(kableExtra)
source("MCSim/function.R")

# Posterior check --------------------------------------------------------------

sim1.1 <- fread("outputs/bupropion_MCMC_vivo_1122.1.out") 
sim2.1 <- fread("outputs/bupropion_MCMC_vivo_2127.1.out") 
sim3.1 <- fread("outputs/bupropion_MCMC_vivo_3365.1.out") 
sim4.1 <- fread("outputs/bupropion_MCMC_vivo_4880.1.out") 
sim5.1 <- fread("outputs/bupropion_MCMC_vivo_5966.1.out") 

sim1.2 <- fread("outputs/bupropion_MCMC_vivo_1122.2.out") 
sim2.2 <- fread("outputs/bupropion_MCMC_vivo_2127.2.out") 
sim3.2 <- fread("outputs/bupropion_MCMC_vivo_3365.2.out") 
sim4.2 <- fread("outputs/bupropion_MCMC_vivo_4880.2.out") 
sim5.2 <- fread("outputs/bupropion_MCMC_vivo_5966.2.out") 

SIM1 <- do.call(rbind, list(sim1.1, sim1.2))
SIM2 <- do.call(rbind, list(sim2.1, sim2.2))
SIM3 <- do.call(rbind, list(sim3.1, sim3.2))
SIM4 <- do.call(rbind, list(sim4.1, sim4.2))
SIM5 <- do.call(rbind, list(sim5.1, sim4.2))

x <-mcmc_array(list(SIM1, SIM2, SIM3, SIM4, SIM5))
#dim(x)

pars_name <- dimnames(x)[[3]]
str <- which(pars_name == "G_Radius(1)")
end <- which(pars_name == "Weibull_slope_ER(1)")
M_pars <- pars_name[str:end]
str <- which(pars_name == "Vr_G_Radius(1)")
end <- which(pars_name == "Vr_Weibull_slope_ER(1)")
Vr_pars <- pars_name[str:end]
str <- which(pars_name == "Va_Peff(1)")
end <- which(pars_name == "Va_Kp2c(1)")
Va_pars <- pars_name[str:end]
mnt <- monitor(x[,,c(M_pars, Vr_pars, Va_pars, "Ve_C_central(1)")], digit=4)

# Create mcmc report sheet
Parameter <- c(M_pars, Vr_pars, Va_pars, "Ve_C_central(1)") %>%
  strsplit(split =  "\\(1\\)") %>% unlist()
Q5 <- mnt$Q5 %>% round(digits = 4)
Q50 <- mnt$Q50 %>% round(digits = 4)
Q95 <- mnt$Q95 %>% round(digits = 4)
Mean <- mnt$mean %>% round(digits = 4)
SD <- mnt$sd %>% round(digits = 4)
n_eff <- mnt$n_eff
Rhat <- mnt$Rhat %>% round(digits = 4)
data.frame(Parameter, Q5, Q50, Q95, Mean, SD, n_eff, Rhat) %>%
  write.csv(file="outputs/mcmc_report_2.csv", row.names=FALSE)

# Probability density - Population mean, inter- & intra-individual variability
j <- 2502:4926
theme_set(theme_bw())

mcmc_dens_overlay(x[j, ,], pars = c(M_pars), facet_args = list(ncol = 4, strip.position = "top")) +
  scale_color_discrete()
ggsave("plots/suppl/population_mean_posterior.pdf", dpi = 300, height = 8, width = 11, units="in")

mcmc_dens_overlay(x[j, ,], pars = c(Vr_pars), facet_args = list(ncol = 4, strip.position = "top")) +
  scale_color_discrete()
ggsave("plots/suppl/population_vr_posterior.pdf", dpi = 300, height = 8, width = 11, units="in")

mcmc_dens_overlay(x[j, ,], pars = c(Va_pars), facet_args = list(ncol = 4, strip.position = "top")) +
  scale_color_discrete()
ggsave("plots/suppl/population_va_posterior.pdf", dpi = 300, height = 4, width = 11, units="in")

# Weibull parameter release - individuals  
k <- 1:15
Weibull_scale_IR <- paste0("Weibull_scale_IR(1.", k, ")")
Weibull_scale_SR <- paste0("Weibull_scale_SR(1.", k, ")")
Weibull_scale_ER <- paste0("Weibull_scale_ER(1.", k, ")")
Weibull_slope_ER <- paste0("Weibull_slope_ER(1.", k, ")")

# Correlation matrix
sum_chains <- length(j)*5

x[j,,M_pars] %>% matrix(nrow = sum_chains) %>%
  `colnames<-`(Parameter[1:13]) %>%
  ggcorr(nbreaks = 5, label = T, label_size = 3, label_round = 2, size = 2.5, label_alpha = TRUE) +
  theme_minimal()
ggsave("plots/suppl/correlation_matrix.pdf", dpi = 300, height = 10, width = 14, units="in")


# Fig 6 -------------------------------------------------------------------

ir.scale.x <- x[,,c("Weibull_scale_IR(1)", paste0("Weibull_scale_IR(1.", c(1:15),")"))]
sr.scale.x <- x[,,c("Weibull_scale_SR(1)", paste0("Weibull_scale_SR(1.", c(1:15),")"))]
er.scale.x <- x[,,c("Weibull_scale_ER(1)", paste0("Weibull_scale_ER(1.", c(1:15),")"))]
er.slope.x <- x[,,c("Weibull_slope_ER(1)", paste0("Weibull_slope_ER(1.", c(1:15),")"))]
dimnames(ir.scale.x)[[3]] <- c("\u03bc", paste0("\u03b8", 1:15)) #Unicode characters

# Fig. 6B
theme_set(theme_pubclean())
p6.1 <- ir.scale.x %>% mcmc_intervals() +
  labs(title="IR scale") + 
  geom_vline(xintercept = 0.17, lty=3, lwd=1) +
  theme(axis.text = element_text (size = 12, colour = "black", face = "bold"))
p6.2 <- sr.scale.x %>% mcmc_intervals() +
  geom_vline(xintercept = 3, lty=3, lwd=1) +
  labs(title="SR scale") + theme(axis.text.y = element_blank()) +
  theme(axis.text = element_text (size = 12, colour = "black", face = "bold"))
p6.3 <- er.scale.x %>% mcmc_intervals() +
  labs(title="ER scale") + theme(axis.text.y = element_blank()) +
  geom_vline(xintercept = 3, lty=3, lwd=1) +
  theme(axis.text = element_text (size = 12, colour = "black", face = "bold"))
p6.4 <- er.slope.x %>% mcmc_intervals() +
  geom_vline(xintercept = 1.8, lty=3, lwd=1) +
  labs(title="ER slope") + theme(axis.text.y = element_blank()) +
  theme(axis.text = element_text (size = 12, colour = "black", face = "bold"))


# Random simulations = 20
testfile <- "MCSim/bupropion_MCMC_ppc.in"

check_post <- function(x, iter, testfile){
  d <- x[j,,] %>% matrix(nrow = 12125) %>%
    as.data.frame() %>% `colnames<-`(dimnames(x)[[3]])
  n <- sample(1:nrow(d), iter)
  d[n,] %>% write.table(file="testo.out", row.names=FALSE, sep="\t")
  for(i in 1:iter){ # last i posterior fit
    data.table::fread("testo.out") %>%
      as.data.frame() %>% tail(i) %>% head(1) %>% 
      write.table(file="testn.out", row.names=FALSE, sep="\t")
    system(paste0("./mcsim.bupropion_CAT.model.exe ", testfile))
    df <- read.delim("check.out")
    df$iter <- i
    if (i == 1) DF <- df else DF <- rbind(DF, df)
  }
  file.remove(c("check.out", "testn.out", "testo.out"))
  return(DF)
}

DF <- check_post(x, 20, testfile)

# Tidy data
all <- DF %>%
  mutate(Data1 = ifelse(Data == -1, NA, Data)) %>%
  na.omit() %>%
  separate(col = Level, c("Pop", "Subj", "Exp")) %>%
  group_by(Subj, Simulation, Time, Data) %>% 
  mutate(Pred=Prediction, select="", ratio = Pred/Data1) %>%
  mutate(Formulation = ifelse(Exp == 1, "IR075", 
                              ifelse(Exp == 2, "IR100",
                                     ifelse(Exp == 3, "SR100",
                                            ifelse(Exp == 4, "SR150",
                                                   ifelse(Exp == 5, "ER150", "ER300"))))))

global.dat <- DF %>% mutate(Data1 = ifelse(Data == -1, NA, Data)) %>%
  na.omit() %>%
  separate(col = Level, c("Pop", "Subj", "Exp")) %>%
  mutate(Formulation = ifelse(Exp == 1, "IR075",
                              ifelse(Exp == 2, "IR100",
                                     ifelse(Exp == 3, "SR100",
                                            ifelse(Exp == 4, "SR150",
                                                   ifelse(Exp == 5, "ER150", "ER300")))))) %>%
  mutate(R = Prediction/Data) %>%
  group_by(Subj, Simulation, Formulation, Time, Data) %>% 
  summarize(ratio = median(R)) %>%
  mutate(select = ifelse(ratio > 5 | ratio < 0.2, Subj, ""))

global.dat$Formulation <- factor(global.dat$Formulation, 
                                 level = c("IR075", "IR100", "SR100",
                                           "SR150", "ER150", "ER300"))

# Fig. 6B
p6.5 <- global.dat %>% ggplot(aes(Data, ratio, label = select)) +    
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(lim = c(1e-3,1e3),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", scales::math_format(10^.x))) +
  geom_hline(yintercept = 0.5,lty = 3, lwd = 1) +
  geom_hline(yintercept = 2, lty = 3, lwd = 1) +
  geom_hline(yintercept = 1, lty = 1, col = "grey", lwd=2) +
  geom_hline(yintercept = 5, lty = 2, lwd = 1) +
  geom_hline(yintercept = 0.2, lty = 2, lwd = 1) +
  geom_text(hjust = - 0.5, size = 5, color = "grey") +
  geom_point(data = all, color = "grey", alpha=0.2) +
  geom_point(aes(group = as.factor(Data), shape = Formulation)) +
  geom_smooth(method=lm, se=FALSE) +
  theme_pubclean() +
  labs(x = "Observed (ng/mL)", y = "Expected/Observed") +
  theme (axis.text = element_text (size   = 15,
                                   colour = "black", face = "bold"),
         axis.title = element_text (size   = 18,
                                    colour = "black", face = "bold"),
         legend.background=element_blank(), 
         legend.position=c(0.8, 0.2))

all$Formulation <- factor(all$Formulation,
                          level = c("IR075", "IR100", "SR100",
                                    "SR150", "ER150", "ER300"))

# Fig. 6C
p6.6 <- all %>% ggplot(aes(x=Formulation, y=ratio)) +
  geom_hline(yintercept = 1, col = "grey", lwd=2) +
  geom_hline(yintercept = 5, lty = 2, lwd = 1) +
  geom_hline(yintercept = 0.2, lty = 2, lwd = 1) +
  geom_hline(yintercept = 0.5, lty = 3, lwd = 1) +
  geom_hline(yintercept = 2, lty = 3, lwd = 1) +
  geom_violin(alpha=0.5) +
  geom_boxplot(width = 0.2, alpha=0.5, outlier.colour = NA) +
  scale_y_log10(lim = c(1e-3,1e3),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_pubclean() +
  theme(axis.text.y=element_blank(),
        axis.text.x = element_text (size = 12,
                                    colour = "black", face = "bold"), 
        axis.title=element_blank(),
        strip.text = element_blank(),
        legend.text=element_blank(),
        legend.title = element_blank(),
        legend.position="none")


# Accuracy (ratio within the range of 0.5-2)
all %>% as.data.frame() %>% 
  filter(Formulation == "IR075") %>% mutate(pass = ifelse(ratio > 2 | ratio < 0.5, 0, 1)) %>% summarise(accuracy = sum(pass)/length(pass))
all %>% as.data.frame() %>% 
  filter(Formulation == "IR100") %>% mutate(pass = ifelse(ratio > 2 | ratio < 0.5, 0, 1)) %>% summarise(accuracy = sum(pass)/length(pass))
all %>% as.data.frame() %>% 
  filter(Formulation == "SR100") %>% mutate(pass = ifelse(ratio > 2 | ratio < 0.5, 0, 1)) %>% summarise(accuracy = sum(pass)/length(pass))
all %>% as.data.frame() %>% 
  filter(Formulation == "SR150") %>% mutate(pass = ifelse(ratio > 2 | ratio < 0.5, 0, 1)) %>% summarise(accuracy = sum(pass)/length(pass))
all %>% as.data.frame() %>% 
  filter(Formulation == "ER150") %>% mutate(pass = ifelse(ratio > 2 | ratio < 0.5, 0, 1)) %>% summarise(accuracy = sum(pass)/length(pass))
all %>% as.data.frame() %>% 
  filter(Formulation == "ER300") %>% mutate(pass = ifelse(ratio > 2 | ratio < 0.5, 0, 1)) %>% summarise(accuracy = sum(pass)/length(pass))

# Precision (differene between highr and lower limit)
all %>% as.data.frame() %>% 
  filter (Formulation == "IR075") %>% summarise(q_ll = min(ratio), q_ul = max(ratio)) %>% 
  mutate(diff = q_ul - q_ll)
all %>% as.data.frame() %>% 
  filter (Formulation == "IR100") %>% summarise(q_ll = min(ratio), q_ul = max(ratio)) %>% 
  mutate(diff = q_ul - q_ll)
all %>% as.data.frame() %>% 
  filter (Formulation == "SR100") %>% summarise(q_ll = min(ratio), q_ul = max(ratio)) %>% 
  mutate(diff = q_ul - q_ll)
all %>% as.data.frame() %>% 
  filter (Formulation == "SR150") %>% summarise(q_ll = min(ratio), q_ul = max(ratio)) %>% 
  mutate(diff = q_ul - q_ll)  
all %>% as.data.frame() %>% 
  filter (Formulation == "ER150") %>% summarise(q_ll = min(ratio), q_ul = max(ratio)) %>% 
  mutate(diff = q_ul - q_ll)  
all %>% as.data.frame() %>% 
  filter (Formulation == "ER300") %>% summarise(q_ll = min(ratio), q_ul = max(ratio)) %>% 
  mutate(diff = q_ul - q_ll)  

#tiff("plots/fig6.tiff", res=600, compression = "lzw", height = 12, width = 15, units="in")
pdf("plots/fig6.pdf", height = 12, width = 15)
plot_grid(
  plot_grid(p6.1, p6.2, plot_grid(p6.3, p6.4, nrow=1), ncol = 3, labels = c('A'), label_size = 22),
  plot_grid(p6.5, p6.6, align = 'h', axis = 'b', rel_widths = c(3/5,2/5), labels = c('B', 'C'), label_size = 22),
  nrow = 2, rel_heights = c(2/5, 3/5)
)
dev.off()
#ggsave("plots/fig6.jpeg", dpi = 300, height = 12, width = 15, units="in")



calib_plot <- function(data, global = F, ratio = F, Exper = 1, ...){
  
  d <- exp(log(data$Data))
  d <- d[!is.na(d)]
  
  rng <- range(d)
  all <- data %>% mutate(Data1 = ifelse(Data == -1, NA, Data)) %>% na.omit() %>%
    separate(col = Level, c("Pop", "Subj", "Exp")) %>%
    group_by(Subj, Simulation, Time, Data) %>% 
    mutate(Pred=Prediction, select="", ratio = Pred/Data1) %>%
    mutate(Formulation = ifelse(Exp == 1, "IR075", 
                                ifelse(Exp == 2, "IR100",
                                       ifelse(Exp == 3, "SR100",
                                              ifelse(Exp == 4, "SR150",
                                                     ifelse(Exp == 5, "ER150", "ER300"))))))
  
  qtl <- data %>% separate(col = Level, c("Pop", "Subj", "Exp")) %>% 
    filter(Exp == Exper) %>% group_by(Time, Subj) %>% 
    summarise(Data = mean(Data), Prediction = median(Prediction),
              UCL = quantile(Prediction, 0.95), LCL = quantile(Prediction, 0.05))  
  
  data %>% separate(col = Level, c("Pop", "Subj", "Exp")) %>% 
    filter(Exp == Exper) %>%
    ggplot(aes(Time, Prediction)) + 
    facet_wrap(~as.numeric(Subj), ncol=3) + 
    #scale_y_log10() + 
    scale_x_log10(limits = c(0.1,99)) +
    geom_line(aes(Time, Prediction, group = iter), color = "grey") + 
    geom_point(data = qtl, aes(Time, exp(log(Data)))) + 
    geom_line(data = qtl, aes(Time, Prediction), color = "black") +
    labs(x = "time", y = "conc") +
    theme(legend.position="none") +
    labs(...)
}

p1 <- calib_plot(DF, Exper = 1, x = "Time (hr)", y = "Conc (ng/ml)", title = "IR 75")
p2 <- calib_plot(DF, Exper = 2, x = "Time (hr)", y = "Conc (ng/ml)", title = "IR 100")
p3 <- calib_plot(DF, Exper = 3, x = "Time (hr)", y = "Conc (ng/ml)", title = "SR 100")
p4 <- calib_plot(DF, Exper = 4, x = "Time (hr)", y = "Conc (ng/ml)", title = "SR 150")
p5 <- calib_plot(DF, Exper = 5, x = "Time (hr)", y = "Conc (ng/ml)", title = "ER 150")
p6 <- calib_plot(DF, Exper = 6, x = "Time (hr)", y = "Conc (ng/ml)", title = "ER 300")

plot_grid(p1, p2, p3, p4, p5, p6, ncol = 2)
ggsave("plots/suppl/calibration_individuals.pdf", dpi = 300, height = 20, width = 16, units="in")

# Validation -------------------------------------------------------------------

j <- 2502:4926
sum_chains <- length(j)*5

str <- which(pars_name == "Ve_C_central(1)")
end <- which(pars_name == "Va_Kp2c(1)")
parms <- pars_name[str:end]

n = 500 # 500 virtual subject
d <- x[j,,] %>% matrix(nrow = sum_chains) %>% as.data.frame() %>% 
  `colnames<-`(dimnames(x)[[3]])
i <- which(d[,"LnPosterior"]==max(d[,"LnPosterior"])) %>%
  rep(n) # Pick the sample based on the maximum a posteriori estimation

tmp.x <- d[i, parms] 
tmp.x$`Vr_Peff(1)` <- tmp.x$`Vr_Peff(1)` + tmp.x$`Va_Peff(1)`
tmp.x$`Vr_Vmax_met_vitro_liver(1)` <- tmp.x$`Vr_Vmax_met_vitro_liver(1)` + tmp.x$`Va_Vmax_met_vitro_liver(1)`
tmp.x$`Vr_V_central_LKg(1)` <- tmp.x$`Vr_V_central_LKg(1)` + tmp.x$`Va_V_central_LKg(1)` 
tmp.x$`Vr_Kc2p(1)` <- tmp.x$`Vr_Kc2p(1)` + tmp.x$`Va_Kc2p(1)`
tmp.x$`Vr_Kp2c(1)` <- tmp.x$`Vr_Kp2c(1)` + tmp.x$`Va_Kp2c(1)`

tmp.x %>% write.table(file="poppred_validation.dat", row.names=T, sep="\t")

vld <- "./mcsim.bupropion_CAT.model.exe MCSim/bupropion_MCMC_validation.in"
system(vld)

df <- fread("poppred_validation.out") %>% as.data.frame()
file.remove(c("poppred_validation.dat", "poppred_validation.out"))
#dim(df)

# 
str_1.1 <- which(names(df) == "C_central_ngml_1.1")
end_3.2 <- which(names(df) == "C_central_ngml_6.14")
calib <- c("3122", "2666", "3135", "1605", "3146", "3156", "2180",
           "3219", "3422", "3353", "3427", "3239", "3457", "3403", "3474")

BUP <- read_csv("datasets/bupro.csv")
BUP %<>% mutate(group = ifelse(Subject %in% calib, "Calibration", "Validation"))
BUP %<>% unite(FormDose, Formulation, Dose, sep = " ", remove = FALSE)
form <- c("IR 75", "IR 100", "SR 100", "SR 150", "ER 150", "ER 300")
BUP$Formulation <- factor(BUP$FormDose, level = form)

t <- c(0.25,  0.50,  0.75,  1.00,  1.50,  2.00, 3.00, 
       4.00,  6.00, 8.00, 24.00, 48.00, 72.00, 96.00)

qt.line <- df[,c(str_1.1:end_3.2)] %>% gather() %>% 
  `colnames<-`(c("var","Conc")) %>%
  add_column(Time = rep(rep(t, each = n), 6)) %>%
  add_column(Subject = rep(c(1:n), 84)) %>%
  add_column(Formulation = rep(form, each = 14*n)) %>%
  group_by(Formulation, Time) %>%
  summarize(med = median(Conc),
            lcl95 = p.interval(Conc)[1,1],
            ucl95 = p.interval(Conc)[1,2],
            lcl99 = p.interval(Conc, prob = 0.99)[1,1],
            ucl99 = p.interval(Conc, prob = 0.99)[1,2])
qt.line$Formulation <- factor(qt.line$Formulation, level = form)

subset.n <- 50
dat.v <- df[c(1:subset.n),c(str_1.1:end_3.2)] %>% gather() %>% 
  `colnames<-`(c("var","Conc")) %>%
  add_column(Time = rep(rep(t, each = subset.n), 6)) %>%
  add_column(Subject = rep(c(1:subset.n), 84)) %>%
  add_column(Formulation = rep(form, each = 14*subset.n))
dat.v$Formulation <- factor(dat.v$Formulation, level = form)

#tiff("plots/fig7.tiff", res=600, compression = "lzw", height = 8, width = 14, units="in")
pdf("plots/fig7.pdf", height = 8, width = 14)
dat.v %>% ggplot(aes(Time, Conc)) +
  facet_wrap(~Formulation, dir="v", nrow=2)+
  scale_x_log10() +
  scale_color_viridis_d(direction = -1, begin = 0.5, end=0.95) +
  geom_line(data = BUP, aes(color=group, group = Subject))+
  geom_point(data = BUP, aes(color=group, group = Subject)) +
  geom_line(data = qt.line, aes(x=Time, y= med), lwd=1.2, col="grey") +
  geom_line(data = qt.line, aes(x=Time, y= ucl99), lty=3, lwd=0.5) +
  geom_line(data = qt.line, aes(x=Time, y= lcl99), lty=3, lwd=0.5) +
  geom_line(data = qt.line, aes(x=Time, y= ucl95), lty=2, lwd=0.5) +
  geom_line(data = qt.line, aes(x=Time, y= lcl95), lty=2, lwd=0.5) +
  theme_bw() + 
  labs(x="Time, hr", y="Concentration, ng/ml") +
  theme(
    legend.title = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_text(size = 12, colour = "black", face = "bold"),
    axis.title = element_text(size = 12, colour = "black", face = "bold"),
    axis.text = element_text(size = 12, colour = "black", face = "bold"))
dev.off()
#ggsave("plots/fig7.jpeg", dpi = 300, height = 8, width = 14, units="in")



# Estimate Cmax, Tmax, and AUC 
## Calibration group
m_cal <- matrix(nrow = 3, ncol = 6)
for (i in  1:6){
  conc_obj <- BUP %>% filter(Formulation == form[i] & group == "Calibration") %>%
    mutate(conc = ifelse(Time == 0, 0, Conc)) %>% 
    filter(!is.na(conc)) %>% 
    PKNCAconc(conc~Time|Subject)
  d_dose <- BUP %>% filter(Formulation == form[i] & group == "Calibration") %>%
    select(Dose, Time, Subject, Conc) %>% filter(Time == 0)
  dose_obj <- PKNCAdose(d_dose, Dose~Time|Subject)
  data_obj_automatic <- PKNCAdata(conc_obj, dose_obj)
  results_obj_automatic <- pk.nca(data_obj_automatic)
  sum <- summary(results_obj_automatic)
  m_cal[1,i] <- sum$cmax[2] # geometric mean and geometric coefficient of variation 
  m_cal[2,i] <- sum$tmax[2] # median and range
  m_cal[3,i] <- sum$aucinf.obs[2] # geometric mean and geometric coefficient of variation   
}

## Validation group
m_val <- matrix(nrow = 3, ncol = 6)
for (i in  1:6){
  conc_obj <- BUP %>% filter(Formulation == form[i] & group == "Validation") %>%
    mutate(conc = ifelse(Time == 0, 0, Conc)) %>% 
    filter(!is.na(conc)) %>% 
    PKNCAconc(conc~Time|Subject)
  d_dose <- BUP %>% filter(Formulation == form[i] & group == "Validation") %>%
    select(Dose, Time, Subject, Conc) %>% filter(Time == 0)
  dose_obj <- PKNCAdose(d_dose, Dose~Time|Subject)
  data_obj_automatic <- PKNCAdata(conc_obj, dose_obj)
  results_obj_automatic <- pk.nca(data_obj_automatic)
  sum <- summary(results_obj_automatic)
  m_val[1,i] <- sum$cmax[2] # geometric mean and geometric coefficient of variation 
  m_val[2,i] <- sum$tmax[2] # median and range
  m_val[3,i] <- sum$aucinf.obs[2] # geometric mean and geometric coefficient of variation   
}

## Simulation group

str_Sim <- c("C_central_ngml_1.1", "C_central_ngml_2.1",
             "C_central_ngml_3.1", "C_central_ngml_4.1",
             "C_central_ngml_5.1", "C_central_ngml_6.1")
l <- list()
for(i in 1:6){
  str_sim_location <- which(names(df) == str_Sim[i])
  end_sim_location <- str_sim_location + 13
  l[[i]] <- cbind(0, df[,c(str_sim_location:end_sim_location)]) %>% gather()
}

sim_BUP <- do.call(rbind, l) %>%
  `colnames<-`(c("var","Conc")) %>%
  add_column(Time = rep(rep(c(0, t), each = n), 6)) %>%
  add_column(Subject = rep(c(1:n), 90)) %>%
  add_column(Formulation = rep(form, each = 15*n)) %>%
  add_column(Dose = rep(c(75, 100, 100, 150, 150, 300), each = 15*n)) %>%
  group_by(Formulation, Time)

m_sim <- matrix(nrow = 3, ncol = 6)
for (i in  1:6){
  conc_obj <- sim_BUP %>% filter(Formulation == form[i]) %>%
    PKNCAconc(Conc~Time|Subject)
  d_dose <- sim_BUP %>% filter(Formulation == form[i]) %>%
    select(Dose, Time, Subject, Conc) %>% 
    filter(Time == 0)
  
  dose_obj <- PKNCAdose(d_dose, Dose~Time|Subject)
  data_obj_automatic <- PKNCAdata(conc_obj, dose_obj)
  results_obj_automatic <- pk.nca(data_obj_automatic)
  sum <- summary(results_obj_automatic)
  m_sim[1,i] <- sum$cmax[2]  
  m_sim[2,i] <- sum$tmax[2] 
  m_sim[3,i] <- sum$aucinf.obs[2] 
}

m_cal
m_val
m_sim

data.frame(rep(c("Cmax", "Tmax", "AUC"), 6), c(m_cal), c(m_val), c(m_sim)) %>% 
  `colnames<-`(c("Formulato", "Calibation", "Validation", "Simulation")) %>%
  saveRDS(file = "outputs/pk_parameters.rds")

# Virtual individuals (reference) ----------------------------------------------

j <- 2502:4926
d <- x[j,,] %>% matrix(ncol = dim(x)[3]) %>% as.data.frame() %>% 
  `colnames<-`(dimnames(x)[[3]])

nV <- 100
nU <- 100

set.seed(55) 
r.sample <- sort(sample.int(1000, 20)) # take 20 from each chain
v.sample <- rep(r.sample+3000, each=nV) # 25*24

str <- which(pars_name == "Ve_C_central(1)")
end <- which(pars_name == "Vr_Weibull_slope_ER(1)")
parms <- pars_name[str:end]
x[v.sample,,parms] %>%
  matrix(ncol = length(parms)) %>% as.data.frame() %>% 
  `colnames<-`(parms) %>%
  write.table(file="poppred_reference.dat", row.names=T, sep="\t")

ref <- "./mcsim.bupropion_CAT.model.exe MCSim/bupropion_MCMC_reference.in"
system(paste0(ref))
ref <- fread("poppred_reference.out") %>% as.data.frame()
dim(ref) # 10,000 individuals

ref.0 <- which(names(ref) == "C_central_ngml_1.1")
ref.1 <- which(names(ref) == "C_central_ngml_1.14")

ref_CMAX <- apply(ref[,(ref.0:ref.1)], 1, max)
ref_AUCt <- ref$AUC_central_1.1

#median(ref_CMAX); p.interval(ref_CMAX)
#median(ref_AUCt); p.interval(ref_AUCt)

nVU <- nV*nU

ref.dat <- tibble(
  subject = c(1:nVU), #nV*nU
  formulation = "Ref",
  group = rep(paste0("G", c(1:100)), each = nV),
  period = rep(rep(c(1,2), each = nV/2), nU),
  CMAX = ref_CMAX,
  AUCt = ref_AUCt # time = inf
)

#dim(ref.dat)

ref.dat$subject = as.factor(ref.dat$subject)
ref.dat$formulation = as.factor(ref.dat$formulation)
ref.dat$group = as.factor(ref.dat$group)
ref.dat$period = as.factor(ref.dat$period)

file.remove(c("poppred_reference.dat", "poppred_reference.out"))

# Ref vs Test (check the requirement of individuals) ----------------------------
ind.parms <-c("G_Radius",
              "Peff",
              "Fu_plasma",
              "Fu_vitro_liver",
              "Kpuu_liver",
              "Vmax_met_vitro_liver", 
              "V_central_LKg",
              "Kc2p",
              "Kp2c",
              "BDM")

va <- c("Va_Peff(1)",
        "Va_Vmax_met_vitro_liver(1)",
        "Va_V_central_LKg(1)",
        "Va_Kc2p(1)",
        "Va_Kp2c(1)") 

# Use the parameters from virtual individual that generate above
va.parms <- x[v.sample,,va] %>% matrix(ncol = length(va)) %>% as.data.frame() %>% `colnames<-`(va)

# Add Weibull parameters
test.df <- cbind(ref[,ind.parms], va.parms)
map_scale <- d[,"Weibull_scale_IR(1)"] %>% tail(500) %>% map_estimate() %>% c()
test.df$Weibull_scale_IR <- rep(map_scale, each = nV*nU)
test.df$Weibull_slope_IR <- 1

# Use Va minimum 
test.df$`Va_Peff(1)` <- p.interval(d$`Va_Peff(1)`)[1]
test.df$`Va_Vmax_met_vitro_liver(1)` <- p.interval(d$`Va_Vmax_met_vitro_liver(1)`)[1]
test.df$`Va_V_central_LKg(1)` <- p.interval(d$`Va_V_central_LKg(1)`)[1]
test.df$`Va_Kc2p(1)` <- p.interval(d$`Va_Kc2p(1)`)[1]
test.df$`Va_Kp2c(1)` <- p.interval(d$`Va_Kp2c(1)`)[1]

# Simulation
test.df %>% write.table(file="poppred_test.dat", row.names=T, sep="\t")
system("./mcsim.bupropion_CAT.model.exe MCSim/bupropion_MCMC_test.in")

Test <- fread("poppred_test.out") %>% as.data.frame()
file.remove(c("poppred_test.dat", "poppred_test.out"))

# Analysis
Test.0 <- which(names(Test) == "C_central_ngml_1.1")
Test.1 <- which(names(Test) == "C_central_ngml_1.14")

Test_CMAX <- apply(Test[,(Test.0:Test.1)], 1, max)
Test_AUCt <- Test$AUC_central_1.1

test.dat <- tibble(
  subject = c(1:nVU),
  formulation = rep("T1", each=nU*nV),
  group = rep(paste0("G", c(1:100)), each = nV),
  period = rep(rep(c(2,1), each = nV/2), nU),
  CMAX = Test_CMAX,
  AUCt = Test_AUCt # time = 92hr
)

test.dat$subject = as.factor(test.dat$subject)
test.dat$formulation = as.factor(test.dat$formulation)
test.dat$group = as.factor(test.dat$group)
test.dat$period = as.factor(test.dat$period)

T_mean <- test.dat$CMAX %>% log() %>% mean() %>% exp()
R_mean <- ref.dat$CMAX %>% log() %>% mean() %>% exp()
#T_mean / R_mean
T_mean <- test.dat$AUCt %>% log() %>% mean() %>% exp()
R_mean <- ref.dat$AUCt %>% log() %>% mean() %>% exp()
#T_mean / R_mean

#
scale <- map_scale
slope <- 1 
dat_tab.1 <- data.frame(scale, slope)
dat_tab.1$test <- "T1"

t1.list <- list("t1")
for(j in 1:1){
  
  test.form <- dat_tab.1$test[j]
  
  for(i in 1:100){
    
    gp <- unique(test.dat$group)[i]
    
    data <- tibble(
      subject = rep(seq(nV), 2),
      period = c(rep(1, nV/2), rep(2, nV/2), rep(2, nV/2), rep(1, nV/2)),
      formulation = c(rep("R", nV), rep("T", nV)),
      CMAX = c(filter(ref.dat, group == gp)$CMAX, filter(test.dat, group == gp & formulation == test.form)$CMAX) %>% log(),
      AUCt = c(filter(ref.dat, group == gp)$AUCt, filter(test.dat, group == gp & formulation == test.form)$AUCt) %>% log()
    )
    data$subject = as.factor(data$subject)
    data$period = as.factor(data$period)
    data$formulation = as.factor(data$formulation)
    
    mean_CMAX_R <- filter(data, formulation == "R")$CMAX %>% mean() %>% exp() 
    mean_CMAX_T <- filter(data, formulation == "T")$CMAX %>% mean() %>% exp() 
    gmr_CMAX <- mean_CMAX_T / mean_CMAX_R  
    modCMAX <- lm(CMAX ~ subject + period + formulation, data = data)
    mseCMAX <- anova(modCMAX)[[3]][[4]]
    cvCMAX <- mse2CV(mseCMAX)
    ci_CMAX <- CI.BE(pe = gmr_CMAX, CV = cvCMAX, n = nV)
    
    mean_AUCt_R <- filter(data, formulation == "R")$AUCt %>% mean() %>% exp() 
    mean_AUCt_T <- filter(data, formulation == "T")$AUCt %>% mean() %>% exp() 
    gmr_AUCt <- mean_AUCt_T / mean_AUCt_R  
    modAUCt <- lm(AUCt ~ subject + period + formulation, data = data)
    mseAUCt <- anova(modAUCt)[[3]][[4]]
    cvAUCt <- mse2CV(mseAUCt)
    ci_AUCt <- CI.BE(pe = gmr_AUCt, CV = cvAUCt, n = nV)
    
    if(i == 1) {
      GMR_CMAX <- c(gmr_CMAX, cvCMAX, ci_CMAX)
      GMR_AUCt <- c(gmr_AUCt, cvAUCt, ci_AUCt)
    } else {
      GMR_CMAX <- c(GMR_CMAX, c(gmr_CMAX, cvCMAX, ci_CMAX))
      GMR_AUCt <- c(GMR_AUCt, c(gmr_AUCt, cvAUCt, ci_AUCt))
    }
  }
  test_CMAX <- matrix(GMR_CMAX, nco=4, byrow = T) %>% as.data.frame() %>%
    mutate(T1 = ifelse(V1>1.25 | V1<0.8, 1, 0)) %>%
    mutate(T2 = ifelse(V3>1.25 | V3<0.8, 1, 0)) %>%
    mutate(T3 = ifelse(V4>1.25 | V4<0.8, 1, 0))
  test_AUCt <- matrix(GMR_AUCt, nco=4, byrow = T) %>% as.data.frame() %>%
    mutate(T4 = ifelse(V1>1.25 | V1<0.8, 1, 0)) %>%
    mutate(T5 = ifelse(V3>1.25 | V3<0.8, 1, 0)) %>%
    mutate(T6 = ifelse(V4>1.25 | V4<0.8, 1, 0))
  t1.list[[j]] <- cbind(test_CMAX, test_AUCt) %>%
    `colnames<-`(c("V1","V2","V3","V4", # GMR, CV, LCL, UCL
                   "T1","T2","T3",
                   "V5","V6","V7","V8",
                   "T4","T5","T6")) %>%
    mutate(fail = ifelse(T1+T2+T3+T4+T5+T6==0, 0, 1))
  
}

for(j in 1:1){
  dat_tab.1$max.cv[j] <- max(c(t1.list[[j]]$V2, t1.list[[j]]$V6))
  dat_tab.1$pass[j] <- (nU - length(which(t1.list[[j]]$fail==1))) * 0.01
}
dat_tab.1
sampleN.TOST(CV=max(dat_tab.1$max.cv)) # required xx individuals

# First-order ------------------------------------------------------------------

# Set the range of Weibull scale
map_scale <- d[, "Weibull_scale_IR(1)"] %>% map_estimate() %>% c()

Weibull_scale <- seq(log(map_scale/10),
                     log(map_scale*10),
                     (log(map_scale*10)-log(map_scale))/10) %>% exp()
#Weibull_scale
Weibull_slope <- 1

temp.df <- cbind(ref[,ind.parms], va.parms)

test.df <- do.call(rbind,
                   list(temp.df, temp.df, temp.df, temp.df,
                        temp.df, temp.df, temp.df, temp.df,
                        temp.df, temp.df, temp.df, temp.df,
                        temp.df, temp.df, temp.df, temp.df,
                        temp.df, temp.df, temp.df, temp.df, temp.df))
#dim(test.df)

n.grid <- 21
test.df$Weibull_scale_IR <- rep(Weibull_scale, each = nV*nU)
test.df$Weibull_slope_IR <- 1
test.df %>% write.table(file="poppred_test.dat", row.names=T, sep="\t")

#Str.t <- Sys.time()
system("./mcsim.bupropion_CAT.model.exe MCSim/bupropion_MCMC_test.in")
#Sys.time() - Str.t
# 48*100*21: Time difference of 4.652673 mins

Test_1 <- fread("poppred_test.out") %>% as.data.frame()
file.remove(c("poppred_test.dat", "poppred_test.out"))
#dim(Test_1)

Test_1.0 <- which(names(Test_1) == "C_central_ngml_1.1")
Test_1.1 <- which(names(Test_1) == "C_central_ngml_1.14")

Test_1_CMAX <- apply(Test_1[,(Test_1.0:Test_1.1)], 1, max)
Test_1_AUCt <- Test_1$AUC_central_1.1

test.dat <- tibble(
  subject = rep(c(1:nVU), 21),
  formulation = rep(paste0("T", c(1:21)), each=nU*nV),
  group = rep(rep(paste0("G", c(1:100)), each = nV), 21),
  period = rep(rep(rep(c(2,1), each = nV/2), nU), 21),
  CMAX = Test_1_CMAX,
  AUCt = Test_1_AUCt # time = 92hr
)

test.dat$subject = as.factor(test.dat$subject)
test.dat$formulation = as.factor(test.dat$formulation)
test.dat$group = as.factor(test.dat$group)
test.dat$period = as.factor(test.dat$period)

#test.dat$CMAX %>% log() %>% density() %>% plot()
#ref.dat$CMAX %>% log() %>% density() %>% lines()

#test.dat$AUCt %>% log() %>% density() %>% plot()
#ref.dat$AUCt %>% log() %>% density() %>% lines()

#T_mean <- test.dat$CMAX %>% log() %>% mean() %>% exp()
#R_mean <- ref.dat$CMAX %>% log() %>% mean() %>% exp()
#T_mean / R_mean
#T_mean <- test.dat$AUCt %>% log() %>% mean() %>% exp()
#R_mean <- ref.dat$AUCt %>% log() %>% mean() %>% exp()
#T_mean / R_mean

#
scale <- rep(round(Weibull_scale, 4))
slope <- 1 
dat_tab.1 <- data.frame(scale, slope)
dat_tab.1$test <- paste0("T", c(1:21))

t1.list <- list("t1")
#Str.t <- Sys.time()
for(j in 1:21){
  
  test.form <- dat_tab.1$test[j]
  
  for(i in 1:100){
    gp <- unique(test.dat$group)[i]
    
    data <- tibble(
      subject = rep(seq(nV), 2),
      period = c(rep(1, nV/2), rep(2, nV/2), rep(2, nV/2), rep(1, nV/2)),
      formulation = c(rep("R", nV), rep("T", nV)),
      CMAX = c(filter(ref.dat, group == gp)$CMAX, filter(test.dat, group == gp & formulation == test.form)$CMAX) %>% log(),
      AUCt = c(filter(ref.dat, group == gp)$AUCt, filter(test.dat, group == gp & formulation == test.form)$AUCt) %>% log()
    )
    data$subject = as.factor(data$subject)
    data$period = as.factor(data$period)
    data$formulation = as.factor(data$formulation)
    
    mean_CMAX_R <- filter(data, formulation == "R")$CMAX %>% mean() %>% exp() 
    mean_CMAX_T <- filter(data, formulation == "T")$CMAX %>% mean() %>% exp() 
    gmr_CMAX <- mean_CMAX_T / mean_CMAX_R  
    modCMAX <- lm(CMAX ~ subject + period + formulation, data = data)
    mseCMAX <- anova(modCMAX)[[3]][[4]]
    cvCMAX <- mse2CV(mseCMAX)
    ci_CMAX <- CI.BE(pe = gmr_CMAX, CV = cvCMAX, n = nV)
    
    mean_AUCt_R <- filter(data, formulation == "R")$AUCt %>% mean() %>% exp() 
    mean_AUCt_T <- filter(data, formulation == "T")$AUCt %>% mean() %>% exp() 
    gmr_AUCt <- mean_AUCt_T / mean_AUCt_R  
    modAUCt <- lm(AUCt ~ subject + period + formulation, data = data)
    mseAUCt <- anova(modAUCt)[[3]][[4]]
    cvAUCt <- mse2CV(mseAUCt)
    ci_AUCt <- CI.BE(pe = gmr_AUCt, CV = cvAUCt, n = nV)
    
    if(i == 1) {
      GMR_CMAX <- c(gmr_CMAX, cvCMAX, ci_CMAX)
      GMR_AUCt <- c(gmr_AUCt, cvAUCt, ci_AUCt)
    } else {
      GMR_CMAX <- c(GMR_CMAX, c(gmr_CMAX, cvCMAX, ci_CMAX))
      GMR_AUCt <- c(GMR_AUCt, c(gmr_AUCt, cvAUCt, ci_AUCt))
    }
  }
  test_CMAX <- matrix(GMR_CMAX, nco=4, byrow = T) %>% as.data.frame() %>%
    mutate(T1 = ifelse(V1>1.25 | V1<0.8, 1, 0)) %>%
    mutate(T2 = ifelse(V3>1.25 | V3<0.8, 1, 0)) %>%
    mutate(T3 = ifelse(V4>1.25 | V4<0.8, 1, 0))
  test_AUCt <- matrix(GMR_AUCt, nco=4, byrow = T) %>% as.data.frame() %>%
    mutate(T4 = ifelse(V1>1.25 | V1<0.8, 1, 0)) %>%
    mutate(T5 = ifelse(V3>1.25 | V3<0.8, 1, 0)) %>%
    mutate(T6 = ifelse(V4>1.25 | V4<0.8, 1, 0))
  t1.list[[j]] <- cbind(test_CMAX, test_AUCt) %>%
    `colnames<-`(c("V1","V2","V3","V4",
                   "T1","T2","T3",
                   "V5","V6","V7","V8",
                   "T4","T5","T6")) %>%
    mutate(fail = ifelse(T1+T2+T3+T4+T5+T6==0, 0, 1))
  
}
#Sys.time() - Str.t
# 48*100*21: Time difference of 1.658254 mins

# 
# save(t1.list, file = "outputs/t1.list.rda")
# load(file = "outputs/t1.list.rda")

for(j in 1:21){
  dat_tab.1$max.cv[j] <- max(c(t1.list[[j]]$V2, t1.list[[j]]$V6))
  dat_tab.1$pass[j] <- (nU - length(which(t1.list[[j]]$fail==1))) * 0.01
}
sampleN.TOST(CV=max(dat_tab.1$max.cv)) # required xx individuals

BE.lab <- c("0-20", "20-40", "40-60", "60-80", "80-100")
dat_tab.1 %<>% 
  mutate(level.pass = cut(pass, breaks=c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf), labels = BE.lab))
dat_tab.1$level.pass <- factor(dat_tab.1$level.pass, level = rev(BE.lab))

dat_tab.1
# save(dat_tab.1, file = "outputs/dat_tab.1.rda")
# load(file = "outputs/dat_tab.1.rda")

hdi <- d[,"Weibull_scale_IR(1)"] %>% hdi() 

p8.1 <- dat_tab.1 %>% 
  mutate(pass.text = ifelse(pass > 0, paste0(pass*100, "%"), NA)) %>%
  ggplot(aes(scale, slope)) + 
  geom_tile(aes(fill = level.pass), colour = "white") +
  geom_text(aes(label=pass.text), size=5, angle = 90) +
  scale_fill_viridis_d(option = "plasma", begin = 0.2, end = 0.95, direction = -1) +
  labs(x = "Scale", y="", title="First-order release rate (slope = 1)") +
  scale_x_log10()+
  geom_point(aes(x=map_scale, y=0.5), shape = 17, cex =8)+
  #geom_point(aes(x=hdi$CI_low, y=0.5), shape = "l", cex =4)+
  #geom_point(aes(x=hdi$CI_high, y=0.5), shape = "l", cex =4)+
  geom_point(aes(x=map_scale*10-1.2, y=1.5), shape = "l", cex =8)+
  geom_point(aes(x=map_scale*10, y=1.5), shape = "l", cex =8)+
  theme_minimal() + 
  theme(legend.position = "none", 
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(size = 10, colour = "black", face = "bold"),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 12, colour = "black", face = "bold"))

# Weibull -------------------------------------------------------------------

# Set group 
temp.df <- cbind(ref[,ind.parms], va.parms)

test.df <- do.call(rbind, list(temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df,
                               temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df,
                               temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df,
                               temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df,
                               temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df,
                               temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df,
                               temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df,
                               temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df,
                               temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df,
                               temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df,
                               temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df, temp.df))

Weibull_scale <- seq(map_scale*10-1.2, map_scale*10, 0.12)
Weibull_slope <- seq(1.3, 2.3, 0.1)
n.grid <- length(Weibull_scale)

test.df$Weibull_scale_IR <- rep(Weibull_scale, each = nV*nU*n.grid)
test.df$Weibull_slope_IR <- rep(rep(Weibull_slope, each = nV*nU), n.grid)
#dim(test.df)
test.df$G_Radius <- round(test.df$G_Radius, 0)

parr.1 <- c(1:c(nrow(test.df)/4))
parr.2 <- c(c(nrow(test.df)/4+1) : c(nrow(test.df)/2))
parr.3 <- c(c(nrow(test.df)/4*2+1) : c(nrow(test.df)/4*3))
parr.4 <- c(c(nrow(test.df)/4*3+1) : nrow(test.df))
test.df[parr.1,] %>% write.table(file="poppred_test_1.dat", row.names=T, sep="\t")
test.df[parr.2,] %>% write.table(file="poppred_test_2.dat", row.names=T, sep="\t")
test.df[parr.3,] %>% write.table(file="poppred_test_3.dat", row.names=T, sep="\t")
test.df[parr.4,] %>% write.table(file="poppred_test_4.dat", row.names=T, sep="\t")

detectCores()
cores <- 4 
cl <- makeCluster(cores)
registerDoParallel(cl)

# Parallel computing
#strt<-Sys.time()
system.time( 
  out <- foreach(i = 1:cores) %dopar% { 
    tx  <- readLines("MCSim/bupropion_MCMC_test.in")
    tx2 <- gsub(pattern = "poppred_test.dat", replace = paste0("poppred_test_", i,".dat"), x = tx)
    input <- paste0("MCSim/bupropion_MCMC_test_", i,".in")
    writeLines(tx2, con = input)
    exec <- paste0("./mcsim.bupropion_CAT.model.exe MCSim/bupropion_MCMC_test_", i,".in", " poppred_test_", i ,".out")
    system(exec)
    file.remove(input)
  }
)
#print(Sys.time()-strt)
stopCluster(cl)   
# 48*100*11*11: Time difference of 31.60807 mins

Test_2 <- do.call(rbind, list(
  fread("poppred_test_1.out") %>% as.data.frame(),
  fread("poppred_test_2.out") %>% as.data.frame(),
  fread("poppred_test_3.out") %>% as.data.frame(),
  fread("poppred_test_4.out") %>% as.data.frame()
))
#dim(Test_2) # check computing error
#file.remove(c(paste0("poppred_test_", 1:4, ".dat"), paste0("poppred_test_", 1:4, ".out")))

Test_2.0 <- which(names(Test_2) == "C_central_ngml_1.1")
Test_2.1 <- which(names(Test_2) == "C_central_ngml_1.14")

Test_2_CMAX <- apply(Test_2[,(Test_2.0:Test_2.1)], 1, max)
Test_2_AUCt <- Test_2$AUC_central_1.1

test.dat <- tibble(
  subject = rep(c(1:nVU), n.grid^2),
  formulation = rep(paste0("T", c(1:n.grid^2)), each=nU*nV),
  group = rep(rep(paste0("G", c(1:100)), each = nV), n.grid^2),
  period = rep(rep(rep(c(2,1), each = nV/2), nU), n.grid^2),
  CMAX = Test_2_CMAX,
  AUCt = Test_2_AUCt # time = 92hr
)

test.dat$subject = as.factor(test.dat$subject)
test.dat$formulation = as.factor(test.dat$formulation)
test.dat$group = as.factor(test.dat$group)
test.dat$period = as.factor(test.dat$period)

#test.dat$CMAX %>% log() %>% density() %>% plot()
#ref.dat$CMAX %>% log() %>% density() %>% lines()
#test.dat$AUCt %>% log() %>% density() %>% plot()
#ref.dat$AUCt %>% log() %>% density() %>% lines()

#T_mean <- test.dat$CMAX %>% log() %>% mean() %>% exp()
#R_mean <- ref.dat$CMAX %>% log() %>% mean() %>% exp()
#T_mean / R_mean
#T_mean <- test.dat$AUCt %>% log() %>% mean() %>% exp()
#R_mean <- ref.dat$AUCt %>% log() %>% mean() %>% exp()
#T_mean / R_mean

#
scale <- rep(round(Weibull_scale, 2), each=n.grid)
slope <- rep(round(Weibull_slope, 2), n.grid)
dat_tab.2 <- data.frame(scale, slope)
dat_tab.2$test <- paste0("T", c(1:n.grid^2))

t2.list <- list("t2.list")
#Str.t <- Sys.time()
for(j in 1:n.grid^2){
  test.form <- dat_tab.2$test[j]
  for(i in 1:100){
    gp <- unique(test.dat$group)[i]
    data <- tibble(
      subject = rep(seq(nV), 2),
      period = c(rep(1, nV/2), rep(2, nV/2), rep(2, nV/2),
                 rep(1, nV/2)),
      formulation = c(rep("R", nV), rep("T", nV)),
      CMAX = c(filter(ref.dat, group == gp)$CMAX, filter(test.dat, group == gp & formulation == test.form)$CMAX) %>% log(),
      AUCt = c(filter(ref.dat, group == gp)$AUCt, filter(test.dat, group == gp & formulation == test.form)$AUCt) %>% log())
    data$subject = as.factor(data$subject)
    data$period = as.factor(data$period)
    data$formulation = as.factor(data$formulation)
    
    mean_CMAX_R <- filter(data, formulation == "R")$CMAX %>%
      mean() %>% exp() 
    mean_CMAX_T <- filter(data, formulation == "T")$CMAX %>%
      mean() %>% exp() 
    gmr_CMAX <- mean_CMAX_T / mean_CMAX_R  
    modCMAX <- lm(CMAX ~ subject + period + formulation, data = data)
    mseCMAX <- anova(modCMAX)[[3]][[4]]
    cvCMAX <- mse2CV(mseCMAX)
    ci_CMAX <- CI.BE(pe = gmr_CMAX, CV = cvCMAX, n = nV)
    
    mean_AUCt_R <- filter(data, formulation == "R")$AUCt %>%
      mean() %>% exp() 
    mean_AUCt_T <- filter(data, formulation == "T")$AUCt %>%
      mean() %>% exp() 
    gmr_AUCt <- mean_AUCt_T / mean_AUCt_R  
    modAUCt <- lm(AUCt ~ subject + period + formulation, data = data)
    mseAUCt <- anova(modAUCt)[[3]][[4]]
    cvAUCt <- mse2CV(mseAUCt)
    ci_AUCt <- CI.BE(pe = gmr_AUCt, CV = cvAUCt, n = nV)
    
    if(i == 1) {
      GMR_CMAX <- c(gmr_CMAX, cvCMAX, ci_CMAX)
      GMR_AUCt <- c(gmr_AUCt, cvAUCt, ci_AUCt)
    } else {
      GMR_CMAX <- c(GMR_CMAX, c(gmr_CMAX, cvCMAX, ci_CMAX))
      GMR_AUCt <- c(GMR_AUCt, c(gmr_AUCt, cvAUCt, ci_AUCt))
    }
  }
  test_CMAX <- matrix(GMR_CMAX, nco=4, byrow = T) %>%
    as.data.frame() %>%
    mutate(T1 = ifelse(V1>1.25 | V1<0.8, 1, 0)) %>%
    mutate(T2 = ifelse(V3>1.25 | V3<0.8, 1, 0)) %>%
    mutate(T3 = ifelse(V4>1.25 | V4<0.8, 1, 0))
  test_AUCt <- matrix(GMR_AUCt, nco=4, byrow = T) %>%
    as.data.frame() %>%
    mutate(T4 = ifelse(V1>1.25 | V1<0.8, 1, 0)) %>%
    mutate(T5 = ifelse(V3>1.25 | V3<0.8, 1, 0)) %>%
    mutate(T6 = ifelse(V4>1.25 | V4<0.8, 1, 0))
  t2.list[[j]] <- cbind(test_CMAX, test_AUCt) %>%
    `colnames<-`(c("V1","V2","V3","V4",
                   "T1","T2","T3",
                   "V5","V6","V7","V8",
                   "T4","T5","T6")) %>%
    mutate(fail = ifelse(T1+T2+T3+T4+T5+T6==0, 0, 1))
}
#Sys.time() - Str.t
# 100*100*11*11: Time difference of 25.83942 mins

# save(t2.list, file = "outputs/t2.list.rda")
# load(file = "outputs/t2.list.rda")

for(j in 1:n.grid^2){
  dat_tab.2$max.cv[j] <- max(c(t2.list[[j]]$V2, t2.list[[j]]$V6))
  dat_tab.2$pass[j] <- (nU - length(which(t2.list[[j]]$fail==1))) * 0.01
}
sampleN.TOST(CV=max(dat_tab.2$max.cv)) # required 44 individuals

BE.lab <- c("0-20", "20-40", "40-60", "60-80", "80-100")
dat_tab.2 %<>% 
  mutate(level.pass = cut(pass, breaks=c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf), labels = BE.lab))
dat_tab.2$level.pass <- factor(dat_tab.2$level.pass, level = rev(BE.lab))

dat_tab.2
# save(dat_tab.2, file = "outputs/dat_tab.2.rda")
# load(file = "outputs/dat_tab.2.rda")

p8.2 <- dat_tab.2 %>% 
  mutate(pass.text = ifelse(pass > 0, paste0(pass*100, "%"), NA)) %>% 
  ggplot(aes(scale, slope)) + 
  geom_tile(aes(fill = level.pass), colour = "white") +
  geom_text(aes(label=pass.text)) +
  scale_fill_viridis_d(option = "plasma", begin = 0.2, end = 0.95, direction = -1) +
  labs(x = "Scale", y="Slope", title="Weibull model release rate") +
  guides(fill=guide_legend("Success rate (%)")) +
  theme_minimal() + 
  geom_point(aes(x=map_scale*10-1.2, y=1.2), shape = "l", cex =8)+
  geom_point(aes(x=map_scale*10, y=1.2), shape = "l", cex =8)+
  theme(axis.ticks = element_blank(),
        panel.border = element_blank(),
        axis.text = element_text(size = 10, colour = "black", face = "bold"),
        axis.title = element_text(size = 12, colour = "black", face = "bold"))

#tiff("plots/fig8.tiff", res=600, compression = "lzw", height = 9, width = 12, units="in")
pdf("plots/fig8.pdf", height = 9, width = 12)
plot_grid(p8.2, p8.1, nrow = 2, rel_heights = c(4/5, 1/5),  align = "v", axis = "l")
dev.off()
#ggsave("plots/fig8.jpeg", dpi = 300, height = 9, width = 12, units="in")

