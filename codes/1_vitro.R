# Load packages ----------------------------------------------------------------
library(foreach)
library(doParallel)
library(data.table) 
library(tidyverse)  
library(rstan)
library(reshape2)
library(bayesplot)
library(bayestestR)
library(LaplacesDemon)
library(cowplot)
source("MCSim/function.R")

# Run MCMC (can remove comment to re-run)---------------------------------------
#model <- "bupropion_CAT.model" 

#current.files <- list.files()
#detectCores()
#cores <- 4    # 4 chains 
#cl <- makeCluster(cores)
#registerDoParallel(cl)

#strt<-Sys.time()
#out <- foreach(i = 1:cores) %dopar% 
#  {
#    set.seed(i+10)
#    mcsim(model = model, input = "bupropion_MCMC_vitro.in", parallel = T)  
#  }
#print(Sys.time()-strt) # Time difference of 25 secs

#new.files <- setdiff(list.files(), current.files)
#to.remove <- new.files[grep('.kernel|.in', new.files)]
#file.remove(to.remove)
#out.files <- setdiff(list.files(), current.files)

#for(i in 1:4){
#  file.copy(out.files[i], paste0("outputs/", out.files[i]))
#  file.remove(out.files[i])
#}

# Posterior check ---------------------------------------------------------
sim1 <- fread("outputs/bupropion_MCMC_vitro_1375.out") %>% as.data.frame()
sim2 <- fread("outputs/bupropion_MCMC_vitro_2896.out") %>% as.data.frame() 
sim3 <- fread("outputs/bupropion_MCMC_vitro_5656.out") %>% as.data.frame() 
sim4 <- fread("outputs/bupropion_MCMC_vitro_7099.out") %>% as.data.frame() 
x <- mcmc_array(list(sim1, sim2,sim3, sim4))
pars_name <- dimnames(x)[[3]]
weibull_param <- pars_name[2:9] 
mnt <- monitor(x[,,weibull_param], digit=4)

# Create mcmc report sheet
Parameter <- c("Weibull_scale_IR75", "Weibull_scale_IR100", 
               "Weibull_scale_SR100", "Weibull_scale_SR150",
               "Weibull_scale_ER150", "Weibull_slope_ER150",
               "Weibull_scale_ER300", "Weibull_slope_ER300")
Q5 <- mnt$Q5 %>% round(digits = 4)
Q50 <- mnt$Q50 %>% round(digits = 4)
Q95 <- mnt$Q95 %>% round(digits = 4)
Mean <- mnt$mean %>% round(digits = 4)
SD <- mnt$sd %>% round(digits = 4)
n_eff <- mnt$n_eff
Rhat <- mnt$Rhat %>% round(digits = 4)
data.frame(Parameter, Q5, Q50, Q95, Mean, SD, n_eff, Rhat) %>%
  write.csv(file="outputs/mcmc_report_1.csv", row.names=FALSE)



nrow(x)
j=502:1001
#mcmc_dens_overlay(x[j,,], pars = weibull_param, facet_args = list(ncol = 4, strip.position = "top")) 
#mcmc_trace(x[j,,], pars = weibull_param, facet_args = list(ncol = 4, strip.position = "top"))

weibull_scale_IR <- c("Weibull_scale_IR(1.1)", "Weibull_scale_IR(1.2)")
weibull_scale_SR <- c("Weibull_scale_SR(2.1)", "Weibull_scale_SR(2.2)")
weibull_scale_ER <- c("Weibull_scale_ER(3.1)", "Weibull_scale_ER(3.2)")
weibull_slope_ER <- c("Weibull_slope_ER(3.1)", "Weibull_slope_ER(3.2)")
#hdi(x[j,,weibull_scale_IR]) 
#hdi(x[j,,weibull_scale_SR]) 
#hdi(x[j,,c("Weibull_scale_ER(3.1)", "Weibull_scale_ER(3.2)")]) 
#hdi(x[j,,c("Weibull_slope_ER(3.1)", "Weibull_slope_ER(3.2)")]) 

rel.df <- x[j,,weibull_param] %>% 
  matrix(ncol=4) %>% melt() %>%
  add_column(param = c(rep("Scale", 10000),
                       rep("Slope", 2000),
                       rep("Scale", 2000),
                       rep("Slope", 2000))) %>%
  add_column(form = c(rep("IR", 4000),
                      rep("SR", 4000),
                      rep("ER", 8000)))
rel.df$form <- factor(rel.df$form, levels = c("IR", "SR", "ER"))

# IR dissolution profile -------------------------------------------------------
t <- seq(0, 1, 0.05)
ir.out <- x[j,,weibull_scale_IR] %>% melt() %>% as.data.frame()
sr.out <- x[j,,weibull_scale_SR] %>% melt() %>% as.data.frame()
for(i in 1:length(t)){
  ir.out[,4+i] <- 1 - exp(-(t[i]/ir.out[,4]))
}
quant <- c(0.5, 0.025, 0.975)
IR.pred <- apply(ir.out[,-c(1:4)], 2, quantile, quant) %>% 
  t() %>% as.data.frame()
for(i in 6:25){
  IR.pred[i-4,c(2,3)] <- p.interval(ir.out[,i])[1,] 
}
IR.pred$time <- t

# SR dissolution profile -------------------------------------------------------
t <- seq(0, 10, 0.5)
for(i in 1:length(t)){
  sr.out[,4+i] <- 1 - exp(-(t[i]/sr.out[,4]))
}
SR.pred <- apply(sr.out[,-c(1:4)], 2, quantile, quant) %>% 
  t() %>% as.data.frame()
for(i in 6:25){
  SR.pred[i-4,c(2,3)] <- p.interval(sr.out[,i])[1,] 
}
SR.pred$time <- t

# ER dissolution profile -------------------------------------------------------
t <- seq(0, 10, 0.5)
er.out <- x[j,,weibull_scale_ER] %>% melt() %>% as.data.frame()
er.out[,5] <- x[j,,weibull_slope_ER] %>% melt() %>% select(value)
for(i in 1:length(t)){
  er.out[,5+i] <- 1 - exp(-(t[i]/er.out[,4])^er.out[,5])
}
ER.pred <- apply(er.out[,-c(1:5)], 2, quantile, quant) %>%
  t() %>% as.data.frame()
for(i in 7:26){
  ER.pred[i-5, c(2,3)] <- p.interval(er.out[,i])[1,] 
}
ER.pred$time <- t

# IR + SR + ER -----------------------------------------------------------------
rel.df2 <- do.call(rbind, list(IR.pred, SR.pred, ER.pred)) %>% 
  add_column(Form = rep(c("IR", "SR", "ER"), each=21)) %>%
  `colnames<-`(c("med", "lcl", "ucl", "time", "form")) 
rel.df2$form <- factor(rel.df2$form, levels = c("IR", "SR", "ER"))

# Dissolution data -------------------------------------------------------------
diss.t <- c(rep(c(0.083, 0.167, 0.25, 0.50, 0.75), 2),
            rep(c(1,2,4,8), 4))
diss.y <- c(0.34, 0.75, 0.85, 0.91, 0.95, 0.32, 0.72, 0.88, 0.94, 0.96,
            0.32, 0.50, 0.73, 0.98, 0.12, 0.37, 0.72, 0.91,
            0.15, 0.43, 0.83, 0.98, 0.12, 0.37, 0.72, 0.91)
vitro.dat <- data.frame(diss.t, diss.y)
vitro.dat$form <- c(rep("IR", 10), rep("SR", 8), rep("ER", 8))
vitro.dat$form <- factor(vitro.dat$form, levels = c("IR", "SR", "ER"))

# Plot -------------------------------------------------------------------------
p4.1 <- rel.df %>% ggplot(aes(x=form, y=value), group=param) +
  geom_violin(aes(fill = param)) +
  facet_wrap(~form, ncol = 3, scales = "free") +
  scale_fill_grey(start=0.6, end=0.9) +
  theme_bw() +
  labs(y = "Weibull parameter") +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(colour = "black", size = 10, face = "bold"),
        legend.title = element_blank(),
        legend.position = c(0.06, 0.12), 
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()
  )

p4.2 <- rel.df2 %>%
  ggplot() +
  facet_wrap(~form, scales = "free") +
  geom_line(aes(x=time, y=med)) +
  geom_ribbon(aes(x = time, ymin=lcl, ymax=ucl), alpha = 0.3) +
  geom_point(data = vitro.dat, aes(x=diss.t, y=diss.y)) + 
  labs(x= "Time, hr", y="%Dissolution") +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.text = element_blank())

tiff("plots/fig4.tiff", res=600, compression = "lzw", height = 6, width = 9, units="in")
plot_grid(p4.1, p4.2, nrow = 2, rel_heights = c(2/3, 1/3),
          align = "v", axis = "l", labels = c("A", "B"))
dev.off()
#ggsave("plots/fig4.jpeg", dpi = 300, height = 6, width = 9, units="in")
