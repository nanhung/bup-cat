library(tidyverse)
library(data.table)
library(cowplot)
source("MCSim/function.R")

# IVIVC -------------------------------------------------------------------

CR_1 <- c(0.05, 1.0)
CR_2 <- c(0.5, 1.0)
CR_3 <- c(1.0, 2.3)

test.name <- c("CR-1 (scale=0.05, slope=1.0)", 
               "CR-2 (scale=0.5, slope=1.0)",
               "CR-3 (scale=1.0, slope=2.3)")

diss.t <- rep(seq(0, 2.5, 0.05), 3)
scale <- rep(c(CR_1[1], CR_2[1], CR_3[1]), each=51)
slope <- rep(c(CR_1[2], CR_2[2], CR_3[2]), each=51)
test <- rep(test.name, each=51)

data.t<-c(0.083, 0.167, 0.25, 0.50, 0.75)
data.y<-c(0.34, 0.75, 0.85, 0.91, 0.95, 0.32, 0.72, 0.88, 0.94, 0.96)*100
diss.dat <- data.frame(data.t, data.y)

diss.pred <- data.frame(diss.t, scale, slope, test)
diss.pred$test <- factor(diss.pred$test, levels = test.name)

safe.space <- data.frame(diss.t, 0.05, 1.0, 2.3) %>%
  `colnames<-`(c("time", "ir.scale", "cr.scale", "cr.slope")) %>%
  mutate(y.max = 100*(1 - exp(-(time/ir.scale)))) %>%
  mutate(y.min = 100*(1 - exp(-(time/cr.scale)^cr.slope)))

p9.1 <- diss.pred %>%
  mutate(diss.y = 100*(1 - exp(-(diss.t/scale)^slope)) ) %>%
  ggplot() +
  geom_hline(yintercept = 80, lty=4) +
  geom_vline(xintercept = 1.24, lty=4) +
  geom_point(data = diss.dat, aes(x=data.t, y=data.y)) +
  geom_ribbon(data = safe.space, aes(x=time, ymin=y.min , ymax=y.max), alpha = .1) +
  annotate(geom="text", x=0.5, y=50, label="Safe space", color="grey50", size = 9, angle=40, fontface="bold") +
  geom_line(aes(x=diss.t, y=diss.y, lty=test, col=test)) +
  scale_color_manual(values=c("grey20", "grey50", "grey80")) +
  theme_bw() +
  labs(x="Time, hr", y="%Dissolved") +
  theme(legend.title = element_blank(),
        legend.position = "top",
        axis.text = element_text(size = 10, colour = "black", face = "bold"),
        axis.title = element_text(size = 12, colour = "black", face = "bold"))

#
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

j <- 2502:4926
sum_chains <- length(j)*5
d <- x[j,,] %>% matrix(nrow = sum_chains) %>% as.data.frame() %>% 
  `colnames<-`(dimnames(x)[[3]])

M.end <-which(names(d)=="Kp2c(1)")
Vr.str <- which(names(d)=="Vr_G_Radius(1)")
Vr.end <- which(names(d)=="Vr_Kp2c(1)")
d$`Vr_Peff(1)` <- d$`Vr_Peff(1)` + d$`Va_Peff(1)`
d$`Vr_Vmax_met_vitro_liver(1)` <- d$`Vr_Vmax_met_vitro_liver(1)` + d$`Va_Vmax_met_vitro_liver(1)`
d$`Vr_V_central_LKg(1)` <- d$`Vr_V_central_LKg(1)` + d$`Va_V_central_LKg(1)` 
d$`Vr_Kc2p(1)` <- d$`Vr_Kc2p(1)` + d$`Va_Kc2p(1)`
d$`Vr_Kp2c(1)` <- d$`Vr_Kp2c(1)` + d$`Va_Kp2c(1)`  

nMC <- 900
d[rep(which(d$LnPosterior == max(d$LnPosterior)), nMC), c(1:M.end, Vr.str:Vr.end)] %>% 
  write.table(file="poppred_be.dat", row.names=F, sep="\t")
system(paste0("./mcsim.bupropion_CAT.model.exe MCSim/bupropion_MCMC_be.in"))

comp <- fread("poppred_be.out") %>% as.data.frame()

comp_1.0 <- which(names(comp) == "C_central_ngml_1.1")
comp_1.1 <- which(names(comp) == "C_central_ngml_1.12")
comp_2.0 <- which(names(comp) == "C_central_ngml_2.1")
comp_2.1 <- which(names(comp) == "C_central_ngml_2.12")
comp_3.0 <- which(names(comp) == "C_central_ngml_3.1")
comp_3.1 <- which(names(comp) == "C_central_ngml_3.12")

test <- rep(test.name, each=nMC) %>% factor(levels = test.name)
CMAX <- c(apply(comp[,(comp_1.0:comp_1.1)], 1, max),
          apply(comp[,(comp_2.0:comp_2.1)], 1, max),
          apply(comp[,(comp_3.0:comp_3.1)], 1, max))
AUCt <- c(comp$AUC_central_1.1, comp$AUC_central_2.1,
          comp$AUC_central_3.1)


p9.2 <- data.frame(rep(test, 2), c(CMAX, AUCt)) %>%
  add_column(parameter = rep(c("Cmax", "AUC"), each = nMC*3)) %>% 
  `colnames<-`(c("Variable","value","parameter")) %>%
  ggplot(aes(x=value)) + 
  facet_wrap(~parameter, scales = "free_x") +
  scale_color_grey() +
  theme_bw() +
  geom_density(aes(color=Variable)) + 
  geom_histogram(aes(y=..density.., color=Variable), position="identity", alpha=0.2, bins=900^0.5)+
  scale_x_log10() +
  labs(x="", y="Density") +
  theme(legend.title = element_blank(),
        legend.position = "top")

tiff("plots/fig9.tiff", res=600, compression = "lzw", height = 8, width = 8, units="in")
plot_grid(p9.1, p9.2, nrow = 2,  align = "v", axis = "l", labels = c("A", "B"))
dev.off()
#ggsave("plots/fig9.jpeg", dpi = 300, height = 8, width = 8, units="in")
