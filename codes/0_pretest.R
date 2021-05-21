# Set working directory (need to be customized) --------------------------------
#setwd("C:/Users/nan_1/Documents/bup-cat") 

if(!dir.exists("outputs")) dir.create("outputs")
if(!dir.exists("plots")) dir.create("plots")
if(!dir.exists("plots/suppl")) dir.create("plots/suppl")

# Create mod.exe ---------------------------------------------------------------
source("MCSim/function.R")
set_PATH()
makemod()

# Compile model code -----------------------------------------------------------
model <- "bupropion_CAT.model" 
makemcsim(model)

# Test in-vitro modeling (single chain) ----------------------------------------
# Setting small iteration in testing
input <- "bupropion_MCMC_vitro.in"
set.seed(1111)
system.time(y <- mcsim(model, input)) 

# Test in-vivo modeling (single chain) -----------------------------------------
# Setting small iteration in testing
input <- "bupropion_MCMC_vivo.in"
set.seed(2222)
system.time(y <- mcsim(model, input)) 

# Test in-vivo modeling (parallel) ---------------------------------------------
library(foreach)
library(doParallel)

current.files <- list.files()
detectCores()
cores <- 4 
cl <- makeCluster(cores)
registerDoParallel(cl)

strt<-Sys.time()
out <- foreach(i = 1:cores) %dopar% 
  {
    set.seed(i+10)
    mcsim(model = model, input = "bupropion_MCMC_vitro.in", parallel = T)  
  }
print(Sys.time()-strt)

new.files <- setdiff(list.files(), current.files)
to.remove <- new.files[grep('.kernel|.in', new.files)]
file.remove(to.remove)
out.files <- setdiff(list.files(), current.files)

for(i in 1:4){
  file.copy(out.files[i], paste0("outputs/", out.files[i]))
  file.remove(out.files[i])
}
