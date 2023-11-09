library(dplyr)
library(foreach)
library(doParallel)
library(boot)
library(INLA)
library(tidyverse)

v <- Sys.getenv("SLURM_ARRAY_TASK_ID", NA)
vn <- as.integer(v)
ncores <- Sys.getenv("SLURM_CPUS_PER_TASK")
ncores <- as.numeric(ncores)

inla.setOption(num.threads = ncores)
pdf_outcome <- readRDS("pdf_outcome.rds")

params <- read.csv("cases.csv")
myparams <- params[params$case==vn,]


OR_TX = myparams$OR_TX #alter this for different treatment effect odds ratio settings
burn.in <- myparams$burn.in #alter this for different initial enrollment before analysis settings
prob_sup <- myparams$prob_sup #alter this for the superiority threshold
prob_futi <- myparams$prob_fut #alter this for the futility threshold
N_patient <- 5000 #alter this for different maximum sample size settings
recruit_rate <- 190

VFDs <- unique(pdf_outcome["VFDs"])
VFDs <- as.numeric(VFDs[,1])
treat_futi <- "TX" 
states <- c("states")
TxARM_n <- c("TX","UC")
p_UC <- as.numeric(pdf_outcome[, "p_prob"])
p_TX <- as.numeric(pdf_outcome[, as.character(OR_TX)])


DB_simulation <- function(N_patient, 
                          burn.in, 
                          prob_sup, 
                          prob_futi,
                          treat_futi,
                          OR_TX, 
                          recruit_rate, 
                          p_UC,
                          p_TX){
  
  Result = array(dim = c(4, 1),
                 dimnames = list(c("Sample.size", "TxARM", "Result", "Tx.effect"), "states"))
  there.has.been.an.error <- FALSE
  #n_state <- cbind(recruit_rate); colnames(n_state) <- states
  np <- n <- num_patient <- recruit_rate
  TxARM <- factor(sample(TxARM_n , n, replace = TRUE))
  #State <- factor(sample(states , n, replace = TRUE))
  #n_state <- cbind(sum(State == "HSE")); colnames(n_state) <- states
  
  outcome <- vector("numeric", length(TxARM))
  #outcome[TxARM == "UC" & State == "LSE"] <- sample(VFDs, sum(TxARM == "UC"& State == "LSE"), replace = TRUE, prob = p_UC_LSE)
  outcome[TxARM == "UC"] <- sample(VFDs, sum(TxARM == "UC"), replace = TRUE, prob = p_UC)
  #outcome[TxARM == "TX"& State == "LSE"] <- sample(VFDs, sum(TxARM == "TX"& State == "LSE"), replace = TRUE, prob = p_TX_LSE)
  outcome[TxARM == "TX"] <- sample(VFDs, sum(TxARM == "TX"), replace = TRUE, prob = p_TX)
  sample_data <- data.frame(outcome, TxARM)
  ref = "UC"
  sample_data$TxARM = relevel(sample_data$TxARM, ref = ref)
  ind <- 0
  
  #----half t prior-----------------------------------
  #https://becarioprecario.bitbucket.io/inla-gitbook/ch-priors.html#sec:priors
  HT.prior = "expression:
  sigma = exp(-theta/2);
  nu = 3;
  tau = 7;
  log_dens = 0 - 0.5 * log(nu * pi) - (-0.1207822);
  log_dens = log_dens - 0.5 * (nu + 1) * log(1 + (sigma * sigma) / (nu * (tau * tau)));
  log_dens = log_dens - log(2) - theta / 2;
  return(log_dens);
  "
  h.t = list(theta = list(prior = HT.prior))
  
  #list(theta=list(prior="loggamma", param=c(alpha, beta)))
  
  #lc <- inla.make.lincombs(State = matrix(c(1,0, 0, 1), nrow = 2), TxARMTX = c(1, 1), StateLSE = c(0, 0))
  
  while (TRUE) {
    if(sum((num_patient > burn.in) * (ind == 0)) >= 1){
      if(length(unique(sample_data$outcome)) != 9){
        
        sample_data_adj <- sample_data
        sample_data_adj$outcome <- as.numeric(factor(sample_data$outcome))
        tryCatch(
          {
            model <- inla(outcome ~  TxARM,
                          family = "pom",
                          data = sample_data_adj)
          },
          error = function(e){
            print("Error in INLA")
            there.has.been.an.error <<- TRUE
          }
        )
      }
      if(length(unique(sample_data$outcome)) == 9){
        
        tryCatch(
          {
            model <- inla(outcome ~  TxARM,
                          family = "pom",
                          data = sample_data,
                          verbose=FALSE, control.compute = list(config=TRUE),
                          control.family = list(hyper = list(theta1 = list(param = 100))))
          },
          error = function(e){
            print("Error in INLA")
            there.has.been.an.error <<- TRUE
          }
        )
      }
      
      if(!there.has.been.an.error){
        m0 <- model$marginals.fixed$TxARMTX
        trt.effect <- exp(model$summary.fixed["TxARMTX",1])
        #marg_HSE <- test$marginals.fixed$TxARMUC
        rm(model)
        Prob.differ.futi <- inla.pmarginal(log(1.2), m0)
        post.p <- cbind("treat.opt" = TxARM_n,
                        "states" = c(1 - inla.pmarginal(log(1), m0), inla.pmarginal(log(1), m0)))
        
        
        #n.state <- as.vector(n_state[, states[j]]) 
        #ind.state <- as.vector(ind[,states[j]])
        
        n.state <- num_patient
        ind.state <- ind
        
        if (n.state > burn.in){
          if (ind.state == 0) {
            probs <- as.numeric(post.p[, states])
            treatments <- post.p[, "treat.opt"]
            
            if  (sum(probs > prob_sup) > 0) {
              treat <- treatments[which(probs > prob_sup)] 
              probs <- probs[which(probs > prob_sup)]
              Result["Sample.size", ] = n.state
              Result["Result", ] = "superiority" 
              Result["TxARM", ]  = treat
              Result["Tx.effect", ]  = trt.effect
              ind = 1
            }
            if  (sum(probs > prob_sup) == 0) {
              if  (Prob.differ.futi > prob_futi) {
                Result["Sample.size", ] = n.state
                Result["Result", ] = "futility" 
                Result["TxARM", ]  =  treat_futi
                Result["Tx.effect", ]  = trt.effect
                ind = 1
              }
              else
              {Result["Sample.size", ] = n.state
              Result["Result", ] = "Nothing" 
              Result["TxARM", ]  = "Nothing"
              Result["Tx.effect", ]  = trt.effect
              }  
            }
            if (n.state >= N_patient){
              Result["Sample.size", ] = n.state
              Result["Result", ] = "Nothing" 
              Result["TxARM", ]  = "Nothing"
              Result["Tx.effect", ]  = trt.effect
              ind = 1
            }  
          }
        } 
        else { 
          Result["Sample.size", ] = n.state
          Result["Result", ] = "NA" 
          Result["TxARM", ]  = "NA"
          Result["Tx.effect", ]  = trt.effect
        }
        
      }
    }
    n_new <- 0
    if (ind == 1){
      n_new <- 0
    }
    else {
      n_new <- recruit_rate
    } 
    
    
    #n_LSE <- as.vector(n_state[, "LSE"])  
    #n_HSE <- as.vector(n_state[, "HSE"])
    #new_LSE <- as.vector(n_new[, "LSE"])  
    #new_HSE <- as.vector(n_new[, "HSE"])
    n <- n_new
    if (num_patient >= N_patient | n == 0){break}
    TxARM <- factor(sample(TxARM_n , n, replace = TRUE))
    #State <- factor(sample(states , n, replace = TRUE, prob = c(new_LSE/n, new_HSE/n)))
    outcome <- vector("numeric", length(TxARM))
    outcome[TxARM == "UC" ] <- sample(VFDs, sum(TxARM == "UC"), replace = TRUE, prob = p_UC)
    #outcome[TxARM == "UC" & State == "HSE"] <- sample(VFDs, sum(TxARM == "UC" & State == "HSE"), replace = TRUE, prob = p_UC_HSE)
    outcome[TxARM == "TX"] <- sample(VFDs, sum(TxARM == "TX"), replace = TRUE, prob = p_TX)
    #outcome[TxARM == "TX"& State == "HSE"] <- sample(VFDs, sum(TxARM == "TX"& State == "HSE"), replace = TRUE, prob = p_TX_HSE)
    new_data <- data.frame(outcome, TxARM)
    sample_data <- rbind(sample_data, new_data)
    sample_data$TxARM = relevel(sample_data$TxARM, ref = ref)
    #n_state <- cbind(sum(sample_data$State == "LSE"),
    #                 sum(sample_data$State == "HSE")); colnames(n_state) <- states
    np <- num_patient <- nrow(sample_data)
    if(there.has.been.an.error){break}
  }
  return(Result)
}

#parallel computing code
n_sim <-2
registerDoParallel(cores=ncores)
Sim_Res <- foreach(k = 1:n_sim, .combine = cbind,
                   .packages = c("boot", "dplyr", "INLA")
) %dopar% { 
  set.seed(k + 123)
  DB_simulation(N_patient = N_patient, 
                burn.in = burn.in, 
                prob_sup = prob_sup, 
                prob_futi = prob_futi,
                treat_futi = treat_futi,
                OR_TX = OR_TX, 
                recruit_rate = recruit_rate, 
                p_UC = p_UC,
                p_TX = p_TX)
}

results <- as.data.frame(Sim_Res)

fname = paste(vn, "OR", OR_TX, "int", burn.in, "sup", prob_sup, "fut",prob_futi, sep="")

write.csv(results,paste0(fname,".csv"))
#write.csv(results, file = "OR1.3_n800_sup995_fut90.csv") 