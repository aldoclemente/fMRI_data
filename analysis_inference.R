if(!require(fdaPDEISCHIA)){
  devtools::install_github(repo="aldoclemente/fdaPDEISCHIA", ref="main") 
}

pacman::p_load("R.matlab", "fdaPDEISCHIA")
source("utils.R")

# mesh -------------------------------------------------------------------------
nodes <- readMat("data/vertices.mat")
faces <- readMat("data/faces.mat")

nodes <- nodes$vertices
faces <- faces$faces

mesh <- create.mesh.2.5D(nodes = nodes, triangles = faces)
# plot(mesh)
FEMbasis <- create.FEM.basis(mesh)

folder.name <- "data/application/"
if(!dir.exists(folder.name))
  dir.create(folder.name)

# BOOTSTRAP --------------------------------------------------------------------
load("data/data.RData") 
m = 30
observations <- FC_maps[,1:m]
idx_na <-unique(which(is.na(observations),arr.ind = TRUE)[,1]) # delle locazioni che andranno eliminate
locations <- mesh$nodes

observations <- observations[-idx_na,]
covariates <- X[-idx_na,1:m]
locations <- locations[-idx_na,]

N = nrow(observations)

nsim <- 1000 # consigliato per costruire IC
betas <- matrix(nrow=1, ncol=nsim)
b <- matrix(nrow=m, ncol=nsim)

lambda=seq(1e-3, 1e-1, length.out=12)
for(i in 1:nsim){
  cat("----------------- ", i ," / ", nsim , "-----------------\n")
  idx <- sample(1:N, size=N, replace = TRUE)
  locs <- locations[idx,]
  obs <- observations[idx,]
  covs<- covariates[idx,]
  
  #lambda=seq(1e-3, 1e-1, length.out=12)    #20) #12 patients
  start_ <- Sys.time()
  invisible(capture.output(output_mixed <- smooth.FEM.mixed(
    observations = as.matrix(obs), locations = locs,
    covariates = as.vector(covs), random_effect = c(1), 
    lambda = lambda[12],
    #lambda.selection.criterion = "grid", 
    #lambda.selection.lossfunction = "GCV",
    #DOF.evaluation = "stochastic",
    FEMbasis = FEMbasis, FLAG_ITERATIVE = TRUE)))
  time_ <- difftime(Sys.time(), start_, units ="mins")
  cat("elapsed time: ", time_, "\n")
  #best_lambda <- output_mixed$bestlambda
  betas[i] <- output_mixed$beta
  b[,i] <- output_mixed$b_i
}

filename <- "data/inference_bootstrap.RData"
save(betas, b, file=filename)
