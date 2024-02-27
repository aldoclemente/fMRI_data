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


plot(mesh)
#rgl.postscript("brain_mesh.pdf", fmt="pdf")
snapshot3d(filename = paste0(folder.name,"brain_mesh.png"),
           fmt = "png", width = 800, height = 750, webshot = rgl.useNULL())
rgl.close()

# analysis --------------------------------------------------------------------- 
load("data/data.RData") 

# NA handling 
idx_na <- which(is.na(FC_maps))
X[idx_na] = 0.

# (1): full mesh, full locs
m = 1

observations <- FC_maps[,1:m]
idx_na <- which(is.na(observations))
locations <- mesh$nodes

observations <- observations[-idx_na]
covariates <- X[-idx_na]
locations <- locations[-idx_na,]

idx <- sample(1:nrow(locations), size = 10000)

observations <- observations[idx]
covariates <- covariates[idx]
locations <- locations[idx,]
#lambda=seq(10^-5, 10^-3, length.out=10) #4 patients
#lambda=seq(10^-4, 10^-2, length.out=20) #8 patients
lambda=seq(1e-5, 1e-3, length.out=10)    #20) #12 patients

output_SR_PDE <- fdaPDE::smooth.FEM(observations = as.matrix(observations),
                            covariates = as.matrix(covariates),
                            lambda = lambda,
                            lambda.selection.criterion = "grid", 
                            lambda.selection.lossfunction = "GCV",
                            DOF.evaluation = "stochastic",
                            FEMbasis = FEMbasis)

output_SR_PDE$solution$beta
output_SR_PDE$optimization$lambda_position
# ------------------------------------------------------------------------------
m = 30
observations <- FC_maps[,1:m]
idx_na <-unique(which(is.na(observations),arr.ind = TRUE)[,1]) # delle locazioni che andranno eliminate
locations <- mesh$nodes

observations <- observations[-idx_na,]
covariates <- X[-idx_na,1:m]
locations <- locations[-idx_na,]

# idx <- sample(1:nrow(locations), size = 10000)
# 
# observations <- observations[idx,]
# covariates <- covariates[idx,]
# locations <- locations[idx,]

lambda=seq(1e-4, 1e-2, length.out=12)    #20) #12 patients
start_ <- Sys.time()
output_mixed <- smooth.FEM.mixed(
                           observations = as.matrix(observations), locations = locations,
                           covariates = as.vector(covariates), random_effect = c(1), 
                           lambda = lambda,
                           lambda.selection.criterion = "grid", 
                           lambda.selection.lossfunction = "GCV",
                           DOF.evaluation = "stochastic",
                           FEMbasis = FEMbasis, FLAG_ITERATIVE = TRUE)
time_ <- difftime(Sys.time(), start_, units ="mins")
output_mixed$beta


date_ = unlist(strsplit(as.character(gsub(":","_",gsub(" ","-",Sys.time()))), 
                        split = "[.]"))[1]

folder.name <- paste0(folder.name, date_, "/")
if(!dir.exists(folder.name))
  dir.create(folder.name)
save(output_mixed, file = paste0(folder.name, date_, ".RData"))

nnodes <- nrow(mesh$nodes)
best_lambda <- output_mixed$bestlambda
min.col = min(output_mixed$fit.FEM.mixed$coeff[, best_lambda])
max.col = max(output_mixed$fit.FEM.mixed$coeff[, best_lambda])
for(i in 1:m){ # (nnodes+1):(2*nnodes)
  FEMobject <- FEM(coeff = output_mixed$fit.FEM.mixed$coeff[((i-1)*nnodes+1):(i*nnodes), best_lambda],
                   FEMbasis = FEMbasis)
  
  plot(FEMobject, m=min.col, M=max.col)
  #rgl.postscript("brain_mesh.pdf", fmt="pdf")
  snapshot3d(filename = paste0(folder.name,"patient_", i,".png"),
             fmt = "png", width = 800, height = 750, webshot = rgl.useNULL())
  rgl.close()
}

# RMSE -------------------------------------------------------------------------
rmse <- 0

COV <- matrix(covariates, nrow=nrow(locations), ncol=m)
for(i in 1:m){
  fitted <- COV[,i]%*% as.matrix(output_mixed$beta[, best_lambda]) + 
    eval.FEM(FEM(output_mixed$fit.FEM.mixed$coeff[((i-1)*nnodes+1):(i*nnodes), best_lambda], 
                 FEMbasis), locations) + 
    COV[,i]%*% as.matrix(output_mixed$b_i[i, best_lambda])
  
  rmse <- rmse + (fitted - observations[,i])^2
}

rmse <- sqrt( sum(rmse) / ( m * nrow(locations) ) )
rmse

# 10 folds CROSS VALIDATION ----------------------------------------------------
COV <- matrix(covariates, nrow=nrow(locations), ncol=m)
K <- 10L
seed <- as.integer(12345)
folds <- list()
num_obs_kFold <- vector(mode="integer", length=K)
data_locs <- data.frame(locations = locations)
data_obs <- data.frame(observations = observations)
data_covs <- data.frame(COV = COV)

names(data_obs) <- 1:m
names(data_covs) <- 1:m
sample_ <- sample(1:nrow(data_locs))

data_locs <- data_locs[sample_,]
data_covs <- data_covs[sample_,]
data_obs <- data_obs[sample_,]

num_data = round(nrow(data_locs)/K)
num_obs_kFold[1:(K-1)] <- rep(num_data, times=(K-1))

for(i in 1:(K-1)){
  folds[[i]] = list(data_locs = data_locs[(1 + num_data*(i-1)):(num_data*i),],
                    data_covs = data_covs[(1 + num_data*(i-1)):(num_data*i),],
                    data_obs  = data_obs[(1 + num_data*(i-1)):(num_data*i),])
  
}
folds[[K]] = list(data_locs = data_locs[(num_data*(K-1) + 1):nrow(data_locs),],
                   data_covs = data_covs[(num_data*(K-1) + 1):nrow(data_locs),],
                   data_obs  = data_obs[(num_data*(K-1) + 1):nrow(data_locs),])
num_obs_kFold[K] <- nrow(folds[[K]]$data_locs)
storage.mode(num_obs_kFold) <- "integer" 

num_obs_kFold
CV_error <- vector(mode="numeric", length = K)
nnodes <- nrow(mesh$nodes)
for(k in 1:K){
  cat("----------------- ", k ," / ", K , "-----------------\n")
  train_data = list()
  
  train_locs <- matrix(nrow=0, ncol=3)
  train_covs <- matrix(nrow=0, ncol=m)
  train_obs <- matrix(nrow=0, ncol=m)
  
  test_locs <- matrix(nrow=0, ncol=3)
  test_covs <- matrix(nrow=0, ncol=m)
  test_obs <- matrix(nrow=0, ncol=m)
  
  # training & test
  for(i in 1:K){
    if( i == k){
      test_locs = folds[[i]]$data_locs
      test_covs = folds[[i]]$data_covs
      test_obs = folds[[i]]$data_obs
    }else{
      train_data = rbind(train_data, folds[[i]])
      
      train_locs = rbind(train_locs, folds[[i]]$data_locs)
      train_covs = rbind(train_covs, folds[[i]]$data_covs)
      train_obs =  rbind(train_obs, folds[[i]]$data_obs)
    }
  }
  
  lambda=seq(1e-4, 1e-2, length.out=12)    #20) #12 patients

  invisible(capture.output(output_mixed <-  smooth.FEM.mixed(
    observations = as.matrix(train_obs), locations = as.matrix(train_locs),
    covariates = as.vector(as.matrix(train_covs)), random_effect = c(1), 
    lambda = lambda,
    lambda.selection.criterion = "grid", 
    lambda.selection.lossfunction = "GCV",
    DOF.evaluation = "stochastic",
    FEMbasis = FEMbasis, FLAG_ITERATIVE = TRUE)))
  
  best_lambda <- output_mixed$bestlambda
  rmse <- 0
  
  # RMSE
  for(j in 1:m){
    fitted <- test_covs[,j]%*% as.matrix(output_mixed$beta[, best_lambda]) + 
      eval.FEM(FEM(output_mixed$fit.FEM.mixed$coeff[((j-1)*nnodes+1):(j*nnodes), best_lambda], 
                   FEMbasis), test_locs) + 
      test_covs[,j]%*% as.matrix(output_mixed$b_i[j, best_lambda])
    
    rmse <- rmse + (fitted - test_obs[,j])^2
  }
  
  CV_error[k] <- sqrt( sum(rmse) / ( m * nrow(test_locs) ) )
}

save(CV_error, file = paste0(folder.name, "CV_error.RData"))

png(paste0(folder.name, "CV_error.png"))
boxplot(CV_error, col="lightgray", main="CV error")
dev.off()