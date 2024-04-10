if(!require(fdaPDEmixed)){
  devtools::install_github(repo="aldoclemente/fdaPDEmixed", ref="main") 
}

pacman::p_load("R.matlab", "fdaPDEmixed")
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

date_ = unlist(strsplit(as.character(gsub(":","_",gsub(" ","-",Sys.time()))), 
                        split = "[.]"))[1]
folder.name <- paste0(folder.name, date_, "/")

if(!dir.exists(folder.name))
  dir.create(folder.name)

#rgl.postscript("brain_mesh.pdf", fmt="pdf")
plot(mesh)
snapshot3d(filename = paste0(folder.name,"brain_mesh.png"),
           fmt = "png", width = 800, height = 750, webshot = rgl.useNULL())
rgl.close()

# Gordon (2016) parcellization, see Lila et al.  
gordon_parcellation = read.csv("data/gordon_parcellation2016.csv")
dim(gordon_parcellation) 
na_ <- as.logical(gordon_parcellation$na)
plot(mesh, NA_ = na_)
snapshot3d(filename = paste0(folder.name,"brain_mesh.png"),
           fmt = "png", width = 800, height = 750, webshot = rgl.useNULL())
rgl.close()

roi <- as.logical(gordon_parcellation$Parcel_1) # precuneo !
plot.mesh.2.5D(mesh, ROI = roi, NA_=na_)
snapshot3d(filename = paste0(folder.name, "brain_mesh_roi.png"),
           fmt = "png", width = 800, height = 750, webshot = rgl.useNULL())
rgl.close()

# analysis --------------------------------------------------------------------- 
load("data/gordon_parcellation/FCmaps.RData")
load("data/gordon_parcellation/thickness.RData")

# NA handling 
idx_na <- which(is.na(FCmaps))
X[idx_na] = 0.

# ------------------------------------------------------------------------------
 
m = ncol(FC_maps) # num subjects
observations <- FC_maps
idx_na <-unique(which(is.na(observations),arr.ind = TRUE)[,1]) # delle locazioni che andranno eliminate
locations <- mesh$nodes

observations <- observations[-idx_na,]
covariates <- X[-idx_na,1:m]
locations <- locations[-idx_na,]

lambda=seq(1e-4, 1e-2, length.out=12)    

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

fill_col <- viridis::viridis(2, begin=0.25, end=0.95)[1]
mai_ = par("mai")
mai_[2] = mai_[2] + 0.075
pdf(paste0(folder.name, "CV_error.png"), family = "serif", width = 7, height = 7)
boxplot(CV_error, col=fill_col, main="CV error", xlab="", ylab="",
        cex.lab = 2, cex.axis = 2, cex.main = 2,
        col=fill_col)
dev.off()
