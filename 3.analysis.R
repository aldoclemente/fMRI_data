if(!require(fdaPDEmixed)){
  devtools::install_github(repo="aldoclemente/fdaPDEmixed", ref="main") 
}

system("whoami > user.txt")
user <- read.table("user.txt", header = FALSE)[,1]
unlink("user.txt")
.libPaths( c( paste0("/home/", user, "/R/"), .libPaths() ) ) 
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
{
plot(mesh)
snapshot3d(filename = paste0(folder.name,"brain_mesh.png"),
           fmt = "png", width = 800, height = 750, webshot = rgl.useNULL())
rgl.close()
}

# Gordon (2016) parcellization, see Lila et al.  
gordon_parcellation = read.csv("data/gordon_parcellation2016.csv")
dim(gordon_parcellation) 
na_ <- as.logical(gordon_parcellation$na)
{
plot(mesh, NA_ = na_)
snapshot3d(filename = paste0(folder.name,"brain_mesh.png"),
           fmt = "png", width = 800, height = 750, webshot = rgl.useNULL())
rgl.close()
}

roi <- as.logical(gordon_parcellation$Parcel_1) # precuneo !

{
plot.mesh.2.5D(mesh, ROI = roi, NA_=na_)
snapshot3d(filename = paste0(folder.name, "brain_mesh_roi.png"),
           fmt = "png", width = 800, height = 750, webshot = rgl.useNULL())
rgl.close()
}
# analysis --------------------------------------------------------------------- 
load("data/gordon2016_FCmaps.RData")
load("data/thickness.RData")

# ------------------------------------------------------------------------------
m = ncol(FCmaps) # num subjects
#m = 2 
observations <- FCmaps
na_idx <-unique(which(is.na(observations[,1]))) 
#locations <- mesh$nodes

observations <- observations[,1:m]
covariates <- thickness[,1:m]

covariates[na_idx, ] <- 1e-8
# locations <- locations[-na_idx,]

#lambda=seq(1e-4, 1e1, length.out=30)  # NO
#lambda=seq(1e-4, 1e-1, length.out=40) # Quasi, siamo sicuri che la GCV funzioni ?
lambda=seq(1e-4, 1e-2, length.out=40)

start_ <- Sys.time()
invisible(capture.output(output_mixed <- smooth.FEM.mixed(
  observations = as.matrix(observations), #locations = locations,
  covariates = as.vector(covariates), random_effect = c(1), 
  lambda = lambda,
  lambda.selection.criterion = "grid", 
  lambda.selection.lossfunction = "GCV",
  DOF.evaluation = "stochastic",
  FEMbasis = FEMbasis, FLAG_ITERATIVE = TRUE)))
time_ <- difftime(Sys.time(), start_, units ="mins")
cat("time elapsed: ", round(as.numeric(time_), 2), " mins \n")
output_mixed$beta

save(output_mixed, lambda, file = paste0(folder.name, date_, ".RData"))

png(paste0(folder.name, "GCV.png" ))
plot(log10(lambda), output_mixed$GCV, pch=16, cex=1.5,
     xlab = expression(log[10](lambda)), ylab="GCV")
dev.off()
nnodes <- nrow(mesh$nodes)
best_lambda <- output_mixed$bestlambda
# best_lambda <- 1
# lambda <- lambda[best_lambda]
min.col = min(output_mixed$fit.FEM.mixed$coeff[, best_lambda])
max.col = max(output_mixed$fit.FEM.mixed$coeff[, best_lambda])

colorscales = list(rainbow = jet.col, viridis = viridis::viridis)
colorbar_limits = list(same_limits = "same_limits", different_limits = "different_limits")

for(p in names(colorscales)){
  imgs.name <- paste0(folder.name, p, "/")
    if(!dir.exists(imgs.name))
      dir.create(imgs.name)
   
  for(lims in names(colorbar_limits)){
    imgs.name <- paste0(folder.name, p, "/", lims, "/")
    if(!dir.exists(imgs.name)) dir.create(imgs.name)
  
    for(i in 1:m){ # (nnodes+1):(2*nnodes)
      coeff <- output_mixed$fit.FEM.mixed$coeff[((i-1)*nnodes+1):(i*nnodes), best_lambda]
      coeff[na_idx] <- NA
      FEMobject <- FEM(coeff = coeff,
                       FEMbasis = FEMbasis)
      if( lims == "different_limits" ){
          MIN.COL = MAX.COL = limits = NULL;
        }else{
          MIN.COL = min.col; MAX.COL = max.col; limits=c(min.col, max.col);
        }
      
      
      if(!dir.exists(paste0(imgs.name, "estimates/"))) dir.create(paste0(imgs.name, "estimates/"))
      plot(FEMobject, m=MIN.COL, M=MAX.COL, colorscale = colorscales[[p]])
        #rgl.postscript("brain_mesh.pdf", fmt="pdf")
      snapshot3d(filename = paste0(paste0(imgs.name, "estimates/"),
                                   "patient_", i,".png"),
                   fmt = "png", width = 800, height = 750, webshot = rgl.useNULL())
      rgl.close()
      
      if(!dir.exists(paste0(imgs.name, "colorbars/"))) dir.create(paste0(imgs.name, "colorbars/"))
      
      if(lims == "same_limits" & i == 1){
        plot.colorbar(FEMobject, limits, colorscale =  colorscales[[p]], 
                    file = paste0(paste0(imgs.name, "colorbars/"), 
                                  "colorbar") )
      }else if(lims == "different_limits"){
        plot.colorbar(FEMobject, limits, colorscale =  colorscales[[p]], 
                    file = paste0(paste0(imgs.name, "colorbars/"), 
                                  "patient_", i ,"_colorbar") )
      }
    }
  }
}

# mean field (?) ---------------------------------------------------------------
mean_coeff <- matrix(output_mixed$fit.FEM.mixed$coeff[,best_lambda], nrow=nnodes, ncol=m)
mean_coeff <- rowMeans(mean_coeff, na.rm = T) 

for(p in names(colorscales)){
    imgs.name <- paste0(folder.name, p, "/")
    mean_coeff[na_idx] <- NA
    FEMobject <- FEM(coeff = mean_coeff,
                     FEMbasis = FEMbasis)
    plot(FEMobject, colorscale = colorscales[[p]])
    #rgl.postscript("brain_mesh.pdf", fmt="pdf")
    snapshot3d(filename = paste0(paste0(imgs.name, "mean_estimate.png")),
               fmt = "png", width = 800, height = 750, webshot = rgl.useNULL())
    rgl.close()
    
    plot.colorbar(FEMobject, colorscale =  colorscales[[p]], 
                  file = paste0(paste0(imgs.name, "mean_colorbar") ))
}
# RMSE -------------------------------------------------------------------------
rmse <- matrix(0, nrow=nrow(mesh$nodes[-na_idx,]), ncol=30)

for(i in 1:m){
  coeff <- output_mixed$fit.FEM.mixed$coeff[((i-1)*nnodes+1):(i*nnodes), best_lambda]
  coeff[na_idx] <- NA
  fitted <- covariates[-na_idx,i]%*% as.matrix(output_mixed$beta[-na_idx, best_lambda]) + 
    eval.FEM(FEM(coeff, FEMbasis), mesh$nodes[-na_idx,]) + 
    covariates[-na_idx,i]%*% as.matrix(output_mixed$b_i[i, best_lambda])
  
  rmse[,i] <- (fitted - observations[-na_idx,i])^2
}

sqrt(sum(colSums(rmse, na.rm = T))/(m*nrow(mesh$nodes[-na_idx,])))

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
