if(!require(fdaPDEISCHIA)){
  devtools::install_github(repo="aldoclemente/fdaPDEISCHIA", ref="main") 
}

pacman::p_load("R.matlab", "fdaPDEISCHIA")

# mesh -------------------------------------------------------------------------
nodes <- readMat("data/vertices.mat")
faces <- readMat("data/faces.mat")

nodes <- nodes$vertices
faces <- faces$faces

mesh <- create.mesh.2.5D(nodes = nodes, triangles = faces)
# plot(mesh)
FEMbasis <- create.FEM.basis(mesh)

# analysis --------------------------------------------------------------------- 
load("data/data.RData") 

# NA handling 
idx_na <- which(is.na(FC_maps))
X[idx_na] = 0.

# (1): full mesh, full locs
n = 1

observations <- FC_maps[,1:n]
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
n = 30
observations <- FC_maps[,1:n]
idx_na <-unique(which(is.na(observations),arr.ind = TRUE)[,1]) # delle locazioni che andranno eliminate
locations <- mesh$nodes

observations <- observations[-idx_na,]
covariates <- X[-idx_na,1:n]
locations <- locations[-idx_na,]

idx <- sample(1:nrow(locations), size = 10000)

observations <- observations[idx,]
covariates <- covariates[idx,]
locations <- locations[idx,]

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

# range(output_SR_PDE_KIM$fit.FEM$coeff)
# range(output_SR_PDE_mixed$fit.FEM.mixed$coeff)

# 1
fitted_KIM <- eval.FEM(output_SR_PDE_KIM$fit.FEM, locations = mesh$nodes) + 
  as.matrix(covariates) %*% output_SR_PDE_KIM$beta
sqrt(sum(fitted_KIM[-idx_na]-observations[-idx_na])^2)
# 2
fitted_mixed <- eval.FEM(FEM(output_SR_PDE_mixed$fit.FEM$coeff, FEMbasis),
                         locations = mesh$nodes) +  as.matrix(covariates) %*% output_SR_PDE_mixed$beta
sqrt(sum(fitted_mixed[-idx_na]-observations[-idx_na])^2)


