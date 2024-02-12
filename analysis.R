if(!require(fdaPDEKIM)){
  devtools::install_github("aldoclemente/fdaPDEKIM", ref="master") 
}
pacman::p_load("R.matlab")

# mesh -------------------------------------------------------------------------
nodes <- readMat("data/vertices.mat")
faces <- readMat("data/faces.mat")

nodes <- nodes$vertices
faces <- faces$faces

mesh <- create.mesh.2.5D(nodes = nodes, triangles = faces)
plot(mesh)
FEMbasis <- create.FEM.basis(mesh)

# analysis --------------------------------------------------------------------- 
load("data/data.RData") 

GCVFLAG=T
GCVMETHODFLAG='Stochastic'

# NA handling 
idx_na <- which(is.na(FC_maps))
X[idx_na] = 0.

# (1): full mesh, full locs
n = 12
observations <- FC_maps[,1:n]
covariates <- as.vector(X[, 1:n]) # attenzione 

#lambda=seq(10^-5, 10^-3, length.out=10) #4 patients
#lambda=seq(10^-4, 10^-2, length.out=20) #8 patients
lambda=seq(0.001, 0.005, length.out=20) #12 patients

output_CPP <- fdaPDEKIM:::smooth.FEM.mixed(observations = observations,
                              covariates = covariates,
                              random_effect =c(1),
                              FEMbasis = FEMbasis,
                              lambda = lambda,
                              GCV = GCVFLAG,
                              GCVmethod = GCVMETHODFLAG,
                              TESTFLAG = FALSE)
