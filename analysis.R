if(!require(fdaPDEKIM)){
  devtools::install_github(repo="aldoclemente/fdaPDEKIM", ref="main") 
}

pacman::p_load("R.matlab", "fdaPDEKIM")

# mesh -------------------------------------------------------------------------
nodes <- readMat("data/vertices.mat")
faces <- readMat("data/faces.mat")

nodes <- nodes$vertices
faces <- faces$faces

mesh <- fdaPDEKIM::create.mesh.2.5D(nodes = nodes, triangles = faces)
#plot(mesh)
FEMbasis <- fdaPDEKIM::create.FEM.basis(mesh)

# mesh fdaPDE (CRAN)
mesh_CRAN <- fdaPDE::create.mesh.2.5D(nodes=nodes, triangles = faces)
FEMbasis_CRAN <- fdaPDE::create.FEM.basis(mesh_CRAN)

# analysis --------------------------------------------------------------------- 
load("data/data.RData") 

# NA handling 
idx_na <- which(is.na(FC_maps))
X[idx_na] = 0.

# (1): full mesh, full locs
n = 2
idx <- sample(1:nrow(FC_maps), size = 100)

observations <- FC_maps[idx,1:n]
covariates <- as.vector(X[idx, 1:n]) # attenzione 
locations <- mesh$nodes[idx,]

#lambda=seq(10^-5, 10^-3, length.out=10) #4 patients
#lambda=seq(10^-4, 10^-2, length.out=20) #8 patients
lambda=seq(0.001, 0.005, length.out=10)    #20) #12 patients

# CRASHA
# output_SR_PDE <- fdaPDE::smooth.FEM(observations = observations,
#                                     covariates = covariates, 
#                                     lambda = lambda[1],
#                                     FEMbasis = FEMbasis_CRAN)

# mi aspetto che output_SR_PDE_KIM e output_SR_PDE_mixed con random_effect = NULL 
# restituiscano stessa stima per beta...

# output_SR_PDE_KIM <- fdaPDEKIM:::smooth.FEM(observations = observations,
#                                                     covariates = covariates,
#                                                     FEMbasis = FEMbasis,
#                                                     lambda = lambda
#                                                       )
# output_SR_PDE_KIM$beta

# n = 2
# observations <- FC_maps[,1:n]
# covariates <- as.vector(X[, 1:n]) # attenzione 
GCVFLAG=T
GCVMETHODFLAG='Stochastic'
TESTFLAG <- TRUE
output_SR_PDE_mixed <- fdaPDEKIM:::smooth.FEM.mixed(observations = observations,
                              covariates = covariates, locations = locations,
                              random_effect = c(1),
                              FEMbasis = FEMbasis,
                              lambda = lambda, GCV = GCVFLAG, GCVmethod = GCVMETHODFLAG, 
                              TESTFLAG = TESTFLAG)
output_SR_PDE_mixed$beta

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

# ------------------------------------------------------------------------------

# (2): simplified mesh (10000), full locs

load("data/10000_mesh.RData")
mesh_simp <- fdaPDEKIM::create.mesh.2.5D(nodes=obj2$mesh$nodes, triangles = obj2$mesh$triangles)
plot(mesh_simp)
FEMbasis_simp <- create.FEM.basis(mesh_simp)
# lento .
#projection.points.2.5D(mesh_simp, locations = mesh$nodes[,1:3])

locations <- obj2$locations

load("data/data.RData")
m = 3
observations = FC_maps[,1:m]
covariates = X[1:(dim(observations)[1]*m)]

# but quite different from (1) full mesh, full locs....
lambda=seq(10^-3, 10^-1, length.out=10) #4 patients
# lambda=seq(10^-3, 10^-1, length.out=20) #8 patients
# lambda=seq(10^-3, 10^-1, length.out=20) #12 patients
# lambda=seq(10^-3, 10^-1, length.out=20) #14 patients
# lambda=seq(10^-3, 10^-1, length.out=20) #18 patients
# lambda=seq(0.03, 0.07, length.out=20) #20 patients
# lambda=seq(0.02, 0.07, length.out=20) #30 patients

output_CPP2 = fdaPDEKIM:::smooth.FEM.mixed(locations = locations,
                               observations = observations,
                               covariates = covariates,
                               random_effect =NULL,
                               FEMbasis = FEMbasis_simp,
                               lambda = lambda, 
                               GCV = TRUE,
                               GCVmethod = "Stochastic",  TESTFLAG = TRUE)


# 30 patients
output_CPP2$bestlambda #16
lambda[output_CPP2$bestlambda] #0.05947368

output_CPP2$beta[,output_CPP2$bestlambda]
# 0.005071381
output_CPP2$b_i[,output_CPP2$bestlambda]
#         b_11          b_21          b_31          b_41          b_51          b_61          b_71          b_81 
# 6.236198e-03  8.868908e-03 -5.366381e-03  6.934197e-03 -1.989608e-03 -1.240056e-03 -3.290210e-03 -2.518831e-03 
#         b_91         b_101         b_111         b_121         b_131         b_141         b_151         b_161 
# 1.207917e-02 -2.765547e-03 -2.929114e-03  7.691843e-03 -7.339814e-03 -7.028336e-04  8.638687e-03 -5.646423e-03 
#         b_171         b_181         b_191         b_201         b_211         b_221         b_231         b_241 
# -6.880857e-S03  4.561869e-03 -9.374962e-04 -1.570609e-05  1.579893e-03  2.512863e-03 -8.233994e-03  3.632962e-03 
#       b_251         b_261         b_271         b_281         b_291         b_301 
# 9.585140e-03 -1.663809e-03 -4.204134e-03 -1.872371e-02  4.033723e-03 -1.906925e-03 


plot(FEM.mixed(output_CPP2$fit.FEM.mixed$coeff[,output_CPP2$bestlambda], 
               num_units=wanted_unit,
               FEMbasis2))


obs_index=c(10,11,12,17,21,27)
test_cov = c(covariates[((10-1)*32492+1):(10*32492)],
             covariates[((11-1)*32492+1):(11*32492)],
             covariates[((12-1)*32492+1):(12*32492)],
             covariates[((17-1)*32492+1):(17*32492)],
             covariates[((21-1)*32492+1):(21*32492)],
             covariates[((27-1)*32492+1):(27*32492)])

output_CPP = smooth.FEM.mixed(observations = observations[,obs_index],
                              covariates = test_cov,
                              random_effect =c(1),
                              FEMbasis = FEMbasis,
                              lambda = lambda,
                              GCV = GCVFLAG,
                              GCVmethod = GCVMETHODFLAG)

output_CPP2 = smooth.FEM.mixed(locations = projected.points2,
                               observations = observations[,obs_index],
                               covariates = test_cov,
                               random_effect =c(1),
                               FEMbasis = FEMbasis2,
                               lambda = lambda,
                               GCV = GCVFLAG,
                               GCVmethod = GCVMETHODFLAG)
##### meshsimplification location vs points.projection.2.5D ####
### they are not same....
load("C:/Users/JIYOUNG KIM/Documents/VirtualMachine/ubuntu/MeshDataSimplification/data/10000_brain_mesh.RData")
dim(obj2$locations) #32492     3
head(obj2$locations)
#             [,1]       [,2]     [,3]
# [1,]  -4.496712 -43.769893 32.04211
# [2,] -18.425068 -41.241051 69.17323
# [3,] -52.468578  -6.555558 46.29697
# [4,]  -7.416639  25.434665 60.87541
# [5,] -17.897900 -91.522003 20.73351
# [6,] -48.365742 -52.942780 48.91669

load("C:/Users/JIYOUNG KIM/Desktop/tests - mixed/RData/old_connectivityMaps_10000.mesh_proj.pts.RData")
dim(projected.points2) #32492     3
head(projected.points2)
#           [,1]       [,2]     [,3]
# [1,]  -4.496712 -43.769893 32.04211
# [2,] -18.425323 -41.241336 69.17256
# [3,] -52.468578  -6.555558 46.29697
# [4,]  -7.410121  25.434063 60.87795
# [5,] -17.897900 -91.522003 20.73351
# [6,] -48.365742 -52.942780 48.91669
sum(abs(obj2$locations - projected.points2)) #2841.204


