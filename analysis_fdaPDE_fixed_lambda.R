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

#load("data/data.RData") 

load("data/FCmaps.RData")
load("data/thickness.RData")
# fixed lambda -----------------------------------------------------------------
lambda =  0.1
m = 30
observations <- FCmaps[,1:m]
#idx_na <-unique(which(is.na(observations),arr.ind = TRUE)[,1]) # delle locazioni che andranno eliminate
#locations <- mesh$nodes

#observations <- observations[-idx_na,]
#covariates <- X[-idx_na,1:m]
#locations <- locations[-idx_na,]

covariates <- thickness[,1:m]
dim(observations)
dim(covariates)

start_ <- Sys.time()
output_mixed <- smooth.FEM.mixed(
  observations = as.matrix(observations), #locations = locations,
  covariates = as.vector(covariates), random_effect = c(1), 
  lambda = lambda,
  FEMbasis = FEMbasis, FLAG_ITERATIVE = TRUE)
time_ <- difftime(Sys.time(), start_, units ="mins")
output_mixed$beta

date_ = unlist(strsplit(as.character(gsub(":","_",gsub(" ","-",Sys.time()))), 
                        split = "[.]"))[1]

folder.name <- paste0(folder.name, date_, "/")
if(!dir.exists(folder.name))
  dir.create(folder.name)
save(output_mixed, file = paste0(folder.name, date_, ".RData"))

# ------------------------------------------------------------------------------
library(fdaPDEmixed)
source("utils.R")
folder.name <- "data/application/2024-03-06-16_56_52/"
load(paste0(folder.name, "2024-03-06-16_56_52.RData"))
m = 30
mesh <- output_mixed$fit.FEM.mixed$FEMbasis$mesh
FEMbasis <- output_mixed$fit.FEM.mixed$FEMbasis

nnodes <- nrow(mesh$nodes)
NA_mask <- which(is.na(observations[,1]))

for(i in 1:m){
  coeff <- output_mixed$fit.FEM.mixed$coeff[((i-1)*nnodes+1):(i*nnodes),]
  coeff[NA_mask] <- NA
  output_mixed$fit.FEM.mixed$coeff[((i-1)*nnodes+1):(i*nnodes),] <- coeff
}

min.col = min(output_mixed$fit.FEM.mixed$coeff, na.rm = T)
max.col = max(output_mixed$fit.FEM.mixed$coeff, na.rm = T)
if(!dir.exists(paste0(folder.name, "estimates/")))
  dir.create(paste0(folder.name, "estimates/"))
for(i in 1:m){ # (nnodes+1):(2*nnodes)
  FEMobject <- FEM(coeff = output_mixed$fit.FEM.mixed$coeff[((i-1)*nnodes+1):(i*nnodes),] ,
                   FEMbasis = FEMbasis)
  
  plot(FEMobject, m=min.col, M=max.col)
  #rgl.postscript("brain_mesh.pdf", fmt="pdf")
  snapshot3d(filename = paste0(folder.name,"estimates/","patient_", i,".png"),
             fmt = "png", width = 800, height = 750, webshot = rgl.useNULL())
  rgl.close()
}

if(!dir.exists(paste0(folder.name, "fc_maps/")))
  dir.create(paste0(folder.name, "fc_maps/"))
for(i in 1:m){ # (nnodes+1):(2*nnodes)
  
  FEMobject <- FEM(coeff = FCmaps[,i],
                   FEMbasis = FEMbasis)
  
  plot(FEMobject, m=min.col, M=max.col)
  #rgl.postscript("brain_mesh.pdf", fmt="pdf")
  snapshot3d(filename = paste0(folder.name,"fc_maps/","patient_", i,".png"),
             fmt = "png", width = 800, height = 750, webshot = rgl.useNULL())
  rgl.close()
}
