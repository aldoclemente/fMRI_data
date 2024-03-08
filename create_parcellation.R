
#install.packages("ciftiTools")
library(ciftiTools)
ciftiTools.setOption("wb_path", val = "/usr/bin/wb_command")
ciftiTools.listOptions()

# 

parc <- read_cifti("/home/aldo/Desktop/Parcels_release/Parcels_LR.dlabel.nii")
plot(parc)

summary(parc)
tmp <- parc$data$cortex_left


labels <- as.data.frame(parc$meta$cifti$labels ) 

dim(labels)
table(labels$X..Alpha)
RGB <- la
head(labels)

cortex_left <- parc$data$cortex_left

gordon_parcellation <- matrix(0, nrow=nrow(cortex_left), ncol=length(rownames(labels)))

for(i in 1:nrow(cortex_left)){
  gordon_parcellation[i, cortex_left[i] + 1] = 1
}

gordon_parcellation <- as.data.frame(gordon_parcellation)
colnames(gordon_parcellation) <- c("na", rownames(labels)[-1])
write.csv(gordon_parcellation, file="data/gordon_parcellation2016.csv")

head(gordon_parcellation)

# -----
source("utils.R")
pacman::p_load("fdaPDEmixed", "R.matlab", "signal")
nodes <- readMat("data/vertices.mat")
faces <- readMat("data/faces.mat")
nodes <- nodes$vertices
faces <- faces$faces
mesh <- create.mesh.2.5D(nodes = nodes, triangles = faces)

na_ <- as.logical(gordon_parcellation$na)
plot.mesh.2.5D(mesh, NA_ = na_ )

roi <- as.logical(gordon_parcellation$Parcel_1) # precuneo !
plot.mesh.2.5D(mesh, ROI = roi, NA_=na_)
snapshot3d(filename = "brain_mesh.png",
           fmt = "png", width = 800, height = 750, webshot = rgl.useNULL())
rgl.close()

head(labels)
