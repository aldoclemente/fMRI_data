if(!require(pacman)) install.packages("pacman")
pacman::p_load("fdaPDEmixed", "R.matlab", "signal")
source("utils.R")
# read disk -----------------------------------------------------------------
system("whoami > user.txt")
user <- read.table("user.txt", header = FALSE)[,1]
unlink("user.txt")
ext_disk <- paste0("/media/", user, "/EXTERNAL_USB/")

# create mesh ------------------------------------------------------------------
nodes <- readMat(paste0(ext_disk, "data/vertices.mat"))
faces <- readMat(paste0(ext_disk, "data/faces.mat"))
nodes <- nodes$vertices
faces <- faces$faces
mesh <- create.mesh.2.5D(nodes = nodes, triangles = faces)

# brain regions (!!!) ----------------------------------------------------------
#parcellation <- "unknown"
#brain_regions = read.csv(paste0(ext_disk, "data/BrainRegions.csv"))
#dim(brain_regions) #64984    69

parcellation <- "gordon2016" # see Lila et al. (2016)
brain_regions = read.csv("data/gordon_parcellation2016.csv")
if(parcellation == "gordon2016")
  roi_idx = as.logical(brain_regions$Parcel_1) # precuneo

# read data --------------------------------------------------------------------
rfMRI_path <- paste0(ext_disk, "data/rfMRI_REST_processed/")
m <- length(list.files(rfMRI_path)) # numero soggetti 80..

m = 30
corr_map <- matrix(nrow=nrow(mesh$nodes), ncol=m)
filter_signal = FALSE
# read BOLD signal -------------------------------------------------------------
i = 1
for(file in list.files(rfMRI_path)[1:m]){
  fMRI_signal <- read.csv(paste0(rfMRI_path, file),
                   nrows=32492, header = FALSE, sep=' ')[,1:1200]
  corr_map[,i] <- correlation_maps(fMRI_signal, roi_idx, filter_signal)
  i = i + 1
}

if(!dir.exists("data/")) dir.create("data/")
save(corr_map, file = ifelse(filter_signal, 
                             paste0("data/", parcellation,"_corr_map_filtered.RData"),
                             paste0("data/", parcellation,"_corr_map.RData")))

FCmaps <- fisher.r2z(corr_map)
save(FCmaps, file = ifelse( filter_signal, 
                            paste0("data/", parcellation,"_FCmaps_filtered.RData"),
                            paste0("data/", parcellation,"_FCmaps.RData")))

# read thickness ---------------------------------------------------------------
thickness <- matrix(nrow=nrow(mesh$nodes), ncol=m)
thickness_path <- paste0(ext_disk, "data/thickness_processed/")
i = 1
for(file in list.files(thickness_path)[1:m]){
  thick <- read.csv(paste0(thickness_path, file),
                        nrows=32492, header = FALSE, sep=' ')[,1]
  thickness[,i] <- thick
  i = i + 1
}

save(thickness, file = paste0("data/", parcellation,"_thickness.RData"))

# FCmaps plots -------------------------------------------------------------------
folder.name <- paste0("data/", parcellation,"_FCmaps_filtered/")
if(!dir.exists(folder.name))
 dir.create(folder.name)

nnodes <- nrow(mesh$nodes)

min.col = min(FCmaps, na.rm = T)
max.col = max(FCmaps, na.rm = T)

FEMbasis <- create.FEM.basis(mesh)
for(i in 1:m){ # (nnodes+1):(2*nnodes)

  FEMobject <- FEM(coeff = FCmaps[,i],
                    FEMbasis = FEMbasis)
  
  plot(FEMobject, m=min.col, M=max.col)
  #rgl.postscript("brain_mesh.pdf", fmt="pdf")
  snapshot3d(filename = paste0(folder.name,"patient_", i,".png"),
             fmt = "png", width = 800, height = 750, webshot = rgl.useNULL())
  rgl.close()
}
# ------------------------------------------------------------------------------

na_idx <- vector(mode="list", length = m)
for(i in 1:m){
  na_idx[[i]] <- which(is.na(corr_map[,i]))
}

for(i in 1:m){
cat( length(na_idx[[i]]), "\n")
}

# FC maps ----------------------------------------------------------------------
folder.name <- paste0("data/", parcellation,"_FCmaps/")
if(!dir.exists(folder.name))
  dir.create(folder.name)

min.col = min(FCmaps, na.rm = T)
max.col = max(FCmaps, na.rm = T)

FEMbasis <- create.FEM.basis(mesh)
for(i in 1:m){ # (nnodes+1):(2*nnodes)
  
  FEMobject <- FEM(coeff = FCmaps[,i],
                   FEMbasis = FEMbasis)
  
  plot(FEMobject, m=min.col, M=max.col)
  #rgl.postscript("brain_mesh.pdf", fmt="pdf")
  snapshot3d(filename = paste0(folder.name,"patient_", i,".png"),
             fmt = "png", width = 800, height = 750, webshot = rgl.useNULL())
  rgl.close()
}
