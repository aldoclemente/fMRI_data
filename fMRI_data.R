# download data ----------------------------------------------------------------
if(!require(pacman)) install.packages("pacman")
pacman::p_load("neurohcp", "dplyr", "fdaPDE", "R.matlab")

ACCESS_KEY = '' 
SECRET_KEY = ''

system("ntpq -p")
set_aws_api_key(access_key = ACCESS_KEY, 
                secret_key = SECRET_KEY, 
                default_region = "eu-south-1")
neurohcp::bucketlist(verbose = FALSE)


### total patient id #####
idxs = hcp_900_scanning_info %>%
  filter(scan_type %in% "dMRI") %>%
  select(id) %>%
  unique

dim(idxs) #488   1 (510   1 without filter)
head(idxs)

idxs <- as.integer(idxs$id)
## fMRI ------------------------------------------------------------------------
ext_disk <- "/media/aldo/EXTERNAL_USB/"

if(!dir.exists(paste0(ext_disk, "data/")))
  dir.create(paste0(ext_disk, "data/"))

dest_path_rest <- paste0(ext_disk, "data/rfMRI_REST/")
dest_path_thick <- paste0(ext_disk, "data/thickness/")

if(!dir.exists(dest_path_rest))
  dir.create(dest_path_rest)

if(!dir.exists(dest_path_thick))
  dir.create(dest_path_thick)

for(i in 1:length(idxs[1:80])){
  
  input_path <- paste0("HCP/",idxs[i], "/MNINonLinear/")
  download_hcp_file(paste0(input_path, "Results/rfMRI_REST1_LR/rfMRI_REST1_LR_Atlas.dtseries.nii"),
                    destfile = paste0(dest_path_rest, idxs[i],".rfMRI_REST1_LR_Atlas.dtseries.nii"),
                    verbose = FALSE)
}

mask_no_thickness <- c()
for(i in 1:length(idxs[1:80])){
  input_path <- paste0("HCP/",idxs[i], "/MNINonLinear/")
  tryCatch({
    download_hcp_file(paste0(input_path, "fsaverage_LR32k/", idxs[i], ".thickness.32k_fs_LR.dscalar.nii"),
                      destfile = paste0(dest_path_thick, idxs[i],".thickness.32k_fs_LR.dscalar.nii"),
                      verbose = FALSE)
  }, error = function(e){
    mask_no_thickness <- append(mask_no_thickness, i)
  })
}

# processing -------------------------------------------------------------------
# run from bash processing_data.sh 

# open atlas data --------------------------------------------------------------
# assuming mesh is the 2.5D mesh of the left hemisphere of the brain (Conte69 atlas) 
ext_disk <- "/media/aldoclemente/EXTERNAL_USB/"
nodes <- readMat(paste0(ext_disk, "data/vertices.mat"))
faces <- readMat(paste0(ext_disk, "data/faces.mat"))

nodes <- nodes$vertices
faces <- faces$faces

mesh <- create.mesh.2.5D(nodes = nodes, triangles = faces)
plot(mesh)
FEMbasis <- create.FEM.basis(mesh)

# FC maps ----------------------------------------------------------------------
brain_regions = read.csv(paste0(ext_disk, "data/BrainRegions.csv"))
dim(brain_regions) #64984    69

ROI_idx = brain_regions$L_precuneus[1:(nrow(mesh$nodes))] #nrow(mesh$nodes): 32492 ????
ROI_idx = as.logical(ROI_idx)
length(which(ROI_idx == TRUE)) #1469

### FC map function ####
FCmap_func<-function(fMRI) {
  
  time_series=t(as.matrix(fMRI))
  ROI_ts=time_series[,ROI_idx]
  # print(dim(ROI_ts)) #1200 1469
  
  cross_sec_avg_ROI_row=as.vector(rowMeans(ROI_ts))
  # print(length(cross_sec_avg_ROI_row)) #1200
  cor_nodes<-NULL
  for(j in 1:ncol(time_series)){
    cor_nodes<-c(cor_nodes,cor(cross_sec_avg_ROI_row, time_series[,j])) #have 1200 length in each nodes
  }
  # print(length(cor_nodes)) #32492
  print(length(which(is.na(cor_nodes)))) #2796
  
  r = cor_nodes
  z = 0.5 * log((1+r)/(1-r)) #fisher transformation
  z
}

# open 1D data -----------------------------------------------------------------
# thickness <- read.csv(paste0(ext_disk, "data/thickness_processed/100307.thickness.32k_fs_LR.dscalar.1D"), 
#                       nrows=32492, header = FALSE, sep=' ')
# head(thickness)
# dim(thickness)
# thickness[1,1]
# thickness[32492,1]
# 
# 
# fMRI <- read.csv(paste0(ext_disk, "data/rfMRI_REST_processed/100307.rfMRI_REST1_LR_Atlas.dtseries.1D"),
#                           nrows=32492, header = FALSE, sep=' ')[,1:1200]
#dim(fMRI)

rfMRI_path <- paste0(ext_disk, "data/rfMRI_REST_processed/")
n <- length(list.files(rfMRI_path)) # numero soggetti 

FC_maps <- matrix(nrow=nrow(mesh$nodes), ncol=0)
for(file in list.files(rfMRI_path)){
  fMRI <- read.csv(paste0(rfMRI_path, file),
                   nrows=32492, header = FALSE, sep=' ')[,1:1200]
  FC_maps <- cbind(FC_maps, FCmap_func(fMRI))
}

X <- matrix(nrow=nrow(mesh$nodes), ncol=0)
thickness_path <- paste0(ext_disk, "data/thickness_processed/")
for(file in list.files(thickness_path)){
  thickness <- read.csv(paste0(thickness_path, file),
                         nrows=32492, header = FALSE, sep=' ')[,1]
  X <- cbind(X, thickness)
}

save(FC_maps, X, file = "data/data.RData")