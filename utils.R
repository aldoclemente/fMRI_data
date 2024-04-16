contour.FEM <- function(x, limits=NULL, ...){
  mesh <- x$FEMbasis$mesh
  plot_data <- data.frame(  X=mesh$nodes[,1], 
                            Y=mesh$nodes[,2],
                            Z=rep(0, times=nrow(mesh$nodes)),
                            coeff=x$coeff[1:nrow(mesh$nodes)])
    
  coeff <- apply(mesh$triangles, MARGIN=1, FUN = function(edge){
    mean(x$coeff[edge,])
  })
  
  if(is.null(limits)) limits = c(min(coeff), max(coeff))
  cmin = limits[1]; cmax=limits[2]
  
  I=mesh$triangles[,1]-1
  J=mesh$triangles[,2]-1
  K=mesh$triangles[,3]-1
  fig<- plot_ly(plot_data, x=~X, y=~Y, z=~Z,
                i = I, j = J, k = K, cmin = limits[1], cmax=limits[2],
                intensity=~coeff, color=~coeff, type="mesh3d", 
                colorbar=list(title=""), ...) %>%
    layout(scene = list(
      aspectmode = "data", 
      xaxis = list(
        title = '', showgrid = F, zeroline = F, showticklabels = F),
      yaxis = list(
        title = '', showgrid = F, zeroline = F, showticklabels = F),
      zaxis = list(
        title = '', showgrid = F, zeroline = F, showticklabels = F),
      camera = list(
        eye = list(x = 0, y = -0.01,  z = 1.25))),  dragmode="zoom") %>%
    colorbar(len = 1, title="")
}

# }else{
#   plot_data <- data.frame(X=rep(mesh$nodes()[,1], times=length(times)), 
#                           Y=rep(mesh$nodes()[,2], times=length(times)),
#                           Z=rep(0, times=length(times)),
#                           coeff=as.vector(x$coefficients()[1:nrow(mesh$nodes()),]),
#                           times = round(rep(times, each=nrow(mesh$nodes())),3))
#   coeff <- matrix(nrow=nrow(mesh$elements()), 
#                   ncol=length(times))
#   for(t in 1:length(times)){
#     coeff[,t] <- apply(mesh$elements(), MARGIN=1, FUN=
#                          function(edge){
#                            mean(x$coefficients()[edge,t])
#                          })
#   }
#   limits = c(min(coeff), max(coeff))
#   cmin = limits[1]; cmax=limits[2]
#   I=mesh$elements()[,1]-1
#   J=mesh$elements()[,2]-1
#   K=mesh$elements()[,3]-1
#   fig<- plot_ly(plot_data, x=~X, y=~Y, z=~Z, frame=~times,
#                 i = I, j = J, k = K, cmin = limits[1], cmax=limits[2],
#                 intensity=~coeff, color=~coeff, type="mesh3d",
#                 colorbar=list(title=""), ...) %>%
#     layout(scene = list(
#       aspectmode = "data", 
#       xaxis = list(
#         title = '', showgrid = F, zeroline = F, showticklabels = F),
#       yaxis = list(
#         title = '', showgrid = F, zeroline = F, showticklabels = F),
#       zaxis = list(
#         title = '', showgrid = F, zeroline = F, showticklabels = F),
#       camera = list(
#         eye = list(x = 0, y = -0.01,  z = 1.25))), dragmode="zoom") %>%
#     colorbar(len = 1, title="") %>%
#     animation_slider(currentvalue = list(prefix ="t = "))
#   #animation_opts(frame=5) %>% 
# }
# fig 

plot_mesh.2.5D <- function(mesh, ...){
  data_plot <- data.frame(X=mesh$nodes[,1], 
                          Y=mesh$nodes[,2],
                          Z=mesh$nodes[,3])
  
  data_edge <- data.frame(X= mesh$nodes[mesh$edges[,1]][,1],
                          Y= mesh$nodes[mesh$edges[,1]][,2],
                          Z= mesh$nodes[mesh$edges[,1]][,3])
  
  I=(mesh$triangles[,1]-1); J=(mesh$triangles[,2]-1); K=(mesh$triangles[,3]-1)
  fig <- plot_ly(data_plot,
          type = 'mesh3d', x = ~X, y = ~Y, z = ~Z,
          i = I, j = J, k = K,
          hoverinfo = 'none', facecolor="lightgray")%>%

    layout(scene = list(
      aspectmode = "data", 
      xaxis = list(title = '',showgrid = F,zeroline = F,showticklabels = F),
      yaxis = list(title = '',showgrid = F,zeroline = F,showticklabels = F),
      zaxis = list(title = '',showgrid = F,zeroline = F,showticklabels = F)))
}

# ------------------------------------------------------------------------------

## GRAPHICAL SETTINGS ----------------------------------------------------------
zoom = 0.7657689 # 0.6
userMatrix = rbind(c(  0.02091786,  0.99873853, -0.04564825,    0),
                   c( -0.13139695,  0.04800860,  0.99016660,    0),
                   c(  0.99110913, -0.01471432,  0.13223548,    0),
                   c(  0.00000000,  0.0000000,   0.0000000,     1))
windowRect = c(70,  106, 1920, 1117)

#pp <- par3d(no.readonly = TRUE)

plot.mesh.2.5D <- function(mesh, M = NULL, m = NULL, ROI=NULL, NA_ = NULL,...){
  
  FEM = FEM(rep(0, nrow(mesh$nodes)), create.FEM.basis(mesh))
  FEM$coeff[ROI,] <- 1 
  FEM$coeff[NA_,] <- 2
  
  
  if (is.null(m)) { m = min(FEM$coeff)}
  if (is.null(M)) { M = max(FEM$coeff)}
  triangles = c(t(FEM$FEMbasis$mesh$triangles))
  ntriangles = nrow(FEM$FEMbasis$mesh$triangles)
  order = FEM$FEMbasis$mesh$order
  nodes = FEM$FEMbasis$mesh$nodes
  edges = matrix(rep(0, 6*ntriangles), ncol = 2)
  for(i in 0:(ntriangles-1)){
    edges[3*i+1,] = c(triangles[3*order*i+1], triangles[3*order*i+2])
    edges[3*i+2,] = c(triangles[3*order*i+1], triangles[3*order*i+3])
    edges[3*i+3,] = c(triangles[3*order*i+2], triangles[3*order*i+3])
  }
  edges = edges[!duplicated(edges),]
  edges <- as.vector(t(edges))
  
  coeff = FEM$coeff
  
  FEMbasis = FEM$FEMbasis
  
  mesh = FEMbasis$mesh
 
  #p = jet.col(n = 1000, alpha = 0.8)
  # alternative color palette: p <- colorRampPalette(c("#0E1E44", "#3E6DD8", "#68D061", "#ECAF53", "#EB5F5F", "#E11F1C"))(1000)
  p = c("lightgray", "red3", "blue3")
  palette(p)
  
  ncolor = length(p)
  nsurf = dim(coeff)[[2]]
  for (isurf in 1:nsurf)
  {
    open3d(zoom = zoom, userMatrix = userMatrix, windowRect = windowRect)
    rgl.pop("lights") 
    light3d(specular = "black") 
    
    diffrange = M - m
    
    col = coeff[triangles,isurf]
    col = (col - min(coeff[,isurf]))/diffrange*(ncolor-1)+1
    
    rgl.triangles(x = nodes[triangles ,1], y = nodes[triangles ,2],
                  z = nodes[triangles,3],
                  color = col,...)
    rgl.lines(x = nodes[edges ,1], y = nodes[edges ,2],
              z = nodes[edges,3],
              color = "black",...)
    aspect3d("iso")
    
    if (nsurf > 1 && isurf<nsurf)
    {readline("Press a button for the next plot...")}
  }
}

plot.FEM <- function(FEM, M = NULL, m = NULL, ...){
  
  if (is.null(m)) { m = min(FEM$coeff, na.rm = T)}
  if (is.null(M)) { M = max(FEM$coeff, na.rm = T)}
  triangles = c(t(FEM$FEMbasis$mesh$triangles))
  ntriangles = nrow(FEM$FEMbasis$mesh$triangles)
  order = FEM$FEMbasis$mesh$order
  nodes = FEM$FEMbasis$mesh$nodes
  edges = matrix(rep(0, 6*ntriangles), ncol = 2)
  for(i in 0:(ntriangles-1)){
    edges[3*i+1,] = c(triangles[3*order*i+1], triangles[3*order*i+2])
    edges[3*i+2,] = c(triangles[3*order*i+1], triangles[3*order*i+3])
    edges[3*i+3,] = c(triangles[3*order*i+2], triangles[3*order*i+3])
  }
  edges = edges[!duplicated(edges),]
  edges <- as.vector(t(edges))
  
  coeff = FEM$coeff
  
  FEMbasis = FEM$FEMbasis
  
  mesh = FEMbasis$mesh
  
  p = jet.col(n = 1000, alpha = 0.8)
  # alternative color palette: p <- colorRampPalette(c("#0E1E44", "#3E6DD8", "#68D061", "#ECAF53", "#EB5F5F", "#E11F1C"))(1000)
  
  palette(p)
  
  ncolor = length(p)
  
  nsurf = dim(coeff)[[2]]
  for (isurf in 1:nsurf)
  {
    open3d(zoom = zoom, userMatrix = userMatrix, windowRect = windowRect)
    rgl.pop("lights") 
    light3d(specular = "black") 
    
    diffrange = M - m
    
    col = coeff[triangles,isurf]
    col = (col - min(coeff[,isurf], na.rm =T))/diffrange*(ncolor-1)+1
    
    rgl.triangles(x = nodes[triangles ,1], y = nodes[triangles ,2],
                  z = nodes[triangles,3],
                  color = col,...)
    # rgl.lines(x = nodes[edges ,1], y = nodes[edges ,2],
    #           z = nodes[edges,3],
    #           color = "black",...)
    aspect3d("iso")
    
    if (nsurf > 1 && isurf<nsurf)
    {readline("Press a button for the next plot...")}
  }
}

# Preprocessing ----------------------------------------------------------------
filter_time_series <- function(Y, interval, filter_signal){
  Y_filtered <- Y
  if(filter_signal){
    nrows <- nrow(Y)
    ncols <- ncol(Y)
    
    Y_filtered <- matrix(nrow=nrows, ncol=ncols)
    bf <- signal::butter(3, interval, type="pass")
    
    Y_filtered <-  t(apply(Y, MARGIN=1, FUN=function(row){
                          
            signal::filter(bf, row)
    }) )
  }
  
  Y_filtered
}

correlation_maps <- function(fMRI_signal, roi_idx, filter_signal = FALSE){
  
  time_series = filter_time_series(fMRI_signal, 
                                   interval = c(0.009, 0.08), filter_signal = filter_signal)
  roi_time_series = time_series[roi_idx,]
  roi_mean_time_series = colMeans(roi_time_series) # == ncol(fmRI_signal) OK!
  
  corr_map <- apply(time_series, MARGIN=1, function(row){
    cor(row, roi_mean_time_series)
  })
  corr_map
}

fisher.r2z <- function(r) { 0.5 * (log(1+r) - log(1-r)) }

# FCmap_func<-function(corr) {
#   
#   time_series=t(as.matrix(fMRI))
#   ROI_ts=time_series[,ROI_idx]
#   # print(dim(ROI_ts)) #1200 1469
#   
#   cross_sec_avg_ROI_row=as.vector(rowMeans(ROI_ts))
#   # print(length(cross_sec_avg_ROI_row)) #1200
#   cor_nodes<-NULL
#   for(j in 1:ncol(time_series)){
#     cor_nodes<-c(cor_nodes,cor(cross_sec_avg_ROI_row, time_series[,j])) #have 1200 length in each nodes
#   }
#   # print(length(cor_nodes)) #32492
#   print(length(which(is.na(cor_nodes)))) #2796
#   
#   r = cor_nodes
#   z = 0.5 * log((1+r)/(1-r)) #fisher transformation
#   z
# }