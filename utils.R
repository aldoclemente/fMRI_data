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

plot.mesh.2.5D <- function(mesh, M = NULL, m = NULL, ...){
  
  FEM = FEM(rep(0, nrow(mesh$nodes)), create.FEM.basis(mesh))
  
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
    col = (col - min(coeff[,isurf]))/diffrange*(ncolor-1)+1
    
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
