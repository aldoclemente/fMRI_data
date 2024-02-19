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