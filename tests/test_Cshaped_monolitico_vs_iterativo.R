if(!require(fdaPDEISCHIA)){
  devtools::install_github(repo ="aldoclemente/fdaPDEISCHIA")
}

if(!require(fdaPDEKIM)){
  devtools::install_github(repo ="aldoclemente/fdaPDEKIM")
}

library(fdaPDEISCHIA)

data(horseshoe2D)
mesh=create.mesh.2D(nodes=horseshoe2D$boundary_nodes, 
                    segments = horseshoe2D$boundary_segments)
mesh = refine.mesh.2D(mesh, maximum_area = 0.025, minimum_angle = 30)
locations <- mesh$nodes[!mesh$nodesmarkers, ]

mesh = refine.mesh.2D(mesh, maximum_area = 0.0075, minimum_angle = 30)
FEMbasis <- create.FEM.basis(mesh)

test.locations <- refine.mesh.2D(mesh, maximum_area = 0.0025, minimum_angle = 30)$nodes
plot(mesh)
points(locations, pch=16, col="blue")
points(test.locations, pch=16, col="red")

GCVFLAG=T
GCVMETHODFLAG='Exact' #For now, Exact because 'Stochastic' changes every time..

Cov1 <- function(x,y){
  sin(2*pi*x) * cos(2*pi*y)
}

fs.test.time<-function(x,y,t=y) {
  K <- (y/0.1*as.double((abs(y)<=0.1 & x>-0.5))+as.double((abs(y)>0.1 | x<=-0.5)))^2
  
  res=numeric(length =length(x))
  
  for(i in 1:length(x)) {
    if(x[i]>=0 && y[i]>0)
      res[i]=cos(t[i])*(0.25*pi+x[i])+(y[i]-0.5)^2
    
    if(x[i]>=0 && y[i]<=0)
      res[i]=cos(2*t[i])*(-0.25*pi-x[i])+(-y[i]-0.5)^2
    
    if(x[i]<0 && y[i]>0)
      res[i]=cos(t[i])*(-atan(y[i]/x[i])*0.5)+(sqrt(x[i]^2+y[i]^2)-0.5)^2*K[i]
    
    if(x[i]<0 && y[i]<=0)
      res[i]=cos(2*t[i])*(-atan(y[i]/x[i])*0.5)+(sqrt(x[i]^2+y[i]^2)-0.5)^2*K[i]
  }
  res
}

nlocs <- nrow(locations)
nnodes <- nrow(mesh$nodes)
betas <- as.matrix(c(3,0.5))
b <- as.matrix( c(-5, 0, 5) )


X1 = cbind( Cov1(locations[,1], locations[,2]), rnorm(nlocs, mean = 0, sd = 2))
X2 = cbind( Cov1(locations[,1], locations[,2]), rnorm(nlocs, mean = 0, sd = 2))
X3 = cbind( Cov1(locations[,1], locations[,2]), rnorm(nlocs, mean = 0, sd = 2))

# image(FEM(Cov1(mesh$nodes[,1], mesh$nodes[,2]), FEMbasis))
# image(FEM(rnorm(nnodes, sd=2), FEMbasis))
# 
# image(FEM(fs.test.time(x=mesh$nodes[,1], y=mesh$nodes[,2],  t=rep(0.5, nnodes)), FEMbasis))
# image(FEM(fs.test.time(x=mesh$nodes[,1], y=mesh$nodes[,2],  t=rep(2, nnodes)), FEMbasis))
# image(FEM(fs.test.time(x=mesh$nodes[,1], y=mesh$nodes[,2],  t=rep(1.5, nnodes)), FEMbasis))

func1 = fs.test.time(x=locations[,1], y=locations[,2],  t=rep(0.5, nlocs))
func2 = fs.test.time(x=locations[,1], y=locations[,2],  t=rep(1, nlocs))
func3 = fs.test.time(x=locations[,1], y=locations[,2],  t=rep(1.5, nlocs))
test.func1 = fs.test.time(x=test.locations[,1], y=test.locations[,2],  t=rep(0.5, nrow(test.locations)))
test.func2 = fs.test.time(x=test.locations[,1], y=test.locations[,2],  t=rep(1, nrow(test.locations)))
test.func3 = fs.test.time(x=test.locations[,1], y=test.locations[,2],  t=rep(1.5, nrow(test.locations)))

n_sim <- 30
results_MONOLITICO <- list(
  beta_1 = matrix(0, nrow=n_sim,ncol=1),
  beta_2 = matrix(0, nrow=n_sim,ncol=1),
  b_1 = matrix(0, nrow=n_sim,ncol=2), # (b11, b12)
  b_2 = matrix(0, nrow=n_sim,ncol=2), # (b21, b22)
  b_3 = matrix(0, nrow=n_sim,ncol=2),  # (b31, b32)
  time = matrix(0, nrow=n_sim,ncol=1)
)

results_ITERATIVO <- results_MONOLITICO

errors_MONOLITICO <- list(
  beta_1 = matrix(0, nrow=n_sim,ncol=1),
  beta_2 = matrix(0, nrow=n_sim,ncol=1),
  b_1 = matrix(0, nrow=n_sim,ncol=2),
  b_2 = matrix(0, nrow=n_sim,ncol=2),
  b_3 = matrix(0, nrow=n_sim,ncol=2),
  f_1 = matrix(0, nrow=n_sim,ncol=1),
  f_2 = matrix(0, nrow=n_sim,ncol=1),
  f_3 = matrix(0, nrow=n_sim,ncol=1),
  response = matrix(0, nrow=n_sim,ncol=1)
)

errors_ITERATIVO <- errors_MONOLITICO

rmse <- function(x,y){
  return(sqrt( mean( x-y )^2 ) ) 
}

for(i in 1:n_sim){
  # V == Cov1
  exact1 <- X1%*% betas + func1 + X1[,1] * b[1] 
  obs1 <- exact1  # + rnorm(nlocs,mean=0,sd=0.05*diff(range(exact1)))
  
  exact2 <- X2%*% betas + func2 + X2[,1] * b[2]
  obs2 <- exact2  # + rnorm(nlocs,mean=0,sd=0.05*diff(range(exact2)))
  
  exact3 <- X3%*% betas + func3 + X3[,1] * b[3] 
  obs3 <- exact3  # + rnorm(nlocs,mean=0,sd=0.05*diff(range(exact3)))
  
  observations <- c(obs1, obs2, obs3)
  observations <- observations + rnorm(nlocs*3, mean=0, sd=0.05*(diff(range(c(func1, func2, func3)))))
  observations <- matrix(observations, nrow=nlocs, ncol=3)
  X = rbind(X1, X2, X3)
  
  #lambda= 10^seq(-2,1,by=0.1) # 31
  lambda= 10^seq(-1,1,length=10)
  start <- Sys.time()
  output_MONOLITICO = smooth.FEM.mixed(observations = observations, locations = locations,
                                       covariates = X, random_effect = c(1,2),
                                       FEMbasis = FEMbasis, lambda = lambda, 
                                       lambda.selection.criterion = "grid", 
                                       lambda.selection.lossfunction = "GCV",
                                       DOF.evaluation = "exact")
  results_MONOLITICO$time[i] <- difftime(Sys.time(), start, units="secs")
  
  start <- Sys.time()
  output_ITERATIVO = smooth.FEM.mixed(observations = observations, locations = locations,
                                      covariates = X, random_effect = c(1,2),
                                      FEMbasis = FEMbasis, lambda = lambda, 
                                      lambda.selection.criterion = "grid", 
                                      lambda.selection.lossfunction = "GCV",
                                      DOF.evaluation = "exact", FLAG_ITERATIVE = TRUE)
  results_ITERATIVO$time[i] <- difftime(Sys.time(), start, units="secs")
  
  best_lambda_MONOLITICO <- output_MONOLITICO$bestlambda
  best_lambda_ITERATIVO <- output_ITERATIVO$bestlambda
  
  #results KIM
  results_MONOLITICO$beta_1[i] <- output_MONOLITICO$beta[1, best_lambda_MONOLITICO]
  results_MONOLITICO$beta_2[i] <- output_MONOLITICO$beta[2, best_lambda_MONOLITICO]
  results_MONOLITICO$b_1[i,] <- output_MONOLITICO$b_i[1:2, best_lambda_MONOLITICO]
  results_MONOLITICO$b_2[i,] <- output_MONOLITICO$b_i[3:4, best_lambda_MONOLITICO]
  results_MONOLITICO$b_3[i,] <- output_MONOLITICO$b_i[5:6, best_lambda_MONOLITICO]
  
  #results ISCHIA
  results_ITERATIVO$beta_1[i] <- output_ITERATIVO$beta[1, best_lambda_ITERATIVO]
  results_ITERATIVO$beta_2[i] <- output_ITERATIVO$beta[2, best_lambda_ITERATIVO]
  results_ITERATIVO$b_1[i,] <- output_ITERATIVO$b_i[1:2, best_lambda_ITERATIVO]
  results_ITERATIVO$b_2[i,] <- output_ITERATIVO$b_i[3:4, best_lambda_ITERATIVO]
  results_ITERATIVO$b_3[i,] <- output_ITERATIVO$b_i[5:6, best_lambda_ITERATIVO]
  
  # errors KIM
  errors_MONOLITICO$beta_1[i] <- rmse(output_MONOLITICO$beta[1, best_lambda_MONOLITICO], betas[1])
  errors_MONOLITICO$beta_2[i] <- rmse(output_MONOLITICO$beta[2, best_lambda_MONOLITICO], betas[2])
  errors_MONOLITICO$b_1[i,] <- c(rmse(results_MONOLITICO$b_1[i,1], b[1]), rmse(results_MONOLITICO$b_1[1,2], 0.))
  errors_MONOLITICO$b_2[i,] <- c(rmse(results_MONOLITICO$b_2[i,1], b[2]), rmse(results_MONOLITICO$b_2[1,2], 0.))
  errors_MONOLITICO$b_3[i,] <- c(rmse(results_MONOLITICO$b_3[i,1], b[3]), rmse(results_MONOLITICO$b_3[1,2], 0.))
  
  errors_MONOLITICO$f_1[i] <- rmse(eval.FEM(FEM(as.matrix(output_MONOLITICO$fit.FEM.mixed$coeff[1:nnodes,best_lambda_MONOLITICO]), FEMbasis), test.locations),
                            test.func1)
  errors_MONOLITICO$f_2[i] <- rmse(eval.FEM(FEM(as.matrix(output_MONOLITICO$fit.FEM.mixed$coeff[(nnodes+1):(2*nnodes),best_lambda_MONOLITICO]), FEMbasis), test.locations),
                            test.func2)
  errors_MONOLITICO$f_3[i] <- rmse(eval.FEM(FEM(as.matrix(output_MONOLITICO$fit.FEM.mixed$coeff[(2*nnodes+1):(3*nnodes),best_lambda_MONOLITICO]), FEMbasis), test.locations),
                            test.func3)
  
  y_hat1 <- X1%*% as.matrix(output_MONOLITICO$beta[, best_lambda_MONOLITICO]) + 
    eval.FEM(FEM(as.matrix(output_MONOLITICO$fit.FEM.mixed$coeff[1:nnodes,best_lambda_MONOLITICO]), FEMbasis), locations) + 
    X1%*% as.matrix(output_MONOLITICO$b_i[1:2, best_lambda_MONOLITICO])
  
  y_hat2 <- X2%*% as.matrix(output_MONOLITICO$beta[, best_lambda_MONOLITICO]) + 
    eval.FEM(FEM(as.matrix(output_MONOLITICO$fit.FEM.mixed$coeff[(nnodes+1):(2*nnodes),best_lambda_MONOLITICO]), FEMbasis), locations) + 
    X2%*%as.matrix(output_MONOLITICO$b_i[3:4, best_lambda_MONOLITICO])
  
  y_hat3 <- X3%*% as.matrix(output_MONOLITICO$beta[, best_lambda_MONOLITICO]) + 
    eval.FEM(FEM(as.matrix(output_MONOLITICO$fit.FEM.mixed$coeff[(2*nnodes+1):(3*nnodes),best_lambda_MONOLITICO]), FEMbasis), locations) + 
    X3%*% as.matrix(output_MONOLITICO$b_i[5:6, best_lambda_MONOLITICO])
  
  errors_MONOLITICO$response[i] <- rmse(c(y_hat1, y_hat2, y_hat3), as.vector(observations))
  
  # errors ISCHIA
  errors_ITERATIVO$beta_1[i] <- rmse(output_ITERATIVO$beta[1, best_lambda_ITERATIVO], betas[1])
  errors_ITERATIVO$beta_2[i] <- rmse(output_ITERATIVO$beta[2, best_lambda_ITERATIVO], betas[2])
  errors_ITERATIVO$b_1[i,] <- c(rmse(results_ITERATIVO$b_1[i,1], b[1]), rmse(results_ITERATIVO$b_1[1,2], 0.))
  errors_ITERATIVO$b_2[i,] <- c(rmse(results_ITERATIVO$b_2[i,1], b[2]), rmse(results_ITERATIVO$b_2[1,2], 0.))
  errors_ITERATIVO$b_3[i,] <- c(rmse(results_ITERATIVO$b_3[i,1], b[3]), rmse(results_ITERATIVO$b_3[1,2], 0.))
  
  errors_ITERATIVO$f_1[i] <- rmse(eval.FEM(FEM(as.matrix(output_ITERATIVO$fit.FEM.mixed$coeff[1:nnodes,best_lambda_ITERATIVO]), FEMbasis), test.locations),
                               test.func1)
  errors_ITERATIVO$f_2[i] <- rmse(eval.FEM(FEM(as.matrix(output_ITERATIVO$fit.FEM.mixed$coeff[(nnodes+1):(2*nnodes),best_lambda_ITERATIVO]), FEMbasis), test.locations),
                               test.func2)
  errors_ITERATIVO$f_3[i] <- rmse(eval.FEM(FEM(as.matrix(output_ITERATIVO$fit.FEM.mixed$coeff[(2*nnodes+1):(3*nnodes),best_lambda_ITERATIVO]), FEMbasis), test.locations),
                               test.func3)
  
  y_hat1 <- X1%*% as.matrix(output_ITERATIVO$beta[, best_lambda_ITERATIVO]) + 
    eval.FEM(FEM(as.matrix(output_ITERATIVO$fit.FEM.mixed$coeff[1:nnodes,best_lambda_ITERATIVO]), FEMbasis), locations) + 
    X1%*% as.matrix(output_ITERATIVO$b_i[1:2, best_lambda_ITERATIVO])
  
  y_hat2 <- X2%*% as.matrix(output_ITERATIVO$beta[, best_lambda_ITERATIVO]) + 
    eval.FEM(FEM(as.matrix(output_ITERATIVO$fit.FEM.mixed$coeff[(nnodes+1):(2*nnodes),best_lambda_ITERATIVO]), FEMbasis), locations) + 
    X2%*%as.matrix(output_ITERATIVO$b_i[3:4, best_lambda_ITERATIVO])
  
  y_hat3 <- X3%*% as.matrix(output_ITERATIVO$beta[, best_lambda_ITERATIVO]) + 
    eval.FEM(FEM(as.matrix(output_ITERATIVO$fit.FEM.mixed$coeff[(2*nnodes+1):(3*nnodes),best_lambda_ITERATIVO]), FEMbasis), locations) + 
    X3%*% as.matrix(output_ITERATIVO$b_i[5:6, best_lambda_ITERATIVO])
  
  errors_ITERATIVO$response[i] <- rmse(c(y_hat1, y_hat2, y_hat3), as.vector(observations))
}


errors <- results_MONOLITICO
errors$beta_1 <- errors$beta_1 - results_ITERATIVO$beta_1
errors$beta_2 <- errors$beta_2 - results_ITERATIVO$beta_2
errors$b_1 <- errors$b_1 - results_ITERATIVO$b_1
errors$b_2 <- errors$b_2 - results_ITERATIVO$b_2
errors$b_3 <- errors$b_3 - results_ITERATIVO$b_3

# Building folders -------------------------------------------------------------
date_ = unlist(strsplit(as.character(gsub(":","_",gsub(" ","-",Sys.time()))), 
                        split = "[.]"))[1]
if(!dir.exists("data/")) {
  dir.create("data/")
}

if( !dir.exists("data/test_Cshaped2D_mono_vs_iter/")){
  dir.create("data/test_Cshaped2D_mono_vs_iter/")
}

folder.name = paste("data/test_Cshaped2D_mono_vs_iter/",date_,"/",sep="")

if(!dir.exists(folder.name)) {
  dir.create(folder.name)
}

save(errors_ITERATIVO , errors_MONOLITICO, results_ITERATIVO, results_MONOLITICO, errors,
     file = paste0(folder.name, "data.RData"))

# ------------------------------------------------------------------------------
pdf(paste0(folder.name, "results.pdf"))
boxplot(data.frame(beta_1 = errors$beta_1, beta_2 = errors$beta_2), main="Errore\nbeta")
abline(h=0., lty=2, col="red")

boxplot(data.frame(b_11 = errors$b_1[,1], b_12 = errors$b_1[,2],
                   b_21 = errors$b_2[,1], b_22 = errors$b_2[,2],
                   b_31 = errors$b_3[,1], b_32 = errors$b_3[,2]), main="Errore\nb")
abline(h=0., lty=2, col="red")

boxplot(data.frame(f_1_monolitico = errors_MONOLITICO$f_1, f_1_iterativo = errors_ITERATIVO$f_1,
                   f_2_monolitico = errors_MONOLITICO$f_2, f_2_iterativo = errors_ITERATIVO$f_2, 
                   f_3_monolitico = errors_MONOLITICO$f_3, f_3_iterativo = errors_ITERATIVO$f_3),
        main="RMSE\nf")

boxplot(data.frame(beta_1_monolitico = results_MONOLITICO$beta_1, beta_1_iterativo = results_ITERATIVO$beta_1))
abline(h=betas[1], lty=2, col="red", main =" Beta 1")

boxplot(data.frame(beta_2_monolitico = results_MONOLITICO$beta_2, beta_2_iterativo = results_ITERATIVO$beta_2))
abline(h=betas[2], lty=2, col="red", main =" Beta 2")

boxplot(data.frame(monolitico=results_MONOLITICO$time,
                   iterativo = results_ITERATIVO$time), main="Time [s]")
abline(h=0, lty=2, col="red")
dev.off()


source("utils.R")

png("Cshaped_mesh.png")
plot(mesh, pch=".", asp=1)
dev.off()

png("Cshaped_locs.png")
plot(mesh, pch=".", asp=1)
points(locations, pch=16, col="red3")
dev.off()

fig <- contour.FEM(FEM(Cov1(mesh$nodes[,1], mesh$nodes[,2]), FEMbasis))
fig %>% layout( scene = list(
                  camera = list(
  eye = list(x = 0, y = -0.01,  z = 2.25))))

fig <- contour.FEM(FEM(rnorm(nnodes, sd=2), FEMbasis))
fig %>% layout( scene = list(
  camera = list(
    eye = list(x = 0, y = -0.01,  z = 2.25))))


fig <-contour.FEM(FEM(fs.test.time(x=mesh$nodes[,1], y=mesh$nodes[,2],  t=rep(0.5, nnodes)), FEMbasis))
fig %>% layout( scene = list(
  camera = list(
    eye = list(x = 0, y = -0.01,  z = 2.25))))

fig <- contour.FEM(FEM(fs.test.time(x=mesh$nodes[,1], y=mesh$nodes[,2],  t=rep(2, nnodes)), FEMbasis))
fig %>% layout( scene = list(
  camera = list(
    eye = list(x = 0, y = -0.01,  z = 2.25))))


fig <-contour.FEM(FEM(fs.test.time(x=mesh$nodes[,1], y=mesh$nodes[,2],  t=rep(1.5, nnodes)), FEMbasis))
fig %>% layout( scene = list(
  camera = list(
    eye = list(x = 0, y = -0.01,  z = 2.25))))

