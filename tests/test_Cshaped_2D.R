library(fdaPDEKIM)

data(horseshoe2D)
mesh=create.mesh.2D(nodes=horseshoe2D$boundary_nodes, 
                    segments = horseshoe2D$boundary_segments)
mesh = refine.mesh.2D(mesh, maximum_area = 0.025, minimum_angle = 30)
locations <- mesh$nodes[!mesh$nodesmarkers, ]

mesh = refine.mesh.2D(mesh, maximum_area = 0.0075, minimum_angle = 30)
FEMbasis <- create.FEM.basis(mesh)

plot(mesh)
points(locations, pch=16, col="blue")
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
Cov2 <- rnorm(nlocs, mean = 0, sd = 2)

X1 = cbind( Cov1(locations[,1], locations[,2]), rnorm(nlocs, mean = 0, sd = 2))
X2 = cbind( Cov1(locations[,1], locations[,2]), rnorm(nlocs, mean = 0, sd = 2))
X3 = cbind( Cov1(locations[,1], locations[,2]), rnorm(nlocs, mean = 0, sd = 2))

image(FEM(Cov1(mesh$nodes[,1], mesh$nodes[,2]), FEMbasis))
image(FEM(rnorm(nnodes, sd=2), FEMbasis))

image(FEM(fs.test.time(x=mesh$nodes[,1], y=mesh$nodes[,2],  t=rep(0.5, nnodes)), FEMbasis))
image(FEM(fs.test.time(x=mesh$nodes[,1], y=mesh$nodes[,2],  t=rep(2, nnodes)), FEMbasis))
image(FEM(fs.test.time(x=mesh$nodes[,1], y=mesh$nodes[,2],  t=rep(1.5, nnodes)), FEMbasis))

func1 = fs.test.time(x=locations[,1], y=locations[,2],  t=rep(0.5, nlocs))
func2 = fs.test.time(x=locations[,1], y=locations[,2],  t=rep(1, nlocs))
func3 = fs.test.time(x=locations[,1], y=locations[,2],  t=rep(1.5, nlocs))

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

lambda= 10^seq(-2,1,by=0.1) # 31
lambda= 10^seq(-1,1,length=10)
output_CPP = fdaPDEKIM:::smooth.FEM.mixed(observations = observations, locations = locations,
                                           covariates = X, random_effect = c(1,2),
                                           FEMbasis = FEMbasis, lambda = lambda,
                                           GCV=GCVFLAG, GCVmethod = GCVMETHODFLAG)

output_CPP$beta
output_CPP$b_i
plot(lambda, output_CPP$GCV)
best_lambda <- output_CPP$bestlambda
lambda[best_lambda]

plot(FEM(output_CPP$fit.FEM.mixed$coeff[1:nnodes,best_lambda], FEMbasis))
plot(FEM(fs.test.time(x=mesh$nodes[,1], y=mesh$nodes[,2],  t=rep(0.5, nnodes)), FEMbasis))
y_hat1 <- X1%*% as.matrix(output_CPP$beta[, best_lambda]) + 
            eval.FEM(FEM(as.matrix(output_CPP$fit.FEM.mixed$coeff[1:nnodes,best_lambda]), FEMbasis), locations) + 
              X1%*% as.matrix(output_CPP$b_i[1:2, best_lambda])

y_hat2 <- X2%*% as.matrix(output_CPP$beta[, best_lambda]) + 
  eval.FEM(FEM(as.matrix(output_CPP$fit.FEM.mixed$coeff[(nnodes+1):(2*nnodes),best_lambda]), FEMbasis), locations) + 
  X2%*%as.matrix(output_CPP$b_i[3:4, best_lambda])

y_hat3 <- X3%*% as.matrix(output_CPP$beta[, best_lambda]) + 
  eval.FEM(FEM(as.matrix(output_CPP$fit.FEM.mixed$coeff[(2*nnodes+1):(3*nnodes),best_lambda]), FEMbasis), locations) + 
  X3%*% as.matrix(output_CPP$b_i[5:6, best_lambda])

rmse_1 <- sqrt( mean((y_hat1 - exact1)^2) )
rmse_1

rmse_2 <- sqrt( mean((y_hat2 - exact2)^2) )
rmse_2

rmse_3 <- sqrt( mean((y_hat3 - exact3)^2) )
rmse_3

evaluations <- fdaPDEKIM:::eval.FEM.mixed(output_CPP$fit.FEM.mixed, locations = locations)
dim(evaluations)
