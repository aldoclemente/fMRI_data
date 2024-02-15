#####WITH random-effect ####

library(fdaPDEKIM)
library(mgcv)

data(horseshoe2D)
mesh=create.mesh.2D(nodes=horseshoe2D$boundary_nodes, 
                    segments = horseshoe2D$boundary_segments)
mesh=refine.mesh.2D(mesh, maximum_area = 0.025, minimum_angle = 30)
FEMbasis <- create.FEM.basis(mesh)
plot(mesh)

GCVFLAG=T
GCVMETHODFLAG='Exact' #For now, Exact because 'Stochastic' changes every time..
nnodes = dim(mesh$nodes)[1]

### lambda####
before_lambda= 10^seq(-2,1,by=0.1) #len: 31

### boundary ####
bnd = horseshoe2D$boundary_nodes
bound <- list(list(x = bnd[,1], y = bnd[,2])) #, f = rep(0, nrow(crds))))
N <- 10
gx <- seq(min(bnd[,1]), max(bnd[,1]), len = N)
gy <- seq(min(bnd[,2]), max(bnd[,2]), len = N)
gp <- expand.grid(gx, gy)
names(gp) <- c("x","y")
knots <- data.frame(v=rep(seq(-.5,3,by=.5),4),
                    w=rep(c(-.6,-.3,.3,.6),rep(8,4)))
names(knots) <- c("x", "y")
names(bound[[1]]) <- c("x", "y") #, "f")
bound[[1]]$y[c(10:28)] = bound[[1]]$y[c(10:28)] +0.005
bound[[1]]$y[c(82:100)] = bound[[1]]$y[c(82:100)] -0.005
bound[[1]]$x[c(1:9,101:108)] = bound[[1]]$x[c(1:9,101:108)]-0.005
bound[[1]]$x[c(29:35,75:81)] = bound[[1]]$x[c(29:35,75:81)]+0.005
bound[[1]]$y[c(36:54)] = bound[[1]]$y[c(36:54)]-0.005
bound[[1]]$x[c(55)]  = bound[[1]]$x[c(55)]+0.005
bound[[1]]$y[c(56:74)] = bound[[1]]$y[c(56:74)]+0.005


# Exact data - Nodes locations (for C shape)
fs.test.time<-function(x,y,t) {
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


## locations ####
output_CPP1 = smooth.FEM(observations = sin(pi*mesh$nodes[,1]), #observations at mesh nodes
                         FEMbasis = FEMbasis,
                         lambda = before_lambda,
                         GCV=GCVFLAG,
                         GCVmethod = GCVMETHODFLAG)

#### uniform Grid Initialization
gridnum = 20
minV1 = min(mesh$nodes[,1])
maxV1 = max(mesh$nodes[,1])
minV2 = min(mesh$nodes[,2])
maxV2 = max(mesh$nodes[,2])

x = seq(minV1+0.1, maxV1-0.1, length.out = gridnum) ### SHOULD ADJUST THIS! (OR ELSE, ERROR IN SOAP FILM)
y= seq(minV2+0.1, maxV2-0.1, length.out = gridnum)
unif_grid <- expand.grid(x = x, y = y)

points1=eval.FEM(output_CPP1$fit.FEM,
                 locations = unif_grid)

loc = unif_grid[-which(is.na(points1)),]

# randomize
set.seed(5847947)
ind = sample.int(dim(loc)[1])[1:100] ##100 locations points as default
loc <- loc[ind,]
nlocs = dim(loc)[1]


### covariates ####
cov1=sin(2*pi*loc[,1])*cos(2*pi*loc[,2])
cov2=rnorm(nlocs, mean=0, sd=2) #previous sd=0.1 was too small for inference...
W=cbind(cov1,cov2)

mesh_cov1=sin(2*pi*mesh$nodes[,1])*cos(2*pi*mesh$nodes[,2])
mesh_cov2=rnorm(nnodes, mean=0, sd=2)
mesh_W=cbind(mesh_cov1, mesh_cov2)

# Fix betas
beta_exact1=c(3-5, 0.5)
beta_exact2=c(3, 0.5)
beta_exact3=c(3+5, 0.5)

mean(5*cov1 + 5*cov2) 

num_units = 3
RMSE_func<-function(f,g) {
  sqrt(mean((f-g)^2))
}


func_evaluation1 = numeric(nlocs)
func_evaluation1 =  fs.test.time(x=loc[,1], y=loc[,2],  t=rep(0.5, nlocs))

func_evaluation2 = numeric(nlocs)
func_evaluation2 = fs.test.time(x=loc[,1], y=loc[,2],  t=rep(1, nlocs))

func_evaluation3 = numeric(nlocs)
func_evaluation3 = fs.test.time(x=loc[,1], y=loc[,2],  t=rep(1.5, nlocs))
  

func1 = numeric(nnodes)
func1 = fs.test.time(x=mesh$nodes[,1], y=mesh$nodes[,2],  t=rep(0.5, nnodes))

func2 = numeric(nnodes)
func2 = fs.test.time(x=mesh$nodes[,1], y=mesh$nodes[,2],  t=rep(1, nnodes))

func3 = numeric(nnodes)
func3 = fs.test.time(x=mesh$nodes[,1], y=mesh$nodes[,2],  t=rep(1.5, nnodes))
  

# try to set as similar range
range(func_evaluation1)[2] - range(func_evaluation1)[1]
range(W%*%beta_exact1)[2] - range(W%*%beta_exact1)[1]

range(func_evaluation2)[2] - range(func_evaluation2)[1]
range(W%*%beta_exact2)[2] - range(W%*%beta_exact2)[1]

range(func_evaluation3)[2] - range(func_evaluation3)[1]
range(W%*%beta_exact3)[2] - range(W%*%beta_exact3)[1]


DF <- data.frame(W%*%beta_exact1, func_evaluation1,
                 W%*%beta_exact2, func_evaluation2,
                 W%*%beta_exact3, func_evaluation3)
x11()
# png('range_cov_func.png')
# par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
boxplot(DF, col = c('red','blue'), at = c(1:2,4:5,7:8), xaxt = "n", main='range',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
axis(side = 1, at = c(1.5,4.5,7.5), labels = c("unit1","unit2","unit3"),
     cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
legend("topleft", fill = c('red','blue'), legend = c('covariate','function'), horiz = F,
       pt.cex=1.5, cex=1.5)
# legend("topright", inset=c(-0.3,0), fill = c('red','blue'), legend = c('covariate','function'), horiz = F,
#        pt.cex=1.5, cex=1.5)

mixed_loc = rbind(loc, loc, loc)
mixed_cov1 = c(cov1, cov1, cov1)
mixed_cov2 = c(cov2, cov2, cov2)
# mixed_obs = c(W%*%beta_exact1+func_evaluation1, 
#               W%*%beta_exact2+func_evaluation2,
#               W%*%beta_exact3+func_evaluation3)
mixed_obs = cbind(W%*%beta_exact1+func_evaluation1,
                   W%*%beta_exact2+func_evaluation2,
                   W%*%beta_exact3+func_evaluation3)
covariates=rbind(W,W,W)
mixed_nlocs= nlocs*3

#### range ####
ran=range(c(func_evaluation1, func_evaluation2, func_evaluation3))

ran[2]-ran[1]
0.05*(ran[2]-ran[1])

# Set number of simulation trials
N=50

simul_observations_before <- array(rep(0, N*nrow(loc)*num_units), 
                                   dim=c(N, nrow(loc), num_units))

RMSE_func<-function(f,g) {
  sqrt(mean((f-g)^2))
}


mod_info = fdaPDEKIM:::smooth.FEM.mixed(locations = loc,
                            observations = mixed_obs, #obs doesn't effect tree, bary, dof
                            covariates = covariates,
                            random_effect = c(1),
                            FEMbasis = FEMbasis,
                            lambda = before_lambda,
                            GCV = GCVFLAG)


mixed_before_betamat=matrix(data=NA, nrow = 2, ncol=N)
mixed_before_rmse_beta1 =NULL
mixed_before_rmse_beta2 =NULL
mixed_before_bimat=matrix(data=NA, nrow = 3, ncol=N)
mixed_before_rmse_b1_1 =NULL
mixed_before_rmse_b2_1 =NULL
mixed_before_rmse_b3_1 =NULL
mixed_before_RMSE_f1 = NULL
mixed_before_RMSE_f2 = NULL
mixed_before_RMSE_f3 = NULL
mixed_before_global_RMSE1 = NULL
mixed_before_global_RMSE2 = NULL
mixed_before_global_RMSE3 = NULL
mixed_before_RMSE =NULL
before_selected_lambda=rep(0,N)

#### before NA simulation ######
for(i in 1:N){
  set.seed(54425 + i)
  data = mixed_obs + rnorm(nlocs*3,mean=0,sd=0.05*(ran[2]-ran[1])) #add column-wise (convertable with matrix and vector)
  simul_observations_before[i,,]=data
  
  output_CPP<-smooth.FEM.mixed(locations = loc,
                                observations=data, 
                                covariates = covariates, #covariates added!
                                random_effect = c(1),
                                # FEMbasis=FEMbasis, 
                                FEMbasis = mod_info$fit.FEM.mixed$FEMbasis, #reuse tree info
                                lambda=before_lambda,
                                GCV=GCVFLAG,
                                GCVmethod = GCVMETHODFLAG,
                                bary.locations = mod_info$bary.locations,
                                DOF_matrix = mod_info$edf) #reuse dof info
  
  before_selected_lambda[i]=output_CPP$bestlambda

  # beta estimates
  mixed_before_betamat[1,i]=output_CPP$beta[,output_CPP$bestlambda][1]
  mixed_before_betamat[2,i]=output_CPP$beta[,output_CPP$bestlambda][2]
  
  # bi estimates
  mixed_before_bimat[1,i]=output_CPP$b_i[,output_CPP$bestlambda][1]
  mixed_before_bimat[2,i]=output_CPP$b_i[,output_CPP$bestlambda][2]
  mixed_before_bimat[3,i]=output_CPP$b_i[,output_CPP$bestlambda][3]
  
  # beta RMSE
  RMSE = RMSE_func(beta_exact2[1],mixed_before_betamat[1,i])
  mixed_before_rmse_beta1=c(mixed_before_rmse_beta1,RMSE)
  
  RMSE = RMSE_func(beta_exact2[2],mixed_before_betamat[2,i])
  mixed_before_rmse_beta2=c(mixed_before_rmse_beta2,RMSE)
  
  # bi MSE
  RMSE = RMSE_func(-5,mixed_before_bimat[1,i])
  mixed_before_rmse_b1_1=c(mixed_before_rmse_b1_1,RMSE)
  
  RMSE = RMSE_func(0, mixed_before_bimat[2,i])
  mixed_before_rmse_b2_1=c(mixed_before_rmse_b2_1,RMSE)
  
  RMSE = RMSE_func(5,mixed_before_bimat[3,i])
  mixed_before_rmse_b3_1=c(mixed_before_rmse_b3_1,RMSE)
  
  # f MSE
  fitted.func = eval.FEM.mixed(output_CPP$fit.FEM.mixed,
                               locations = loc)[,output_CPP$bestlambda]
  fitted.func1 = fitted.func[1:nlocs]
  fitted.func2 = fitted.func[(nlocs+1):(2*nlocs)]
  fitted.func3 = fitted.func[(2*nlocs+1):(3*nlocs)]
  
  RMSE=RMSE_func(func_evaluation1,fitted.func1)
  mixed_before_RMSE_f1=c(mixed_before_RMSE_f1, RMSE)
  
  RMSE=RMSE_func(func_evaluation2,fitted.func2)
  mixed_before_RMSE_f2=c(mixed_before_RMSE_f2, RMSE)
  
  RMSE=RMSE_func(func_evaluation3,fitted.func3)
  mixed_before_RMSE_f3=c(mixed_before_RMSE_f3, RMSE)
  
  #global MSE
  RMSE=RMSE_func(W%*%beta_exact1 + func_evaluation1,
           W%*%mixed_before_betamat[,i] + cov1*mixed_before_bimat[1,i] + fitted.func1)
  mixed_before_global_RMSE1=c(mixed_before_global_RMSE1, RMSE)
  
  RMSE=RMSE_func(W%*%beta_exact2 + func_evaluation2,
           W%*%mixed_before_betamat[,i] + cov1*mixed_before_bimat[2,i] + fitted.func2)
  mixed_before_global_RMSE2=c(mixed_before_global_RMSE2, RMSE)
  
  RMSE=RMSE_func(W%*%beta_exact3 + func_evaluation3,
           W%*%mixed_before_betamat[,i] + cov1*mixed_before_bimat[3,i] + fitted.func3)
  mixed_before_global_RMSE3=c(mixed_before_global_RMSE3, RMSE)
  
  # total
  value = c(W%*%beta_exact1 + func_evaluation1, 
            W%*%beta_exact2 + func_evaluation2, 
            W%*%beta_exact3 + func_evaluation3)
  fitted.value = c(W%*%mixed_before_betamat[,i] + cov1*mixed_before_bimat[1,i] + fitted.func1,
                   W%*%mixed_before_betamat[,i] + cov1*mixed_before_bimat[2,i] + fitted.func2,
                   W%*%mixed_before_betamat[,i] + cov1*mixed_before_bimat[3,i] + fitted.func3)
  mixed_before_RMSE=c(mixed_before_RMSE, RMSE_func(value, fitted.value))
}

before_selected_lambda


### with NA #####
mixed_after_betamat=matrix(data=NA, nrow = 2, ncol=N)
mixed_after_rmse_beta1 =NULL
mixed_after_rmse_beta2 =NULL
mixed_after_bimat=matrix(data=NA, nrow = 3, ncol=N)
mixed_after_rmse_b1_1 =NULL
mixed_after_rmse_b2_1 =NULL
mixed_after_rmse_b3_1 =NULL
mixed_after_RMSE_f1 = NULL
mixed_after_RMSE_f2 = NULL
mixed_after_RMSE_f3 = NULL
mixed_after_global_RMSE1 = NULL
mixed_after_global_RMSE2 = NULL
mixed_after_global_RMSE3 = NULL
mixed_after_RMSE =NULL
after_selected_lambda=rep(0,N)

### num of NA ####
keep_orig = simul_observations_before
simul_observations_before = keep_orig

tot_len = length(mixed_obs)
ratio = 0.30
NA_len = tot_len*ratio
num_NA = NA_len/num_units # num of NA in each statistical unit
NA_len/tot_len
num_NA

### after NA lambda ####
### after NA, lambda changes.... needed to be fixed accordingly...
after_lambda= seq(2,8, length.out=30) #len: 30

length(which(is.na(simul_observations_before[i,,1])))
length(which(is.na(simul_observations_before[i,,2])))
length(which(is.na(simul_observations_before[i,,3])))

### after NA simulation #####
for(i in 1:N){
  for (j in 1:num_units) {
    set.seed(452185+i+j)
    NAIndx=sample(1:nlocs, num_NA)
    simul_observations_before[i,NAIndx,j] = NA
  }

  output_CPP2<-fdaPDEKIM:::smooth.FEM.mixed(locations = loc,
                                observations=simul_observations_before[i,,], 
                                covariates = covariates, #covariates added!
                                random_effect = c(1),
                                # FEMbasis=FEMbasis,
                                FEMbasis = mod_info$fit.FEM$FEMbasis, #reuse tree info
                                lambda=after_lambda,
                                GCV=GCVFLAG,
                                GCVmethod = GCVMETHODFLAG,
                                bary.locations = mod_info$bary.locations)
  
  after_selected_lambda[i]=output_CPP2$bestlambda
  
  # beta estimates
  mixed_after_betamat[1,i]=output_CPP2$beta[,output_CPP2$bestlambda][1]
  mixed_after_betamat[2,i]=output_CPP2$beta[,output_CPP2$bestlambda][2]
  
  # bi estimates
  mixed_after_bimat[1,i]=output_CPP2$b_i[,output_CPP2$bestlambda][1]
  mixed_after_bimat[2,i]=output_CPP2$b_i[,output_CPP2$bestlambda][2]
  mixed_after_bimat[3,i]=output_CPP2$b_i[,output_CPP2$bestlambda][3]
  
  # beta RMSE
  RMSE = RMSE_func(beta_exact2[1],mixed_after_betamat[1,i])
  mixed_after_rmse_beta1=c(mixed_after_rmse_beta1,RMSE)
  
  RMSE = RMSE_func(beta_exact2[2],mixed_after_betamat[2,i])
  mixed_after_rmse_beta2=c(mixed_after_rmse_beta2,RMSE)
  
  # bi MSE
  RMSE = RMSE_func(-5, mixed_after_bimat[1,i])
  mixed_after_rmse_b1_1=c(mixed_after_rmse_b1_1,RMSE)
  
  RMSE = RMSE_func(0, mixed_after_bimat[2,i])
  mixed_after_rmse_b2_1=c(mixed_after_rmse_b2_1,RMSE)
  
  RMSE = RMSE_func(5, mixed_after_bimat[3,i])
  mixed_after_rmse_b3_1=c(mixed_after_rmse_b3_1,RMSE)
  
  # f MSE
  fitted.func = eval.FEM.mixed(output_CPP2$fit.FEM.mixed,
                               locations = loc)[,output_CPP2$bestlambda]
  
  fitted.func1 = fitted.func[1:nlocs]
  fitted.func2 = fitted.func[(nlocs+1):(2*nlocs)]
  fitted.func3 = fitted.func[(2*nlocs+1):(3*nlocs)]
  
  RMSE=RMSE_func(func_evaluation1, fitted.func1)
  mixed_after_RMSE_f1=c(mixed_after_RMSE_f1, RMSE)
  
  RMSE=RMSE_func(func_evaluation2, fitted.func2)
  mixed_after_RMSE_f2=c(mixed_after_RMSE_f2, RMSE)
  
  RMSE=RMSE_func(func_evaluation3,fitted.func3)
  mixed_after_RMSE_f3=c(mixed_after_RMSE_f3, RMSE)
  
  
  #global MSE
  RMSE=RMSE_func(W%*%beta_exact1 + func_evaluation1,
           W%*%mixed_after_betamat[,i] + cov1*mixed_after_bimat[1,i] + fitted.func1)
  mixed_after_global_RMSE1=c(mixed_after_global_RMSE1, RMSE)
  
  RMSE=RMSE_func(W%*%beta_exact2 + func_evaluation2,
           W%*%mixed_after_betamat[,i] + cov1*mixed_after_bimat[2,i] + fitted.func2)
  mixed_after_global_RMSE2=c(mixed_after_global_RMSE2, RMSE)
  
  RMSE=RMSE_func(W%*%beta_exact3 + func_evaluation3,
           W%*%mixed_after_betamat[,i] + cov1*mixed_after_bimat[3,i] + fitted.func3)
  mixed_after_global_RMSE3=c(mixed_after_global_RMSE3, RMSE)
  
  # total
  value = c(W%*%beta_exact1 + func_evaluation1, 
            W%*%beta_exact2 + func_evaluation2, 
            W%*%beta_exact3 + func_evaluation3)
  fitted.value = c(W%*%mixed_after_betamat[,i] + cov1*mixed_after_bimat[1,i] + fitted.func1,
                   W%*%mixed_after_betamat[,i] + cov1*mixed_after_bimat[2,i] + fitted.func2,
                   W%*%mixed_after_betamat[,i] + cov1*mixed_after_bimat[3,i] + fitted.func3)
  mixed_after_RMSE=c(mixed_after_RMSE, RMSE_func(value, fitted.value))
}

before_lambda[before_selected_lambda]
after_lambda[after_selected_lambda]

##### boxplot ####
prefix = paste0('NA_', num_NA)
prefix = paste0(prefix,'_')

png(paste0(prefix,'mixed_RMSE.png'))
boxplot(mixed_before_RMSE, mixed_after_RMSE, names = c('No NA','30% NA'), main='High randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)


png(paste0(prefix,'beta1.png'))
boxplot(cbind(mixed_before_betamat[1,],mixed_after_betamat[1,]), names=c('No NA','30% NA'), main='High randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
abline(h=beta_exact2[1], col='red')

png(paste0(prefix,'beta2.png'))
boxplot(cbind(mixed_before_betamat[2,],mixed_after_betamat[2,]), names=c('No NA','30% NA'), main='High randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
abline(h=beta_exact2[2], col='red')


png(paste0(prefix,'b1_1.png'))
boxplot(cbind(mixed_before_bimat[1,],mixed_after_bimat[1,]), names=c('No NA','30% NA'), main='High randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
abline(h=-5, col='red')

png(paste0(prefix,'b2_1.png'))
boxplot(cbind(mixed_before_bimat[2,],mixed_after_bimat[2,]), names=c('No NA','30% NA'), main='High randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
abline(h=0, col='red')

png(paste0(prefix,'b3_1.png'))
boxplot(cbind(mixed_before_bimat[3,],mixed_after_bimat[3,]), names=c('No NA','30% NA'), main='High randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
abline(h=5, col='red')


png(paste0(prefix,'f1 RMSE.png'))
boxplot(cbind(mixed_before_RMSE_f1, mixed_after_RMSE_f1), names = c('No NA','30% NA'), main='High randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)

png(paste0(prefix,'f2 RMSE.png'))
boxplot(cbind(mixed_before_RMSE_f2, mixed_after_RMSE_f2), names = c('No NA','30% NA'), main='High randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)

png(paste0(prefix,'f3 RMSE.png'))
boxplot(cbind(mixed_before_RMSE_f3, mixed_after_RMSE_f3), names = c('No NA','30% NA'), main='High randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)

DF <- data.frame(mixed_before_RMSE_f1, mixed_after_RMSE_f1, 
                 mixed_before_RMSE_f2, mixed_after_RMSE_f2, 
                 mixed_before_RMSE_f3, mixed_after_RMSE_f3)
png(paste0(prefix,'f RMSE.png'))
boxplot(DF, col = c('red','blue'), at = c(1:2,4:5,7:8), xaxt = "n", main='High randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
axis(side = 1, at = c(1.5,4.5,7.5), labels = c("f1","f2","f3"),
     cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
title(ylab="RMSE", line=2.3, cex.lab=2)
legend("topleft", fill = c('red','blue'), legend = c('No NA','30% NA'), horiz = F,
       pt.cex=1.5, cex=1.5)

png(paste0(prefix,'beta1 RMSE.png'))
boxplot(cbind(mixed_before_rmse_beta1, mixed_after_rmse_beta1), names = c('No NA','30% NA'), main='High randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)

png(paste0(prefix,'beta2 RMSE.png'))
boxplot(cbind(mixed_before_rmse_beta2, mixed_after_rmse_beta2), names = c('No NA','30% NA'), main='High randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)

DF <- data.frame(mixed_before_rmse_beta1, mixed_after_rmse_beta1, 
                 mixed_before_rmse_beta2, mixed_after_rmse_beta2)
png(paste0(prefix,'beta RMSE.png'))
boxplot(DF, col = c('red','blue'), at = c(1:2,4:5), main='High randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, xaxt = "n")
axis(side = 1, at = c(1.5,4.5), labels = c("beta1","beta2"),
     cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
title(ylab="RMSE", line=2.3, cex.lab=2)
legend("topright", fill = c('red','blue'), legend = c('No NA','30% NA'), horiz = F,
       pt.cex=1.5, cex=1.5)


png(paste0(prefix,'b1_1 RMSE.png'))
boxplot(cbind(mixed_before_rmse_b1_1, mixed_after_rmse_b1_1), names = c('No NA','30% NA'), main='High randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)

png(paste0(prefix,'b2_1 RMSE.png'))
boxplot(cbind(mixed_before_rmse_b2_1, mixed_after_rmse_b2_1), names = c('No NA','30% NA'), main='High randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)

png(paste0(prefix,'b3_1 RMSE.png'))
boxplot(cbind(mixed_before_rmse_b3_1, mixed_after_rmse_b3_1), names = c('No NA','30% NA'), main='High randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)

DF <- data.frame(mixed_before_rmse_b1_1, mixed_after_rmse_b1_1, 
                 mixed_before_rmse_b2_1, mixed_after_rmse_b2_1, 
                 mixed_before_rmse_b3_1, mixed_after_rmse_b3_1)
png(paste0(prefix,'bi_1 RMSE.png'))
boxplot(DF, col = c('red','blue'), at = c(1:2,4:5,7:8), xaxt = "n", main='High randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
axis(side = 1, at = c(1.5,4.5,7.5), labels = c("b1_1","b2_1","b3_1"),
     cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
title(ylab="RMSE", line=2.3, cex.lab=2)
legend("topleft", fill = c('red','blue'), legend = c('No NA','30% NA'), horiz = F,
       pt.cex=1.5, cex=1.5)


png(paste0(prefix,'unit1 global RMSE.png'))
boxplot(cbind(mixed_before_global_RMSE1, mixed_after_global_RMSE1), names = c('No NA','30% NA'), main='High randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)

png(paste0(prefix,'unit2 global RMSE.png'))
boxplot(cbind(mixed_before_global_RMSE2, mixed_after_global_RMSE2), names = c('No NA','30% NA'), main='High randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)

png(paste0(prefix,'unit3 global RMSE.png'))
boxplot(cbind(mixed_before_global_RMSE3, mixed_after_global_RMSE3), names = c('No NA','30% NA'), main='High randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)


DF <- data.frame(mixed_before_global_RMSE1, mixed_after_global_RMSE1, 
                 mixed_before_global_RMSE2, mixed_after_global_RMSE2, 
                 mixed_before_global_RMSE3, mixed_after_global_RMSE3)
png(paste0(prefix,'unit global RMSE.png'))
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
boxplot(DF, col = c('red','blue'), at = c(1:2,4:5,7:8), xaxt = "n", main='High randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
axis(side = 1, at = c(1.5,4.5,7.5), labels = c("unit1","unit2","unit3"),
     cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
legend("topright", inset=c(-0.3,0), fill = c('red','blue'), legend = c('No NA','30% NA'), horiz = F,
       pt.cex=1, cex=1)

#### image ####
#### before NA function image #####
library(fdaPDE)
# library(fields)

mx = seq(min(mesh$nodes[,1]),max(mesh$nodes[,1]), length=500)
my = seq(min(mesh$nodes[,2]),max(mesh$nodes[,2]), length=500)
unif_grid = expand.grid(mx, my)
image_len = length(mx)

image_cov1=sin(2*pi*unif_grid[,1])*cos(2*pi*unif_grid[,2])
image_cov2=rnorm(dim(unif_grid)[1], mean=0, sd=2)
image_W=cbind(image_cov1, image_cov2)


eval_sol = eval.FEM.mixed(output_CPP$fit.FEM.mixed,
                          locations = unif_grid)[,output_CPP$bestlambda]

#### NEEED TO BE SAME FOR ALL! (CHECK AGAIN!)
min_before_func = floor(min(eval_sol,na.rm=TRUE))
max_before_func = ceiling(max(eval_sol,na.rm=TRUE))
# -3, 4

#### CHOSEN IN THE LAST MIN
func_min = -3
func_max = 4

num_units = output_CPP$fit.FEM.mixed$num_units
for (i in 1:num_units) {
  evalmat = matrix(eval_sol[((i-1)*image_len*image_len+1):(i*image_len*image_len)], 
                   nrow=image_len, 
                   ncol=image_len, byrow=F)
  
  title = 'before_NA_f'
  title =paste0(title,i)
  title = paste0(title,'.png')
  png(title)
  # image(mx, my, evalmat, axes = TRUE, xlab = 'x', ylab = 'y', col=hcl.colors(100, "YlOrRd", rev = TRUE))
  # contour(mx, my, evalmat, levels = seq(min, max, by=1), add = TRUE, labcex=1, col = "black")
  fields::image.plot(evalmat, zlim=c(func_min,func_max), col = hcl.colors(100, "YlOrRd", rev = TRUE))
  contour(evalmat, add=TRUE, levels = seq(func_min, func_max, by=1), labcex=1.5, col = "black")
}


#### after NA function image #####
eval_sol2 = eval.FEM.mixed(output_CPP2$fit.FEM.mixed,
                          locations = unif_grid)[,output_CPP2$bestlambda]

#### NEEED TO BE SAME FOR ALL! (CHECK AGAIN!)
min_after_func = floor(min(eval_sol2,na.rm=TRUE))
max_after_func = ceiling(max(eval_sol2,na.rm=TRUE))
# -2, 4


num_units = output_CPP2$fit.FEM.mixed$num_units
for (i in 1:num_units) {
  evalmat = matrix(eval_sol2[((i-1)*image_len*image_len+1):(i*image_len*image_len)], 
                   nrow=image_len, 
                   ncol=image_len, byrow=F)
  
  title = 'after_NA_f'
  title =paste0(title,i)
  title = paste0(title,'.png')
  png(title)
  # image(mx, my, evalmat, axes = TRUE, xlab = 'x', ylab = 'y', col=hcl.colors(100, "YlOrRd", rev = TRUE))
  # contour(mx, my, evalmat, levels = seq(min, max, by=1), add = TRUE, labcex=1, col = "black")
  fields::image.plot(evalmat, zlim=c(func_min,func_max), col = hcl.colors(100, "YlOrRd", rev = TRUE))
  contour(evalmat, add=TRUE, levels = seq(func_min, func_max, by=1), labcex=1.5, col = "black")
}

####.####
##### LESS random-effect ####
rm(list=ls())
setwd("C:/Users/JIYOUNG KIM/Desktop/tests - mixed")
library(fdaPDE)
library(mgcv)

data(horseshoe2D)
mesh=create.mesh.2D(nodes=boundary_nodes, segments = boundary_segments)
mesh=refine.mesh.2D(mesh, maximum_area = 0.025, minimum_angle = 30)
FEMbasis <- create.FEM.basis(mesh)

GCVFLAG=T
GCVMETHODFLAG='Exact' #For now, Exact because 'Stochastic' changes every time..
nnodes = dim(mesh$nodes)[1]

### lambda####
before_lambda= 10^seq(-2,1,by=0.1) #len: 31

### boundary ####
bnd = boundary_nodes
bound <- list(list(x = bnd[,1], y = bnd[,2])) #, f = rep(0, nrow(crds))))
N <- 10
gx <- seq(min(bnd[,1]), max(bnd[,1]), len = N)
gy <- seq(min(bnd[,2]), max(bnd[,2]), len = N)
gp <- expand.grid(gx, gy)
names(gp) <- c("x","y")
knots <- data.frame(v=rep(seq(-.5,3,by=.5),4),
                    w=rep(c(-.6,-.3,.3,.6),rep(8,4)))
names(knots) <- c("x", "y")
names(bound[[1]]) <- c("x", "y") #, "f")
bound[[1]]$y[c(10:28)] = bound[[1]]$y[c(10:28)] +0.005
bound[[1]]$y[c(82:100)] = bound[[1]]$y[c(82:100)] -0.005
bound[[1]]$x[c(1:9,101:108)] = bound[[1]]$x[c(1:9,101:108)]-0.005
bound[[1]]$x[c(29:35,75:81)] = bound[[1]]$x[c(29:35,75:81)]+0.005
bound[[1]]$y[c(36:54)] = bound[[1]]$y[c(36:54)]-0.005
bound[[1]]$x[c(55)]  = bound[[1]]$x[c(55)]+0.005
bound[[1]]$y[c(56:74)] = bound[[1]]$y[c(56:74)]+0.005


# Exact data - Nodes locations (for C shape)
fs.test.time<-function(x,y,t) {
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


## locations ####
output_CPP1 = smooth.FEM(observations = sin(pi*mesh$nodes[,1]), #observations at mesh nodes
                         FEMbasis = FEMbasis,
                         lambda = before_lambda,
                         GCV=GCVFLAG,
                         GCVmethod = GCVMETHODFLAG)

#### uniform Grid Initialization
gridnum = 20
minV1 = min(mesh$nodes[,1])
maxV1 = max(mesh$nodes[,1])
minV2 = min(mesh$nodes[,2])
maxV2 = max(mesh$nodes[,2])

x = seq(minV1+0.1, maxV1-0.1, length.out = gridnum) ### SHOULD ADJUST THIS! (OR ELSE, ERROR IN SOAP FILM)
y= seq(minV2+0.1, maxV2-0.1, length.out = gridnum)
unif_grid <- expand.grid(x = x, y = y)

points1=eval.FEM(output_CPP1$fit.FEM,
                 locations = unif_grid)

loc = unif_grid[-which(is.na(points1)),]

# randomize
set.seed(5847947)
ind = sample.int(dim(loc)[1])[1:100] ##100 locations points as default
loc <- loc[ind,]
nlocs = dim(loc)[1]


### covariates ####
cov1=sin(2*pi*loc[,1])*cos(2*pi*loc[,2])
cov2=rnorm(nlocs, mean=0, sd=2) #previous sd=0.1 was too small for inference...
W=cbind(cov1,cov2)

mesh_cov1=sin(2*pi*mesh$nodes[,1])*cos(2*pi*mesh$nodes[,2])
mesh_cov2=rnorm(nnodes, mean=0, sd=2)
mesh_W=cbind(mesh_cov1, mesh_cov2)

# Fix betas
beta_exact1=c(3-1, 0.5)
beta_exact2=c(3, 0.5)
beta_exact3=c(3+1, 0.5)

mean(5*cov1 + 5*cov2) 

num_units = 3
RMSE_func<-function(f,g) {
  sqrt(mean((f-g)^2))
}

func_evaluation1 = numeric(nlocs)
func_evaluation1 =  fs.test.time(x=loc[,1], y=loc[,2],  t=rep(0.5, nlocs))

func_evaluation2 = numeric(nlocs)
func_evaluation2 = fs.test.time(x=loc[,1], y=loc[,2],  t=rep(1, nlocs))

func_evaluation3 = numeric(nlocs)
func_evaluation3 = fs.test.time(x=loc[,1], y=loc[,2],  t=rep(1.5, nlocs))


func1 = numeric(nnodes)
func1 = fs.test.time(x=mesh$nodes[,1], y=mesh$nodes[,2],  t=rep(0.5, nnodes))

func2 = numeric(nnodes)
func2 = fs.test.time(x=mesh$nodes[,1], y=mesh$nodes[,2],  t=rep(1, nnodes))

func3 = numeric(nnodes)
func3 = fs.test.time(x=mesh$nodes[,1], y=mesh$nodes[,2],  t=rep(1.5, nnodes))


# try to set as similar range
range(func_evaluation1)[2] - range(func_evaluation1)[1]
range(W%*%beta_exact1)[2] - range(W%*%beta_exact1)[1]

range(func_evaluation2)[2] - range(func_evaluation2)[1]
range(W%*%beta_exact2)[2] - range(W%*%beta_exact2)[1]

range(func_evaluation3)[2] - range(func_evaluation3)[1]
range(W%*%beta_exact3)[2] - range(W%*%beta_exact3)[1]

DF <- data.frame(W%*%beta_exact1, func_evaluation1,
                 W%*%beta_exact2, func_evaluation2,
                 W%*%beta_exact3, func_evaluation3)
x11()
# png('range_cov_func.png')
# par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
boxplot(DF, col = c('red','blue'), at = c(1:2,4:5,7:8), xaxt = "n", main='range',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
axis(side = 1, at = c(1.5,4.5,7.5), labels = c("unit1","unit2","unit3"),
     cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
legend("topleft", fill = c('red','blue'), legend = c('covariate','function'), horiz = F,
       pt.cex=1.5, cex=1.5)
# legend("topright", inset=c(-0.3,0), fill = c('red','blue'), legend = c('covariate','function'), horiz = F,
#        pt.cex=1.5, cex=1.5)

mixed_loc = rbind(loc, loc, loc)
mixed_cov1 = c(cov1, cov1, cov1)
mixed_cov2 = c(cov2, cov2, cov2)
# mixed_obs = c(W%*%beta_exact1+func_evaluation1, 
#               W%*%beta_exact2+func_evaluation2,
#               W%*%beta_exact3+func_evaluation3)
mixed_obs = cbind(W%*%beta_exact1+func_evaluation1,
                  W%*%beta_exact2+func_evaluation2,
                  W%*%beta_exact3+func_evaluation3)
covariates=rbind(W,W,W)
mixed_nlocs= nlocs*3

#### range ####
ran=range(c(func_evaluation1, func_evaluation2, func_evaluation3))

ran[2]-ran[1]
0.05*(ran[2]-ran[1])

# Set number of simulation trials
N=50

simul_observations_before <- array(rep(0, N*nrow(loc)*num_units), 
                                   dim=c(N, nrow(loc), num_units))

RMSE_func<-function(f,g) {
  sqrt(mean((f-g)^2))
}


mod_info = fdaPDEKIM:::smooth.FEM.mixed(locations = loc,
                            observations = mixed_obs, #obs doesn't effect tree, bary, dof
                            covariates = covariates,
                            random_effect = c(1),
                            FEMbasis = FEMbasis,
                            lambda = before_lambda,
                            GCV = GCVFLAG)


mixed_before_betamat=matrix(data=NA, nrow = 2, ncol=N)
mixed_before_rmse_beta1 =NULL
mixed_before_rmse_beta2 =NULL
mixed_before_bimat=matrix(data=NA, nrow = 3, ncol=N)
mixed_before_rmse_b1_1 =NULL
mixed_before_rmse_b2_1 =NULL
mixed_before_rmse_b3_1 =NULL
mixed_before_RMSE_f1 = NULL
mixed_before_RMSE_f2 = NULL
mixed_before_RMSE_f3 = NULL
mixed_before_global_RMSE1 = NULL
mixed_before_global_RMSE2 = NULL
mixed_before_global_RMSE3 = NULL
mixed_before_RMSE =NULL
before_selected_lambda=rep(0,N)

#### before NA simulation ######
for(i in 1:N){
  set.seed(54425 + i)
  data = mixed_obs + rnorm(nlocs*3,mean=0,sd=0.05*(ran[2]-ran[1])) #add column-wise (convertable with matrix and vector)
  simul_observations_before[i,,]=data
  
  output_CPP<-fdaPDEKIM:::smooth.FEM.mixed(locations = loc,
                               observations=data, 
                               covariates = covariates, #covariates added!
                               random_effect = c(1),
                               # FEMbasis=FEMbasis, 
                               FEMbasis = mod_info$fit.FEM.mixed$FEMbasis, #reuse tree info
                               lambda=before_lambda,
                               GCV=GCVFLAG,
                               GCVmethod = GCVMETHODFLAG,
                               bary.locations = mod_info$bary.locations,
                               DOF_matrix = mod_info$edf) #reuse dof info
  
  before_selected_lambda[i]=output_CPP$bestlambda
  
  # beta estimates
  mixed_before_betamat[1,i]=output_CPP$beta[,output_CPP$bestlambda][1]
  mixed_before_betamat[2,i]=output_CPP$beta[,output_CPP$bestlambda][2]
  
  # bi estimates
  mixed_before_bimat[1,i]=output_CPP$b_i[,output_CPP$bestlambda][1]
  mixed_before_bimat[2,i]=output_CPP$b_i[,output_CPP$bestlambda][2]
  mixed_before_bimat[3,i]=output_CPP$b_i[,output_CPP$bestlambda][3]
  
  # beta RMSE
  RMSE = RMSE_func(beta_exact2[1],mixed_before_betamat[1,i])
  mixed_before_rmse_beta1=c(mixed_before_rmse_beta1,RMSE)
  
  RMSE = RMSE_func(beta_exact2[2],mixed_before_betamat[2,i])
  mixed_before_rmse_beta2=c(mixed_before_rmse_beta2,RMSE)
  
  # bi MSE
  RMSE = RMSE_func(-1,mixed_before_bimat[1,i])
  mixed_before_rmse_b1_1=c(mixed_before_rmse_b1_1,RMSE)
  
  RMSE = RMSE_func(0, mixed_before_bimat[2,i])
  mixed_before_rmse_b2_1=c(mixed_before_rmse_b2_1,RMSE)
  
  RMSE = RMSE_func(1,mixed_before_bimat[3,i])
  mixed_before_rmse_b3_1=c(mixed_before_rmse_b3_1,RMSE)
  
  # f MSE
  fitted.func = eval.FEM.mixed(output_CPP$fit.FEM.mixed,
                               locations = loc)[,output_CPP$bestlambda]
  fitted.func1 = fitted.func[1:nlocs]
  fitted.func2 = fitted.func[(nlocs+1):(2*nlocs)]
  fitted.func3 = fitted.func[(2*nlocs+1):(3*nlocs)]
  
  RMSE=RMSE_func(func_evaluation1,fitted.func1)
  mixed_before_RMSE_f1=c(mixed_before_RMSE_f1, RMSE)
  
  RMSE=RMSE_func(func_evaluation2,fitted.func2)
  mixed_before_RMSE_f2=c(mixed_before_RMSE_f2, RMSE)
  
  RMSE=RMSE_func(func_evaluation3,fitted.func3)
  mixed_before_RMSE_f3=c(mixed_before_RMSE_f3, RMSE)
  
  #global MSE
  RMSE=RMSE_func(W%*%beta_exact1 + func_evaluation1,
                 W%*%mixed_before_betamat[,i] + cov1*mixed_before_bimat[1,i] + fitted.func1)
  mixed_before_global_RMSE1=c(mixed_before_global_RMSE1, RMSE)
  
  RMSE=RMSE_func(W%*%beta_exact2 + func_evaluation2,
                 W%*%mixed_before_betamat[,i] + cov1*mixed_before_bimat[2,i] + fitted.func2)
  mixed_before_global_RMSE2=c(mixed_before_global_RMSE2, RMSE)
  
  RMSE=RMSE_func(W%*%beta_exact3 + func_evaluation3,
                 W%*%mixed_before_betamat[,i] + cov1*mixed_before_bimat[3,i] + fitted.func3)
  mixed_before_global_RMSE3=c(mixed_before_global_RMSE3, RMSE)
  
  # total
  value = c(W%*%beta_exact1 + func_evaluation1, 
            W%*%beta_exact2 + func_evaluation2, 
            W%*%beta_exact3 + func_evaluation3)
  fitted.value = c(W%*%mixed_before_betamat[,i] + cov1*mixed_before_bimat[1,i] + fitted.func1,
                   W%*%mixed_before_betamat[,i] + cov1*mixed_before_bimat[2,i] + fitted.func2,
                   W%*%mixed_before_betamat[,i] + cov1*mixed_before_bimat[3,i] + fitted.func3)
  mixed_before_RMSE=c(mixed_before_RMSE, RMSE_func(value, fitted.value))
}

before_selected_lambda


### with NA #####
mixed_after_betamat=matrix(data=NA, nrow = 2, ncol=N)
mixed_after_rmse_beta1 =NULL
mixed_after_rmse_beta2 =NULL
mixed_after_bimat=matrix(data=NA, nrow = 3, ncol=N)
mixed_after_rmse_b1_1 =NULL
mixed_after_rmse_b2_1 =NULL
mixed_after_rmse_b3_1 =NULL
mixed_after_RMSE_f1 = NULL
mixed_after_RMSE_f2 = NULL
mixed_after_RMSE_f3 = NULL
mixed_after_global_RMSE1 = NULL
mixed_after_global_RMSE2 = NULL
mixed_after_global_RMSE3 = NULL
mixed_after_RMSE =NULL
after_selected_lambda=rep(0,N)

### num of NA ####
keep_orig = simul_observations_before
simul_observations_before = keep_orig

tot_len = length(mixed_obs)
ratio = 0.30
NA_len = tot_len*ratio
num_NA = NA_len/num_units # num of NA in each statistical unit
NA_len/tot_len
num_NA

### after NA lambda ####
### after NA, lambda changes.... needed to be fixed accordingly...
after_lambda= seq(2,8, length.out=30) #len: 30

length(which(is.na(simul_observations_before[i,,1])))
length(which(is.na(simul_observations_before[i,,2])))
length(which(is.na(simul_observations_before[i,,3])))

### after NA simulation #####
for(i in 1:N){
  for (j in 1:num_units) {
    set.seed(452185+i+j)
    NAIndx=sample(1:nlocs, num_NA)
    simul_observations_before[i,NAIndx,j] = NA
  }
  
  output_CPP2<-fdaPDEKIM:::smooth.FEM.mixed(locations = loc,
                                observations=simul_observations_before[i,,], 
                                covariates = covariates, #covariates added!
                                random_effect = c(1),
                                # FEMbasis=FEMbasis,
                                FEMbasis = mod_info$fit.FEM$FEMbasis, #reuse tree info
                                lambda=after_lambda,
                                GCV=GCVFLAG,
                                GCVmethod = GCVMETHODFLAG,
                                bary.locations = mod_info$bary.locations)
  
  after_selected_lambda[i]=output_CPP2$bestlambda
  
  # beta estimates
  mixed_after_betamat[1,i]=output_CPP2$beta[,output_CPP2$bestlambda][1]
  mixed_after_betamat[2,i]=output_CPP2$beta[,output_CPP2$bestlambda][2]
  
  # bi estimates
  mixed_after_bimat[1,i]=output_CPP2$b_i[,output_CPP2$bestlambda][1]
  mixed_after_bimat[2,i]=output_CPP2$b_i[,output_CPP2$bestlambda][2]
  mixed_after_bimat[3,i]=output_CPP2$b_i[,output_CPP2$bestlambda][3]
  
  # beta RMSE
  RMSE = RMSE_func(beta_exact2[1],mixed_after_betamat[1,i])
  mixed_after_rmse_beta1=c(mixed_after_rmse_beta1,RMSE)
  
  RMSE = RMSE_func(beta_exact2[2],mixed_after_betamat[2,i])
  mixed_after_rmse_beta2=c(mixed_after_rmse_beta2,RMSE)
  
  # bi MSE
  RMSE = RMSE_func(-1, mixed_after_bimat[1,i])
  mixed_after_rmse_b1_1=c(mixed_after_rmse_b1_1,RMSE)
  
  RMSE = RMSE_func(0, mixed_after_bimat[2,i])
  mixed_after_rmse_b2_1=c(mixed_after_rmse_b2_1,RMSE)
  
  RMSE = RMSE_func(1, mixed_after_bimat[3,i])
  mixed_after_rmse_b3_1=c(mixed_after_rmse_b3_1,RMSE)
  
  # f MSE
  fitted.func = eval.FEM.mixed(output_CPP2$fit.FEM.mixed,
                               locations = loc)[,output_CPP2$bestlambda]
  
  fitted.func1 = fitted.func[1:nlocs]
  fitted.func2 = fitted.func[(nlocs+1):(2*nlocs)]
  fitted.func3 = fitted.func[(2*nlocs+1):(3*nlocs)]
  
  RMSE=RMSE_func(func_evaluation1, fitted.func1)
  mixed_after_RMSE_f1=c(mixed_after_RMSE_f1, RMSE)
  
  RMSE=RMSE_func(func_evaluation2, fitted.func2)
  mixed_after_RMSE_f2=c(mixed_after_RMSE_f2, RMSE)
  
  RMSE=RMSE_func(func_evaluation3,fitted.func3)
  mixed_after_RMSE_f3=c(mixed_after_RMSE_f3, RMSE)
  
  
  #global MSE
  RMSE=RMSE_func(W%*%beta_exact1 + func_evaluation1,
                 W%*%mixed_after_betamat[,i] + cov1*mixed_after_bimat[1,i] + fitted.func1)
  mixed_after_global_RMSE1=c(mixed_after_global_RMSE1, RMSE)
  
  RMSE=RMSE_func(W%*%beta_exact2 + func_evaluation2,
                 W%*%mixed_after_betamat[,i] + cov1*mixed_after_bimat[2,i] + fitted.func2)
  mixed_after_global_RMSE2=c(mixed_after_global_RMSE2, RMSE)
  
  RMSE=RMSE_func(W%*%beta_exact3 + func_evaluation3,
                 W%*%mixed_after_betamat[,i] + cov1*mixed_after_bimat[3,i] + fitted.func3)
  mixed_after_global_RMSE3=c(mixed_after_global_RMSE3, RMSE)
  
  # total
  value = c(W%*%beta_exact1 + func_evaluation1, 
            W%*%beta_exact2 + func_evaluation2, 
            W%*%beta_exact3 + func_evaluation3)
  fitted.value = c(W%*%mixed_after_betamat[,i] + cov1*mixed_after_bimat[1,i] + fitted.func1,
                   W%*%mixed_after_betamat[,i] + cov1*mixed_after_bimat[2,i] + fitted.func2,
                   W%*%mixed_after_betamat[,i] + cov1*mixed_after_bimat[3,i] + fitted.func3)
  mixed_after_RMSE=c(mixed_after_RMSE, RMSE_func(value, fitted.value))
}

before_lambda[before_selected_lambda]
after_lambda[after_selected_lambda]

##### boxplot ####
prefix = paste0('NA_less_random_', num_NA)
prefix = paste0(prefix,'_')

png(paste0(prefix,'mixed_RMSE.png'))
boxplot(mixed_before_RMSE, mixed_after_RMSE, names = c('No NA','30% NA'), main='Low randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)


png(paste0(prefix,'beta1.png'))
boxplot(cbind(mixed_before_betamat[1,],mixed_after_betamat[1,]), names=c('No NA','30% NA'), main='Low randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
abline(h=beta_exact2[1], col='red')

png(paste0(prefix,'beta2.png'))
boxplot(cbind(mixed_before_betamat[2,],mixed_after_betamat[2,]), names=c('No NA','30% NA'), main='Low randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
abline(h=beta_exact2[2], col='red')


png(paste0(prefix,'b1_1.png'))
boxplot(cbind(mixed_before_bimat[1,],mixed_after_bimat[1,]), names=c('No NA','30% NA'), main='Low randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
abline(h=-1, col='red')

png(paste0(prefix,'b2_1.png'))
boxplot(cbind(mixed_before_bimat[2,],mixed_after_bimat[2,]), names=c('No NA','30% NA'), main='Low randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
abline(h=0, col='red')

png(paste0(prefix,'b3_1.png'))
boxplot(cbind(mixed_before_bimat[3,],mixed_after_bimat[3,]), names=c('No NA','30% NA'), main='Low randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
abline(h=1, col='red')


png(paste0(prefix,'f1 RMSE.png'))
boxplot(cbind(mixed_before_RMSE_f1, mixed_after_RMSE_f1), names = c('No NA','30% NA'), main='Low randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)

png(paste0(prefix,'f2 RMSE.png'))
boxplot(cbind(mixed_before_RMSE_f2, mixed_after_RMSE_f2), names = c('No NA','30% NA'), main='Low randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)

png(paste0(prefix,'f3 RMSE.png'))
boxplot(cbind(mixed_before_RMSE_f3, mixed_after_RMSE_f3), names = c('No NA','30% NA'), main='Low randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)

DF <- data.frame(mixed_before_RMSE_f1, mixed_after_RMSE_f1, 
                 mixed_before_RMSE_f2, mixed_after_RMSE_f2, 
                 mixed_before_RMSE_f3, mixed_after_RMSE_f3)
png(paste0(prefix,'f RMSE.png'))
boxplot(DF, col = c('red','blue'), at = c(1:2,4:5,7:8), xaxt = "n", main='Low randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
axis(side = 1, at = c(1.5,4.5,7.5), labels = c("f1","f2","f3"),
     cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
title(ylab="RMSE", line=2.3, cex.lab=2)
legend("topright", fill = c('red','blue'), legend = c('No NA','30% NA'), horiz = F,
       pt.cex=1.5, cex=1.5)

png(paste0(prefix,'beta1 RMSE.png'))
boxplot(cbind(mixed_before_rmse_beta1, mixed_after_rmse_beta1), names = c('No NA','30% NA'), main='Low randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)

png(paste0(prefix,'beta2 RMSE.png'))
boxplot(cbind(mixed_before_rmse_beta2, mixed_after_rmse_beta2), names = c('No NA','30% NA'), main='Low randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)

DF <- data.frame(mixed_before_rmse_beta1, mixed_after_rmse_beta1, 
                 mixed_before_rmse_beta2, mixed_after_rmse_beta2)
png(paste0(prefix,'beta RMSE.png'))
boxplot(DF, col = c('red','blue'), at = c(1:2,4:5), main='Low randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, xaxt = "n")
axis(side = 1, at = c(1.5,4.5), labels = c("beta1","beta2"),
     cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
title(ylab="RMSE", line=2.3, cex.lab=2)
legend("topright", fill = c('red','blue'), legend = c('No NA','30% NA'), horiz = F,
       pt.cex=1.5, cex=1.5)


png(paste0(prefix,'b1_1 RMSE.png'))
boxplot(cbind(mixed_before_rmse_b1_1, mixed_after_rmse_b1_1), names = c('No NA','30% NA'), main='Low randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)

png(paste0(prefix,'b2_1 RMSE.png'))
boxplot(cbind(mixed_before_rmse_b2_1, mixed_after_rmse_b2_1), names = c('No NA','30% NA'), main='Low randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)

png(paste0(prefix,'b3_1 RMSE.png'))
boxplot(cbind(mixed_before_rmse_b3_1, mixed_after_rmse_b3_1), names = c('No NA','30% NA'), main='Low randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)

DF <- data.frame(mixed_before_rmse_b1_1, mixed_after_rmse_b1_1, 
                 mixed_before_rmse_b2_1, mixed_after_rmse_b2_1, 
                 mixed_before_rmse_b3_1, mixed_after_rmse_b3_1)
png(paste0(prefix,'bi_1 RMSE.png'))
boxplot(DF, col = c('red','blue'), at = c(1:2,4:5,7:8), xaxt = "n", main='Low randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
axis(side = 1, at = c(1.5,4.5,7.5), labels = c("b1_1","b2_1","b3_1"),
     cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
title(ylab="RMSE", line=2.3, cex.lab=2)
legend("topleft", fill = c('red','blue'), legend = c('No NA','30% NA'), horiz = F,
       pt.cex=1.5, cex=1.5)

png(paste0(prefix,'unit1 global RMSE.png'))
boxplot(cbind(mixed_before_global_RMSE1, mixed_after_global_RMSE1), names = c('No NA','30% NA'), main='Low randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)

png(paste0(prefix,'unit2 global RMSE.png'))
boxplot(cbind(mixed_before_global_RMSE2, mixed_after_global_RMSE2), names = c('No NA','30% NA'), main='Low randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)

png(paste0(prefix,'unit3 global RMSE.png'))
boxplot(cbind(mixed_before_global_RMSE3, mixed_after_global_RMSE3), names = c('No NA','30% NA'), main='Low randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)


DF <- data.frame(mixed_before_global_RMSE1, mixed_after_global_RMSE1, 
                 mixed_before_global_RMSE2, mixed_after_global_RMSE2, 
                 mixed_before_global_RMSE3, mixed_after_global_RMSE3)
png(paste0(prefix,'unit global RMSE.png'))
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
boxplot(DF, col = c('red','blue'), at = c(1:2,4:5,7:8), xaxt = "n", main='Low randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
axis(side = 1, at = c(1.5,4.5,7.5), labels = c("unit1","unit2","unit3"),
     cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
legend("topright", inset=c(-0.3,0), fill = c('red','blue'), legend = c('No NA','30% NA'), horiz = F,
       pt.cex=1, cex=1)

#### before NA function image #####
library(fdaPDE)
library(fields)

mx = seq(min(mesh$nodes[,1]),max(mesh$nodes[,1]), length=500)
my = seq(min(mesh$nodes[,2]),max(mesh$nodes[,2]), length=500)
unif_grid = expand.grid(mx, my)
image_len = length(mx)

image_cov1=sin(2*pi*unif_grid[,1])*cos(2*pi*unif_grid[,2])
image_cov2=rnorm(dim(unif_grid)[1], mean=0, sd=2)
image_W=cbind(image_cov1, image_cov2)


eval_sol = eval.FEM.mixed(output_CPP$fit.FEM.mixed,
                          locations = unif_grid)[,output_CPP$bestlambda]

#### NEEED TO BE SAME FOR ALL! (CHECK AGAIN!)
min_before_func = floor(min(eval_sol,na.rm=TRUE))
max_before_func = ceiling(max(eval_sol,na.rm=TRUE))
# -3, 4

#### CHOSEN IN THE LAST MIN
func_min = -3
func_max = 4

num_units = output_CPP$fit.FEM.mixed$num_units
for (i in 1:num_units) {
  evalmat = matrix(eval_sol[((i-1)*image_len*image_len+1):(i*image_len*image_len)], 
                   nrow=image_len, 
                   ncol=image_len, byrow=F)
  
  title = 'before_NA_less_random_f'
  title =paste0(title,i)
  title = paste0(title,'.png')
  png(title)
  # x11()
  # image(mx, my, evalmat, axes = TRUE, xlab = 'x', ylab = 'y', col=hcl.colors(100, "YlOrRd", rev = TRUE))
  # contour(mx, my, evalmat, levels = seq(min, max, by=1), add = TRUE, labcex=1, col = "black")
  fields::image.plot(evalmat, zlim=c(func_min,func_max), col = hcl.colors(100, "YlOrRd", rev = TRUE))
  contour(evalmat, add=TRUE, levels = seq(func_min, func_max, by=1), labcex=1.5, col = "black")
}


#### after NA function image #####
eval_sol2 = eval.FEM.mixed(output_CPP2$fit.FEM.mixed,
                           locations = unif_grid)[,output_CPP2$bestlambda]

#### NEEED TO BE SAME FOR ALL! (CHECK AGAIN!)
min_after_func = floor(min(eval_sol2,na.rm=TRUE))
max_after_func = ceiling(max(eval_sol2,na.rm=TRUE))
# -2, 4 


num_units = output_CPP2$fit.FEM.mixed$num_units
for (i in 1:num_units) {
  evalmat = matrix(eval_sol2[((i-1)*image_len*image_len+1):(i*image_len*image_len)], 
                   nrow=image_len, 
                   ncol=image_len, byrow=F)
  
  title = 'after_NA_less_random_f'
  title =paste0(title,i)
  title = paste0(title,'.png')
  png(title)
  # x11()
  # image(mx, my, evalmat, axes = TRUE, xlab = 'x', ylab = 'y', col=hcl.colors(100, "YlOrRd", rev = TRUE))
  # contour(mx, my, evalmat, levels = seq(min, max, by=1), add = TRUE, labcex=1, col = "black")
  fields::image.plot(evalmat, zlim=c(func_min,func_max), col = hcl.colors(100, "YlOrRd", rev = TRUE))
  contour(evalmat, add=TRUE, levels = seq(func_min, func_max, by=1), labcex=1.5, col = "black")
}

####.####
##### NO random-effect ####
rm(list=ls())
setwd("C:/Users/JIYOUNG KIM/Desktop/tests - mixed")
library(fdaPDE)
library(mgcv)

data(horseshoe2D)
mesh=create.mesh.2D(nodes=boundary_nodes, segments = boundary_segments)
mesh=refine.mesh.2D(mesh, maximum_area = 0.025, minimum_angle = 30)
FEMbasis <- create.FEM.basis(mesh)

GCVFLAG=T
GCVMETHODFLAG='Exact' #For now, Exact because 'Stochastic' changes every time..
nnodes = dim(mesh$nodes)[1]

### lambda####
before_lambda= 10^seq(-2,1,by=0.1) #len: 31

### boundary ####
bnd = boundary_nodes
bound <- list(list(x = bnd[,1], y = bnd[,2])) #, f = rep(0, nrow(crds))))
N <- 10
gx <- seq(min(bnd[,1]), max(bnd[,1]), len = N)
gy <- seq(min(bnd[,2]), max(bnd[,2]), len = N)
gp <- expand.grid(gx, gy)
names(gp) <- c("x","y")
knots <- data.frame(v=rep(seq(-.5,3,by=.5),4),
                    w=rep(c(-.6,-.3,.3,.6),rep(8,4)))
names(knots) <- c("x", "y")
names(bound[[1]]) <- c("x", "y") #, "f")
bound[[1]]$y[c(10:28)] = bound[[1]]$y[c(10:28)] +0.005
bound[[1]]$y[c(82:100)] = bound[[1]]$y[c(82:100)] -0.005
bound[[1]]$x[c(1:9,101:108)] = bound[[1]]$x[c(1:9,101:108)]-0.005
bound[[1]]$x[c(29:35,75:81)] = bound[[1]]$x[c(29:35,75:81)]+0.005
bound[[1]]$y[c(36:54)] = bound[[1]]$y[c(36:54)]-0.005
bound[[1]]$x[c(55)]  = bound[[1]]$x[c(55)]+0.005
bound[[1]]$y[c(56:74)] = bound[[1]]$y[c(56:74)]+0.005


# Exact data - Nodes locations (for C shape)
fs.test.time<-function(x,y,t) {
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


## locations ####
output_CPP1 = smooth.FEM(observations = sin(pi*mesh$nodes[,1]), #observations at mesh nodes
                         FEMbasis = FEMbasis,
                         lambda = before_lambda,
                         GCV=GCVFLAG,
                         GCVmethod = GCVMETHODFLAG)

#### uniform Grid Initialization
gridnum = 20
minV1 = min(mesh$nodes[,1])
maxV1 = max(mesh$nodes[,1])
minV2 = min(mesh$nodes[,2])
maxV2 = max(mesh$nodes[,2])

x = seq(minV1+0.1, maxV1-0.1, length.out = gridnum) ### SHOULD ADJUST THIS! (OR ELSE, ERROR IN SOAP FILM)
y= seq(minV2+0.1, maxV2-0.1, length.out = gridnum)
unif_grid <- expand.grid(x = x, y = y)

points1=eval.FEM(output_CPP1$fit.FEM,
                 locations = unif_grid)

loc = unif_grid[-which(is.na(points1)),]

# randomize
set.seed(5847947)
ind = sample.int(dim(loc)[1])[1:100] ##100 locations points as default
loc <- loc[ind,]
nlocs = dim(loc)[1]


### covariates ####
cov1=sin(2*pi*loc[,1])*cos(2*pi*loc[,2])
cov2=rnorm(nlocs, mean=0, sd=2)
W=cbind(cov1,cov2)

mesh_cov1=sin(2*pi*mesh$nodes[,1])*cos(2*pi*mesh$nodes[,2])
mesh_cov2=rnorm(nnodes, mean=0, sd=2)
mesh_W=cbind(mesh_cov1, mesh_cov2)

# Fix betas
beta_exact1=c(3, 0.5)
beta_exact2=c(3, 0.5)
beta_exact3=c(3, 0.5)

mean(5*cov1 + 5*cov2) 

num_units = 3
RMSE_func<-function(f,g) {
  sqrt(mean((f-g)^2))
}

func_evaluation1 = numeric(nlocs)
func_evaluation1 =  fs.test.time(x=loc[,1], y=loc[,2],  t=rep(0.5, nlocs))

func_evaluation2 = numeric(nlocs)
func_evaluation2 = fs.test.time(x=loc[,1], y=loc[,2],  t=rep(1, nlocs))

func_evaluation3 = numeric(nlocs)
func_evaluation3 = fs.test.time(x=loc[,1], y=loc[,2],  t=rep(1.5, nlocs))


func1 = numeric(nnodes)
func1 = fs.test.time(x=mesh$nodes[,1], y=mesh$nodes[,2],  t=rep(0.5, nnodes))

func2 = numeric(nnodes)
func2 = fs.test.time(x=mesh$nodes[,1], y=mesh$nodes[,2],  t=rep(1, nnodes))

func3 = numeric(nnodes)
func3 = fs.test.time(x=mesh$nodes[,1], y=mesh$nodes[,2],  t=rep(1.5, nnodes))


# try to set as similar range
range(func_evaluation1)[2] - range(func_evaluation1)[1]
range(W%*%beta_exact1)[2] - range(W%*%beta_exact1)[1]

range(func_evaluation2)[2] - range(func_evaluation2)[1]
range(W%*%beta_exact2)[2] - range(W%*%beta_exact2)[1]

range(func_evaluation3)[2] - range(func_evaluation3)[1]
range(W%*%beta_exact3)[2] - range(W%*%beta_exact3)[1]

DF <- data.frame(W%*%beta_exact1, func_evaluation1,
                 W%*%beta_exact2, func_evaluation2,
                 W%*%beta_exact3, func_evaluation3)
x11()
# png('range_cov_func.png')
# par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
boxplot(DF, col = c('red','blue'), at = c(1:2,4:5,7:8), xaxt = "n", main='range',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
axis(side = 1, at = c(1.5,4.5,7.5), labels = c("unit1","unit2","unit3"),
     cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)

mixed_loc = rbind(loc, loc, loc)
mixed_cov1 = c(cov1, cov1, cov1)
mixed_cov2 = c(cov2, cov2, cov2)
# mixed_obs = c(W%*%beta_exact1+func_evaluation1, 
#               W%*%beta_exact2+func_evaluation2,
#               W%*%beta_exact3+func_evaluation3)
mixed_obs = cbind(W%*%beta_exact1+func_evaluation1,
                  W%*%beta_exact2+func_evaluation2,
                  W%*%beta_exact3+func_evaluation3)
covariates=rbind(W,W,W)
mixed_nlocs= nlocs*3

#### range ####
ran=range(c(func_evaluation1, func_evaluation2, func_evaluation3))

ran[2]-ran[1]
0.05*(ran[2]-ran[1])

# Set number of simulation trials
N=50

simul_observations_before <- array(rep(0, N*nrow(loc)*num_units), 
                                   dim=c(N, nrow(loc), num_units))

RMSE_func<-function(f,g) {
  sqrt(mean((f-g)^2))
}


mod_info = fdaPDEKIM:::smooth.FEM.mixed(locations = loc,
                            observations = mixed_obs, #obs doesn't effect tree, bary, dof
                            covariates = covariates,
                            FEMbasis = FEMbasis,
                            lambda = before_lambda,
                            GCV = GCVFLAG)


mixed_before_betamat=matrix(data=NA, nrow = 2, ncol=N)
mixed_before_rmse_beta1 =NULL
mixed_before_rmse_beta2 =NULL
mixed_before_RMSE_f1 = NULL
mixed_before_RMSE_f2 = NULL
mixed_before_RMSE_f3 = NULL
mixed_before_global_RMSE1 = NULL
mixed_before_global_RMSE2 = NULL
mixed_before_global_RMSE3 = NULL
mixed_before_RMSE =NULL
before_selected_lambda=rep(0,N)

#### before NA simulation ######
for(i in 1:N){
  set.seed(54425 + i)
  data = mixed_obs + rnorm(nlocs*3,mean=0,sd=0.05*(ran[2]-ran[1])) #add column-wise (convertable with matrix and vector)
  simul_observations_before[i,,]=data
  
  output_CPP<-fdaPDEKIM:::smooth.FEM.mixed(locations = loc,
                               observations=data, 
                               covariates = covariates, #covariates added!
                               # FEMbasis=FEMbasis, 
                               FEMbasis = mod_info$fit.FEM.mixed$FEMbasis, #reuse tree info
                               lambda=before_lambda,
                               GCV=GCVFLAG,
                               GCVmethod = GCVMETHODFLAG,
                               DOF_matrix = mod_info$edf) #reuse dof info
  
  before_selected_lambda[i]=output_CPP$bestlambda
  
  # beta estimates
  mixed_before_betamat[1,i]=output_CPP$beta[,output_CPP$bestlambda][1]
  mixed_before_betamat[2,i]=output_CPP$beta[,output_CPP$bestlambda][2]
  
  # beta RMSE
  RMSE = RMSE_func(beta_exact2[1],mixed_before_betamat[1,i])
  mixed_before_rmse_beta1=c(mixed_before_rmse_beta1,RMSE)
  
  RMSE = RMSE_func(beta_exact2[2],mixed_before_betamat[2,i])
  mixed_before_rmse_beta2=c(mixed_before_rmse_beta2,RMSE)
  
  # f MSE
  fitted.func = eval.FEM.mixed(output_CPP$fit.FEM.mixed,
                               locations = loc)[,output_CPP$bestlambda]
  fitted.func1 = fitted.func[1:nlocs]
  fitted.func2 = fitted.func[(nlocs+1):(2*nlocs)]
  fitted.func3 = fitted.func[(2*nlocs+1):(3*nlocs)]
  
  RMSE=RMSE_func(func_evaluation1,fitted.func1)
  mixed_before_RMSE_f1=c(mixed_before_RMSE_f1, RMSE)
  
  RMSE=RMSE_func(func_evaluation2,fitted.func2)
  mixed_before_RMSE_f2=c(mixed_before_RMSE_f2, RMSE)
  
  RMSE=RMSE_func(func_evaluation3,fitted.func3)
  mixed_before_RMSE_f3=c(mixed_before_RMSE_f3, RMSE)
  
  #global MSE
  RMSE=RMSE_func(W%*%beta_exact1 + func_evaluation1,
                 W%*%mixed_before_betamat[,i] + fitted.func1)
  mixed_before_global_RMSE1=c(mixed_before_global_RMSE1, RMSE)
  
  RMSE=RMSE_func(W%*%beta_exact2 + func_evaluation2,
                 W%*%mixed_before_betamat[,i]  + fitted.func2)
  mixed_before_global_RMSE2=c(mixed_before_global_RMSE2, RMSE)
  
  RMSE=RMSE_func(W%*%beta_exact3 + func_evaluation3,
                 W%*%mixed_before_betamat[,i]  + fitted.func3)
  mixed_before_global_RMSE3=c(mixed_before_global_RMSE3, RMSE)
  
  # total
  value = c(W%*%beta_exact1 + func_evaluation1, 
            W%*%beta_exact2 + func_evaluation2, 
            W%*%beta_exact3 + func_evaluation3)
  fitted.value = c(W%*%mixed_before_betamat[,i] + fitted.func1,
                   W%*%mixed_before_betamat[,i] + fitted.func2,
                   W%*%mixed_before_betamat[,i] + fitted.func3)
  mixed_before_RMSE=c(mixed_before_RMSE, RMSE_func(value, fitted.value))
}

lambda
before_selected_lambda

### with NA #####
mixed_after_betamat=matrix(data=NA, nrow = 2, ncol=N)
mixed_after_rmse_beta1 =NULL
mixed_after_rmse_beta2 =NULL
mixed_after_RMSE_f1 = NULL
mixed_after_RMSE_f2 = NULL
mixed_after_RMSE_f3 = NULL
mixed_after_global_RMSE1 = NULL
mixed_after_global_RMSE2 = NULL
mixed_after_global_RMSE3 = NULL
mixed_after_RMSE =NULL
after_selected_lambda=rep(0,N)

### num of NA ####
keep_orig = simul_observations_before
simul_observations_before = keep_orig

tot_len = length(mixed_obs)
ratio = 0.30
NA_len = tot_len*ratio
num_NA = NA_len/num_units # num of NA in each statistical unit
NA_len/tot_len
num_NA

### after NA lambda ####
### after NA, lambda changes.... needed to be fixed accordingly...
after_lambda= 10^seq(0,1, length.out=30) #len: 30


### after NA simulation #####
for(i in 1:N){
  for (j in 1:num_units) {
    set.seed(452185+i+j)
    NAIndx=sample(1:nlocs, num_NA)
    simul_observations_before[i,NAIndx,j] = NA
  }
  
  output_CPP2<-fdaPDEKIM:::smooth.FEM.mixed(locations = loc,
                                observations=simul_observations_before[i,,], 
                                covariates = covariates, #covariates added!
                                # FEMbasis=FEMbasis, 
                                FEMbasis = mod_info$fit.FEM$FEMbasis, #reuse tree info
                                lambda=after_lambda,
                                GCV=GCVFLAG,
                                GCVmethod = GCVMETHODFLAG,
                                bary.locations = mod_info$bary.locations)
  
  after_selected_lambda[i]=output_CPP2$bestlambda
  
  # beta estimates
  mixed_after_betamat[1,i]=output_CPP2$beta[,output_CPP2$bestlambda][1]
  mixed_after_betamat[2,i]=output_CPP2$beta[,output_CPP2$bestlambda][2]
  
  # beta RMSE
  RMSE = RMSE_func(beta_exact2[1],mixed_after_betamat[1,i])
  mixed_after_rmse_beta1=c(mixed_after_rmse_beta1,RMSE)
  
  RMSE = RMSE_func(beta_exact2[2],mixed_after_betamat[2,i])
  mixed_after_rmse_beta2=c(mixed_after_rmse_beta2,RMSE)
  
  
  # f MSE
  fitted.func = eval.FEM.mixed(output_CPP2$fit.FEM.mixed,
                               locations = loc)[,output_CPP2$bestlambda]
  
  fitted.func1 = fitted.func[1:nlocs]
  fitted.func2 = fitted.func[(nlocs+1):(2*nlocs)]
  fitted.func3 = fitted.func[(2*nlocs+1):(3*nlocs)]
  
  RMSE=RMSE_func(func_evaluation1, fitted.func1)
  mixed_after_RMSE_f1=c(mixed_after_RMSE_f1, RMSE)
  
  RMSE=RMSE_func(func_evaluation2, fitted.func2)
  mixed_after_RMSE_f2=c(mixed_after_RMSE_f2, RMSE)
  
  RMSE=RMSE_func(func_evaluation3,fitted.func3)
  mixed_after_RMSE_f3=c(mixed_after_RMSE_f3, RMSE)
  
  
  #global MSE
  RMSE=RMSE_func(W%*%beta_exact1 + func_evaluation1,
                 W%*%mixed_after_betamat[,i]  + fitted.func1)
  mixed_after_global_RMSE1=c(mixed_after_global_RMSE1, RMSE)
  
  RMSE=RMSE_func(W%*%beta_exact2 + func_evaluation2,
                 W%*%mixed_after_betamat[,i] + fitted.func2)
  mixed_after_global_RMSE2=c(mixed_after_global_RMSE2, RMSE)
  
  RMSE=RMSE_func(W%*%beta_exact3 + func_evaluation3,
                 W%*%mixed_after_betamat[,i]  + fitted.func3)
  mixed_after_global_RMSE3=c(mixed_after_global_RMSE3, RMSE)
  
  # total
  value = c(W%*%beta_exact1 + func_evaluation1, 
            W%*%beta_exact2 + func_evaluation2, 
            W%*%beta_exact3 + func_evaluation3)
  fitted.value = c(W%*%mixed_after_betamat[,i] + fitted.func1,
                   W%*%mixed_after_betamat[,i] + fitted.func2,
                   W%*%mixed_after_betamat[,i]  + fitted.func3)
  mixed_after_RMSE=c(mixed_after_RMSE, RMSE_func(value, fitted.value))
}

before_selected_lambda
after_selected_lambda

##### boxplot ####
prefix = paste0('NA_no_random_', num_NA)
prefix = paste0(prefix,'_')

png(paste0(prefix,'mixed_RMSE.png'))
boxplot(mixed_before_RMSE, mixed_after_RMSE, names = c('No NA','30% NA'), main='No randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)


png(paste0(prefix,'beta1.png'))
boxplot(cbind(mixed_before_betamat[1,],mixed_after_betamat[1,]), names=c('No NA','30% NA'), main='No randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
abline(h=beta_exact2[1], col='red')

png(paste0(prefix,'beta2.png'))
boxplot(cbind(mixed_before_betamat[2,],mixed_after_betamat[2,]), names=c('No NA','30% NA'), main='No randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
abline(h=beta_exact2[2], col='red')


png(paste0(prefix,'f1 RMSE.png'))
boxplot(cbind(mixed_before_RMSE_f1, mixed_after_RMSE_f1), names = c('No NA','30% NA'), main='No randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)

png(paste0(prefix,'f2 RMSE.png'))
boxplot(cbind(mixed_before_RMSE_f2, mixed_after_RMSE_f2), names = c('No NA','30% NA'), main='No randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)

png(paste0(prefix,'f3 RMSE.png'))
boxplot(cbind(mixed_before_RMSE_f3, mixed_after_RMSE_f3), names = c('No NA','30% NA'), main='No randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)

DF <- data.frame(mixed_before_RMSE_f1, mixed_after_RMSE_f1, 
                 mixed_before_RMSE_f2, mixed_after_RMSE_f2, 
                 mixed_before_RMSE_f3, mixed_after_RMSE_f3)
png(paste0(prefix,'f RMSE.png'))
boxplot(DF, col = c('red','blue'), at = c(1:2,4:5,7:8), xaxt = "n", main='No randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
axis(side = 1, at = c(1.5,4.5,7.5), labels = c("f1","f2","f3"),
     cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
title(ylab="RMSE", line=2.3, cex.lab=2)
legend("topright", fill = c('red','blue'), legend = c('No NA','30% NA'), horiz = F,
       pt.cex=1.5, cex=1.5)

png(paste0(prefix,'beta1 RMSE.png'))
boxplot(cbind(mixed_before_rmse_beta1, mixed_after_rmse_beta1), names = c('No NA','30% NA'), main='No randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)

png(paste0(prefix,'beta2 RMSE.png'))
boxplot(cbind(mixed_before_rmse_beta2, mixed_after_rmse_beta2), names = c('No NA','30% NA'), main='No randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)

DF <- data.frame(mixed_before_rmse_beta1, mixed_after_rmse_beta1, 
                 mixed_before_rmse_beta2, mixed_after_rmse_beta2)
png(paste0(prefix,'beta RMSE.png'))
boxplot(DF, col = c('red','blue'), at = c(1:2,4:5), main='No randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, xaxt = "n")
axis(side = 1, at = c(1.5,4.5), labels = c("beta1","beta2"),
     cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
title(ylab="RMSE", line=2.3, cex.lab=2)
legend("topleft", fill = c('red','blue'), legend = c('No NA','30% NA'), horiz = F,
       pt.cex=1.5, cex=1.5)


png(paste0(prefix,'unit1 global RMSE.png'))
boxplot(cbind(mixed_before_global_RMSE1, mixed_after_global_RMSE1), names = c('No NA','30% NA'), main='No randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)

png(paste0(prefix,'unit2 global RMSE.png'))
boxplot(cbind(mixed_before_global_RMSE2, mixed_after_global_RMSE2), names = c('No NA','30% NA'), main='No randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)

png(paste0(prefix,'unit3 global RMSE.png'))
boxplot(cbind(mixed_before_global_RMSE3, mixed_after_global_RMSE3), names = c('No NA','30% NA'), main='No randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)


DF <- data.frame(mixed_before_global_RMSE1, mixed_after_global_RMSE1, 
                 mixed_before_global_RMSE2, mixed_after_global_RMSE2, 
                 mixed_before_global_RMSE3, mixed_after_global_RMSE3)
png(paste0(prefix,'unit global RMSE.png'))
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
boxplot(DF, col = c('red','blue'), at = c(1:2,4:5,7:8), xaxt = "n", main='No randomness',
        cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
axis(side = 1, at = c(1.5,4.5,7.5), labels = c("unit1","unit2","unit3"),
     cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
legend("topright", inset=c(-0.3,0), fill = c('red','blue'), legend = c('No NA','30% NA'), horiz = F,
       pt.cex=1, cex=1)

#### image ####
#### before NA function image #####
library(fdaPDE)
library(fields)

mx = seq(min(mesh$nodes[,1]),max(mesh$nodes[,1]), length=500)
my = seq(min(mesh$nodes[,2]),max(mesh$nodes[,2]), length=500)
unif_grid = expand.grid(mx, my)
image_len = length(mx)

image_cov1=sin(2*pi*unif_grid[,1])*cos(2*pi*unif_grid[,2])
image_cov2=rnorm(dim(unif_grid)[1], mean=0, sd=2)
image_W=cbind(image_cov1, image_cov2)


eval_sol = eval.FEM.mixed(output_CPP$fit.FEM.mixed,
                          locations = unif_grid)[,output_CPP$bestlambda]

#### NEEED TO BE SAME FOR ALL! (CHECK AGAIN!)
min_before_func = floor(min(eval_sol,na.rm=TRUE))
max_before_func = ceiling(max(eval_sol,na.rm=TRUE))
# -3, 4

#### CHOSEN IN THE LAST MIN
func_min = -3
func_max = 4


num_units = output_CPP$fit.FEM.mixed$num_units
for (i in 1:num_units) {
  evalmat = matrix(eval_sol[((i-1)*image_len*image_len+1):(i*image_len*image_len)], 
                   nrow=image_len, 
                   ncol=image_len, byrow=F)
  
  title = 'before_NA_no_random_f'
  title =paste0(title,i)
  title = paste0(title,'.png')
  png(title)
  # x11()
  # image(mx, my, evalmat, axes = TRUE, xlab = 'x', ylab = 'y', col=hcl.colors(100, "YlOrRd", rev = TRUE))
  # contour(mx, my, evalmat, levels = seq(min, max, by=1), add = TRUE, labcex=1, col = "black")
  fields::image.plot(evalmat, zlim=c(func_min,func_max), col = hcl.colors(100, "YlOrRd", rev = TRUE))
  contour(evalmat, add=TRUE, levels = seq(func_min, func_max, by=1), labcex=1.5, col = "black")
}


#### after NA function image #####
eval_sol2 = eval.FEM.mixed(output_CPP2$fit.FEM.mixed,
                          locations = unif_grid)[,output_CPP2$bestlambda]

#### NEEED TO BE SAME FOR ALL! (CHECK AGAIN!)
min_after_func = floor(min(eval_sol2,na.rm=TRUE))
max_after_func = ceiling(max(eval_sol2,na.rm=TRUE))
#-2, 4

num_units = output_CPP2$fit.FEM.mixed$num_units
for (i in 1:num_units) {
  evalmat = matrix(eval_sol2[((i-1)*image_len*image_len+1):(i*image_len*image_len)], 
                   nrow=image_len, 
                   ncol=image_len, byrow=F)
  
  title = 'after_NA_no_random_f'
  title =paste0(title,i)
  title = paste0(title,'.png')
  png(title)
  # x11()
  # image(mx, my, evalmat, axes = TRUE, xlab = 'x', ylab = 'y', col=hcl.colors(100, "YlOrRd", rev = TRUE))
  # contour(mx, my, evalmat, levels = seq(min, max, by=1), add = TRUE, labcex=1, col = "black")
  fields::image.plot(evalmat, zlim=c(func_min,func_max), col = hcl.colors(100, "YlOrRd", rev = TRUE))
  contour(evalmat, add=TRUE, levels = seq(func_min, func_max, by=1), labcex=1.5, col = "black")
}


