seizures = read.table('data.txt', skip=36
            ,col.names = c('id', 'treatment', 'age', 'baseline',
            'count1', 'count2', 'count3', 'count4'))
seizures[, 'id'] = as.factor(seizures[, 'id'])

# check data types
sapply(seizures, class)

# check for overdispersion for fixed t
sapply(seizures[, 4:8], var) / sapply(seizures[, 4:8], mean)

# check for overdispersion for fixed i
apply(seizures[,4:8], 1, var) / apply(seizures[,4:8], 1, mean)

pivoter = function(i, dframe=seizures){
  nx = 4  # number of fixed columns
  out = dframe[, 1:nx]
  out['period'] = i - nx
  out['count'] = dframe[, i]
  out
}

longseizure = do.call(rbind, lapply(5:8, pivoter))

data.transformed = data.frame(
  y = longseizure$count
  , Base = log((1/4) * longseizure$baseline)
  , Trt = longseizure$treatment
  , Age = log(longseizure$age)
  , V4 = longseizure$period == 4
  , vj = longseizure$period
  , ID = longseizure$id
)

data.transformed$vj[data.transformed$vj == 1] = -0.3
data.transformed$vj[data.transformed$vj == 2] = -0.1
data.transformed$vj[data.transformed$vj == 3] = 0.1
data.transformed$vj[data.transformed$vj == 4] = 0.3

data.transformed = data.transformed[order(data.transformed$ID,data.transformed$vj),]

library(nlme)
library(MASS)
fit.PQL = glmmPQL(fixed = y ~ Base * Trt + Age + V4 + vj,
                  random = ~ vj | ID,
                  family = poisson,
                  data = data.transformed,
                  verbose = FALSE)

b = fixef(fit.PQL)
tau = fit.PQL$sigma
s1 = sqrt(getVarCov(fit.PQL)[1,1])
s2 = sqrt(getVarCov(fit.PQL)[2,2])
p = getVarCov(fit.PQL)[1,2]/(s1*s2)

phi_0 = c(b,s1,s2,p,tau)

# 59*4 x 7 matrix
X = model.matrix(fit.PQL, data = data.transformed)
y = data.transformed$y

### testing ###
uij = function(phi, X){
  b = phi[1:7]
  s1 = phi[8]
  s2 = phi[9]
  p = phi[10]
  tau = phi[11]
  
  # 59*4 x 1 vector
  exp(1)^{X%*%b + 0.5*(s1^2 + X[,6]^2*s2^2 + 2*X[,6]*p*s1*s2 + 
      tau^2)}
}

Vi = function(phi, u_ij, X){
  b = phi[1:7]
  s1 = phi[8]
  s2 = phi[9]
  p = phi[10]
  tau = phi[11]
  
  Vdiag = list()
  
  for(i in 1:(dim(X)[1]/4)){
    seq = 1:4 + (i-1)*4
    # covariance
    V = (u_ij[seq] %*% t(u_ij[seq])) * 
        (exp(1)^{s1^2 + v_add*p*s1*s2 + v_times*s2^2} - 1)

    # variance
    diag(V) = u_ij[seq] + u_ij[seq]^2*(exp(1)^{s1^2 + 
        2*p*s1*s2*X[seq,6] + X[seq,6]^2*s2^2 + tau^2} - 1)

    Vdiag[[i]] = V
  }
  bdiag(Vdiag)
}

max(log(uij(rnorm(11,mean=0,sd=1),y)))
################

### needed to calculate Vi's in G function ###
# will be a 59*4 x 59*4 matrix
v_add = outer(c(-0.3,-0.1,0.1,0.3),c(-0.3,-0.1,0.1,0.3),"+")
v_times = outer(c(-0.3,-0.1,0.1,0.3),c(-0.3,-0.1,0.1,0.3),"*")
library(Matrix)

### calculating the G function ###
G = function(phi, y, X){
  b = phi[1:7]
  s1 = phi[8]
  s2 = phi[9]
  p = phi[10]
  tau = phi[11]
  
  # 59*4 x 1 vector
  u_ij = exp(1)^{X%*%b + 0.5*(s1^2 + X[,6]^2*s2^2 +
         2*X[,6]*p*s1*s2 + tau^2)}
  
  ### u dot matrix ###
  # columns 1 - 7
  u1 = u_ij * X[,1]
  u2 = u_ij * X[,2]
  u3 = u_ij * X[,3]
  u4 = u_ij * X[,4]
  u5 = u_ij * X[,5]
  u6 = u_ij * X[,6]
  u7 = u_ij * X[,7]

  # column 8
  u8 = u_ij * (s1 + X[,6]*p*s2)

  # column 9
  u9 = u_ij * (X[,6]^2*s2 + X[,6]*p*s1)

  # column 10
  u10 = u_ij * X[,6]*s1*s2

  # column 11
  u11 = u_ij * tau

  # 59*4 x 11 matrix
  udot = cbind(u1,u2,u3,u4,u5,u6,u7,u8,u9,u10,u11)
  # udot = do.call(cbind, lapply(paste0("u",1:11),as.name))
  
  # setup empty list to hold blocks
  Vdiag = list()
  
  for(i in 1:(dim(X)[1]/4)){
    seq = 1:4 + (i-1)*4
    # covariance
    V = (u_ij[seq] %*% t(u_ij[seq])) * 
      (exp(1)^{s1^2 + v_add*p*s1*s2 + v_times*s2^2} - 1)
    
    # variance
    diag(V) = u_ij[seq] + u_ij[seq]^2*(exp(1)^{s1^2 + 
        2*p*s1*s2*X[seq,6] + X[seq,6]^2*s2^2 + tau^2} - 1)
    
    Vdiag[[i]] = V
  }
  Vi = bdiag(Vdiag)

  G = t(udot)%*%solve(Vi)%*%(y - u_ij)
  sum(G^2)
}

### testing ###
phi1 = optim(rnorm(11), G, gr = NULL, y, X,
       method = "L-BFGS-B", lower = c(rep(-Inf,7),0,0,-1,0),
       upper = c(rep(Inf,9),1,Inf),
       control = list(factr=10^4))$par
G(phi1,y,X)
################

iterateG = function(phi, y, X){
  phi_new = optim(phi, G, gr = NULL, y, X,
              method = "L-BFGS-B", lower = c(rep(-Inf,7),0,0,-1,0),
              upper = c(rep(Inf,9),1,Inf))$par
  return(phi_new)
}



### bootstrap ###
# setup
set.seed(231)
library(mvtnorm)
B = 100
result = matrix(0, 11, B)

# initial vector
phi_0 = c(b,s1,s2,p,tau)
phi = phi_0

for(i in 1:B) {
  V_a = matrix(c(phi[8]^2,phi[10]*phi[8]*phi[9],
        phi[10]*phi[8]*phi[9],phi[9]^2),2,2)
  a = rmvnorm(59,c(0,0), V_a)
  e = rnorm(59*4, mean = 0, sd = phi[11])
  
  a = a[rep(1:59,rep(4,59)),]

  u_b = exp(1)^{X%*%phi[1:7] + a[,1] + a[,2]*X[,6] + e}
  y_b = rpois(59*4,u_b)
  
  cat(sprintf("Iteration %i\n", i))
  
  # update step
  phi = iterateG(phi,y_b,X)

  # each column is one bootstrap instance of the vector phi
  # each row is the collection of bootstrap instances of a single 
  #   component of phi
  result[,i] = phi
}

# get bootstrap standard errors for each variable
sqrt(apply(result,1,var))

# seed 232
# bootstrap standard error results
# [1] 0.014007445 0.036007544 0.012816919 0.043830689 0.004026900 0.003789310
# [7] 0.035229762 0.049753974 0.001470748 0.009205787 0.129399240



### jack-knife ###
# setup
results = matrix(0,11,59)

for(i in 1:59){
  seq = 1:4 + 4*(i-1)
  cat(sprintf("Iteration %i\n", i))
  results[,i] = iterateG(phi_0,y[-seq],X[-seq,])
}

# get phi estimate from complete data
phi_full = iterateG(phi_0,y,X)

#
# (Intercept)         Base          Trt          Age       V4TRUE 
# -1.527195847  0.555442898 -1.011529533 -1.232648564 -0.225133293 
# vj     Base:Trt                                        
# -0.114341214  0.394006390  1.450268837  0.520737319  0.008849646 
# 3.373077950

phi_full_mat = matrix(rep(phi_full,59),11,59)

# calculate jack-knife standard errors
jackvar = (58/59) * (results - phi_full_mat)%*%t(results - phi_full_mat)
sqrt(diag(jackvar))

# jack-knife standard error results
# [1] 0.92788336798 3.17239505363 1.23139922765 3.28275041756 0.36997545557
# [6] 0.19146177353 3.18197048754 1.51689456139 0.02836353322 0.12804179799
# [11] 2.56711686156

# setting up plotting area
par(oma = c(0,0,0,0), mar=c(3.1, 3.1, 3.1, 2.1), mgp=c(3,0.7,0))

# jack-knife boxplot
boxplot(data.frame(t(results)),names = c("Int","Base","Trt",
      "Age","V4","vj","Base:Trt","s1","s2","p","tau"),
      main = "Jack-Knife Estimtes",cex.axis=0.5,cex=0.5)

# bootstrap boxplot
boxplot(data.frame(t(result)),names = c("Int","Base","Trt",
      "Age","V4","vj","Base:Trt","s1","s2","p","tau"),
      main = "Bootstrap Estimates",cex.axis=0.5,cex=0.5)
