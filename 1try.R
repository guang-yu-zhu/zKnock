#  ---
set.seed(2022)
p=30; n=100; k=15
mu = rep(0,p); Sigma = diag(p)
X = matrix(rnorm(n*p),n)
nonzero = 1:k
beta = 3.5 * (1:p %in% nonzero)
y = X %*% beta + rnorm(n)
Ac = 1:k
Ic = setdiff(1:p,Ac)
metric_fun = function(selected){
  power = sum(selected%in%Ac)/max(1, length(Ac))
  fdp = sum(Ic%in%selected)/max(1, length(selected))
  c(power,fdp)
}

# Knockoff Procedure
Xk = create.knockoff(X = X, type = 'shrink', num = 1)
res1 = knockoff.filter(X,y,Xk,statistic = stat.glmnet_coefdiff,family='gaussian')
res2 = knockoff.filter(X,y,Xk,statistic = stat.random_forest,family='gaussian')
res3 = knockoff.filter(X,y,Xk,statistic = stat.xgboost)
metric_fun(res1$s)
metric_fun(res2$s)
metric_fun(res3$s)


#--------
load('KOBT-results.Rdata')
# Z <- create.knockoff(X = as.matrix(X), type = "pc", num = 1,num.comp = 4)
result <- vector(mode = "list", length = length(Z))
# for (i in 1:length(Z)) {
#   cat('--Calculate',i,'knockoff statistics.\n')
#   x <- cbind(X, Z[[i]])
#   dtrain <- xgboost::xgb.DMatrix(as.matrix(x), label = y)
#   fit.model <- xgboost::xgb.train(data = dtrain, nrounds = 2)
#   result[[i]] <- importance.score(fit = fit.model, Y = as.factor(y),
#                                   X = as.matrix(x))$shap
# }

result1 <- vector(mode = "list", length = length(Z))
for (i in 1:length(Z)) {
  cat('--Calculate',i,'knockoff statistics.\n')
  result[[i]] <- stat.glmnet_coefdiff(X, Z[[i]], y,family = 'binomial')
}
select1 = knockoff.select(Ws1,fdr = 0.1,offset = 1)
select$index
length(select$index)


library(tidyverse)
result2 <- vector(mode = "list", length = length(Z))
for (i in 1:length(Z)) {
  cat('--Calculate',i,'knockoff statistics.\n')
  result2[[i]] <- stat.SHAP(X, Z[[i]], y)
}
select2 = knockoff.select(result2,fdr = 0.1,offset = 1)
select2$index
length(select2$index)

colnames(X)[select2$index]
